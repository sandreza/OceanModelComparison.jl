using Test
using ClimateMachine
ClimateMachine.init()
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.Mesh.Filters
using ClimateMachine.VariableTemplates
using ClimateMachine.Ocean
using ClimateMachine.Ocean.Domains
using ClimateMachine.Ocean.Fields
using ClimateMachine.Mesh.Grids: DiscontinuousSpectralElementGrid
using ClimateMachine.GenericCallbacks: EveryXSimulationTime
using ClimateMachine.GenericCallbacks: EveryXSimulationSteps
using ClimateMachine.Ocean: current_step, Δt, current_time
using ClimateMachine.Ocean: JLD2Writer, OutputTimeSeries, write!
using CLIMAParameters: AbstractEarthParameterSet, Planet

polynomialorders(::DiscontinuousSpectralElementGrid{T, dim, N}) where {T, dim, N} = Tuple([N for i in 1:dim])
using ClimateMachine.Ocean

include(pwd() * "/unstable_bickley/periodic/TwoDimensionalCompressibleNavierStokesEquations.jl")

using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.BalanceLaws: vars_state, Prognostic, Auxiliary
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.MPIStateArrays
using ClimateMachine.VTK

using MPI
using LinearAlgebra
using StaticArrays
using Logging, Printf, Dates

import ClimateMachine.Ocean: ocean_init_state!, ocean_init_aux!

function ocean_init_state!(
    ::TwoDimensionalCompressibleNavierStokes.TwoDimensionalCompressibleNavierStokesEquations,
    state,
    aux,
    localgeo,
    t,
)
    ϵ = 0.1 # perturbation magnitude
    l = 0.5 # Gaussian width
    k = 0.5 # Sinusoidal wavenumber

    x = aux.x
    y = aux.y

    # The Bickley jet
    U = cosh(y)^(-2)

    # Slightly off-center vortical perturbations
    Ψ = exp(-(y + l / 10)^2 / (2 * (l^2))) * cos(k * x) * cos(k * y)

    # Vortical velocity fields (ũ, ṽ) = (-∂ʸ, +∂ˣ) ψ̃
    u = Ψ * (k * tan(k * y) + y / (l^2))
    v = -Ψ * k * tan(k * x)

    ρ = 1
    state.ρ = ρ
    state.ρu = ρ * @SVector [U + ϵ * u, ϵ * v]
    state.ρθ = ρ * sin(k * y)

    return nothing
end

function ocean_init_aux!(
    ::TwoDimensionalCompressibleNavierStokes.CNSE2D,
    aux,
    geom,
)
    @inbounds begin
        aux.x = geom.coord[1]
        aux.y = geom.coord[2]
    end

    return nothing
end

function run_bickley_jet(params)
    mpicomm = MPI.COMM_WORLD
    ArrayType = ClimateMachine.array_type()

    xrange = range(-params.Lˣ / 2; length = params.Nˣ + 1, stop = params.Lˣ / 2)
    yrange = range(-params.Lʸ / 2; length = params.Nʸ + 1, stop = params.Lʸ / 2)

    brickrange = (xrange, yrange)
    topl = BrickTopology(
        mpicomm,
        brickrange,
        periodicity = (true, true),
        boundary = ((0, 0), (0, 0)),
    )
    grid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = ArrayType,
        polynomialorder = params.N,
    )

    model = TwoDimensionalCompressibleNavierStokes.CNSE2D{FT}(
        (params.Lˣ, params.Lʸ),
        ClimateMachine.Ocean.NonLinearAdvectionTerm(),
        TwoDimensionalCompressibleNavierStokes.ConstantViscosity{FT}(
            ν = 0, # 1e-6,   # m²/s
            κ = 0, # 1e-6,   # m²/s
        ),
        nothing,
        nothing,
        ClimateMachine.Ocean.OceanBC(Impenetrable(FreeSlip()), Insulating());
        g = 10, # m/s²
        c = 2, # m/s
    )

    dg = DGModel(
        model,
        grid,
        RusanovNumericalFlux(),
        CentralNumericalFluxSecondOrder(),
        CentralNumericalFluxGradient(),
    )

    Q = init_ode_state(dg, FT(0); init_on_cpu = true)

    lsrk = LSRK54CarpenterKennedy(dg, Q, dt = params.dt, t0 = 0)

    odesolver = lsrk

    vtkstep = 0
    cbvector = make_callbacks(
        vtkpath,
        vtkstep,
        params,
        mpicomm,
        odesolver,
        dg,
        model,
        Q,
    )

    eng0 = norm(Q)
    @info @sprintf """Starting
    norm(Q₀) = %.16e
    ArrayType = %s""" eng0 ArrayType

    solve!(Q, odesolver; timeend = params.timeend, callbacks = cbvector)

    return nothing
end

function make_callbacks(
    vtkpath,
    vtkstep,
    params,
    mpicomm,
    odesolver,
    dg,
    model,
    Q,
)
    if isdir(vtkpath)
        rm(vtkpath, recursive = true)
    end
    mkpath(vtkpath)

    function do_output(vtkstep, model, dg, Q)
        outprefix = @sprintf(
            "%s/mpirank%04d_step%04d",
            vtkpath,
            MPI.Comm_rank(mpicomm),
            vtkstep
        )
        @info "doing VTK output" outprefix
        statenames = flattenednames(vars_state(model, Prognostic(), eltype(Q)))
        auxnames = flattenednames(vars_state(model, Auxiliary(), eltype(Q)))
        writevtk(outprefix, Q, dg, statenames, dg.state_auxiliary, auxnames)

        vtkstep += 1

        return vtkstep
    end

    vtkstep = do_output(vtkstep, model, dg, Q)
    cbvtk =
        GenericCallbacks.EveryXSimulationSteps(params.nout) do (init = false)
            vtkstep = do_output(vtkstep, model, dg, Q)
            return nothing
        end

    starttime = Ref(now())
    cbinfo = GenericCallbacks.EveryXWallTimeSeconds(60, mpicomm) do (s = false)
        if s
            starttime[] = now()
        else
            energy = norm(Q)
            @info @sprintf(
                """Update
                simtime = %8.2f / %8.2f
                runtime = %s
                norm(Q) = %.16e""",
                ODESolvers.gettime(odesolver),
                params.timeend,
                Dates.format(
                    convert(Dates.DateTime, Dates.now() - starttime[]),
                    Dates.dateformat"HH:MM:SS",
                ),
                energy
            )
        end
    end

    return (cbinfo, cbvtk)
end

#################
# RUN THE TESTS #
#################
FT = Float64
vtkpath =
    abspath(joinpath(ClimateMachine.Settings.output_dir, "vtk_bickley_jet"))
tic = time()
let
    timeend = FT(200) # s
    nout = FT(100)
    dt = FT(0.02) # s

    N = 2
    Nˣ = 43
    Nʸ = 43
    Lˣ = 4 * FT(π)  # m
    Lʸ = 4 * FT(π)  # m

    params = (; N, Nˣ, Nʸ, Lˣ, Lʸ, dt, nout, timeend)

    run_bickley_jet(params)
end
toc = time()
println("The amount of time for the simulation is ", toc - tic)
