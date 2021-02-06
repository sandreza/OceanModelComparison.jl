using Test, JLD2
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
using ClimateMachine.BalanceLaws:
    vars_state, Prognostic, Auxiliary, number_states

polynomialorders(::DiscontinuousSpectralElementGrid{T, dim, N}) where {T, dim, N} = N
using ClimateMachine.Ocean

include(pwd() * "/unstable_bickley/periodic/ThreeDimensionalCompressibleNavierStokesEquations.jl")

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
##
function ocean_init_state!(
    ::ThreeDimensionalCompressibleNavierStokes.ThreeDimensionalCompressibleNavierStokesEquations,
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
    state.ρu = ρ * @SVector [U + ϵ * u, ϵ * v, 0]
    state.ρθ = ρ * sin(k * y)

    return nothing
end

function ocean_init_aux!(
    ::ThreeDimensionalCompressibleNavierStokes.CNSE3D,
    aux,
    geom,
)
    @inbounds begin
        aux.x = geom.coord[1]
        aux.y = geom.coord[2]
    end

    return nothing
end

function run_bickley_jet(params; filename = "example")
    mpicomm = MPI.COMM_WORLD
    ArrayType = ClimateMachine.array_type()

    xrange = range(-params.Lˣ / 2; length = params.Nˣ + 1, stop = params.Lˣ / 2)
    yrange = range(-params.Lʸ / 2; length = params.Nʸ + 1, stop = params.Lʸ / 2)
    zrange = range(-params.Lᶻ / 2; length = params.Nᶻ + 1, stop = params.Lᶻ / 2)

    brickrange = (xrange, yrange, zrange)
    topl = StackedBrickTopology(
        mpicomm,
        brickrange,
        periodicity = (true, true, true),
        boundary = ((0, 0), (0, 0), (0, 0)),
    )
    grid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = ArrayType,
        polynomialorder = params.Nint,
    )
    f = jldopen(filename * ".jld2", "w")
    f["grid"] = grid

    model = ThreeDimensionalCompressibleNavierStokes.CNSE3D{FT}(
        (params.Lˣ, params.Lʸ, params.Lᶻ),
        ClimateMachine.Ocean.NonLinearAdvectionTerm(),
        ThreeDimensionalCompressibleNavierStokes.ConstantViscosity{FT}(
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
        # RoeNumericalFlux(),
        CentralNumericalFluxSecondOrder(),
        CentralNumericalFluxGradient(),
    )

    Q = init_ode_state(dg, FT(0); init_on_cpu = true)
    # LSRK144NiegemannDiehlBusch
    # SSPRK22Heuns
    # LSRK54CarpenterKennedy
    if params.Nint > params.N
        cutoff = CutoffFilter(grid, params.N + 1)
        num_state_prognostic = number_states(model, Prognostic())
        Filters.apply!(Q, 1:num_state_prognostic, grid, cutoff)
    end
    function custom_tendency(tendency, x...; kw...)
        dg(tendency, x...; kw...)
        if params.Nint > params.N
            Filters.apply!(tendency, 1:num_state_prognostic, grid, cutoff)
        end
    end

    lsrk = SSPRK22Heuns(custom_tendency, Q, dt = params.dt, t0 = 0)

    # lsrk = SSPRK22Heuns(dg, Q, dt = params.dt, t0 = 0)

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
        filename = filename
    )

    eng0 = norm(Q)
    @info @sprintf """Starting
    norm(Q₀) = %.16e
    ArrayType = %s""" eng0 ArrayType

    solve!(Q, odesolver; timeend = params.timeend, callbacks = cbvector)
    close(f)
    return dg
end

function make_callbacks(
    vtkpath,
    vtkstep,
    params,
    mpicomm,
    odesolver,
    dg,
    model,
    Q;
    filename = " "
)
    if isdir(vtkpath)
        rm(vtkpath, recursive = true)
    end
    mkpath(vtkpath)

    function do_output(vtkstep, model, dg, Q)
        @info "doing JLD2 output" vtkstep
        f = jldopen(filename * ".jld2", "a+")
        f[string(vtkstep)] = Q.realdata
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


## Run
FT = Float64
vtkpath = abspath(joinpath(ClimateMachine.Settings.output_dir, "vtk_bickley_jet"))
effective_node_spacing(Ne, Np, Lx=4π) = Lx / (Ne * (Np + 1)^2)
for N in [3]
    for DOF in [32]
# N = 4
Nint = N + 0
# DOF = 128
Ne = round(Int, DOF / (N+1))
Nˣ = Ne
Nʸ = Ne
Nᶻ = 1
Lˣ = 4 * FT(π)  # m
Lʸ = 4 * FT(π)  # m
# Lᶻ = 4 * FT(π)  # m
Lᶻ = 4 * FT(π)
# grid
xrange = range(-Lˣ / 2; length = Nˣ + 1, stop = Lˣ / 2)
yrange = range(-Lʸ / 2; length = Nʸ + 1, stop = Lʸ / 2)
zrange = range(-Lᶻ / 2; length = Nᶻ + 1, stop = Lᶻ / 2)
mpicomm = MPI.COMM_WORLD
brickrange = (xrange, yrange, zrange)
topl = StackedBrickTopology(
    mpicomm,
    brickrange,
    periodicity = (true, true, true),
    boundary = ((0, 0), (0, 0), (0,0)),
)
grid = DiscontinuousSpectralElementGrid(
    topl,
    FloatType = FT,
    DeviceArray = Array,
    polynomialorder = Nint,
)
Δx =  min_node_distance(grid)
cfl = 0.3 /2
dt = cfl * Δx / √10
# run
timeend = FT(200) # s
nout = round(Int, 2 / dt)
dt = 2 / nout

params = (; N, Nˣ, Nʸ, Nᶻ, Lˣ, Lʸ, Lᶻ, dt, nout, timeend, Nint)

filename = "3D_overint_p" * string(N) * "_N" * string(Ne)
# filename = "3D_compare_p" * string(N) * "_N" * string(Ne)
# filename = "3D_roe_overint_p" * string(N) * "_N" * string(Ne)
# filename = "3D_roe_p" * string(N) * "_N" * string(Ne) * "_Nint" * string(Nint)
# filename = "deletemeagain"
# filename = "compare_p" * string(N) * "_N" * string(Ne)
tic = time()
dgmodel = run_bickley_jet(params, filename = filename)
toc = time()
println("The amount of time for the simulation is ", toc - tic)
f = jldopen(filename * ".jld2", "a+")
f["simulationtime"] = toc - tic
f["threads"] = Threads.nthreads()
f["iop"] = Nint
close(f)
    end
end

##
N = 3
DOF = 32
Nint = N+0
Ne = round(Int, DOF / (N+1))
filename = "3D_overint_p" * string(N) * "_N" * string(Ne)
#filename = "3D_compare_p" * string(N) * "_N" * string(Ne)
# filename = "3D_roe_p" * string(N) * "_N" * string(Ne) * "_Nint" * string(Nint)
# filename = "overint_p" * string(N) * "_N" * string(Ne)
f = jldopen(filename * ".jld2", "r+")
include(pwd() * "/unstable_bickley/periodic/imperohooks.jl")
include(pwd() * "/unstable_bickley/periodic/vizinanigans2.jl")

dg_grid = f["grid"]
gridhelper = GridHelper(dg_grid)   
x, y, z = coordinates(dg_grid)
xC, yC, zC = cellcenters(dg_grid)
ϕ =  ScalarField(copy(x), gridhelper)

newx = range(-2π, 2π, length = DOF * 2 )
newy = range(-2π, 2π, length = DOF * 2 )
newz = range(-2π, 2π, length = 6 )
##
ρ  = zeros(length(newx), length(newy), length(newz), 101)
ρu = zeros(length(newx), length(newy), length(newz), 101)
ρv = zeros(length(newx), length(newy), length(newz), 101)
ρw = zeros(length(newx), length(newy), length(newz), 101)
ρθ = zeros(length(newx), length(newy), length(newz), 101)
tic = time()
for i in 0:100
    Q = f[string(i)]
    ϕ .= Q[:,1,:]
    ρ[:,:,:, i+1]  = ϕ(newx, newy, newz)
    ϕ .= Q[:,2,:]
    ρu[:,:,:, i+1] = ϕ(newx, newy, newz)
    ϕ .= Q[:,3,:]
    ρv[:,:,:, i+1] = ϕ(newx, newy, newz)
    ϕ .= Q[:,4,:]
    ρw[:,:,:, i+1] = ϕ(newx, newy, newz)
    ϕ .= Q[:,5,:]
    ρθ[:,:,:, i+1] = ϕ(newx, newy, newz)
end
toc = time()
close(f)
println("time to interpolate is $(toc-tic)")
##
ind = 100
states = [ρ[:,:,:,ind], ρu[:,:,:,ind], ρv[:,:,:,ind], ρw[:,:,:,ind], ρθ[:,:,:,ind]]

# states = [ρ[:,:,1,:], ρu[:,:,1,:], ρv[:,:,1,:], ρw[:,:,1,:], ρθ[:,:,1,:]]
statenames = ["ρ", "ρu", "ρv", "ρw", "ρθ"]
scene = volumeslice(states, statenames = statenames)

##
