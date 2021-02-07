#!/usr/bin/env julia --project

include("bickley_jet.jl")
include("convenience.jl")

FT = Float64

#################
# RUN THE TESTS #
#################

# Domain size
Lˣ = 4 * FT(π)  # m
Lʸ = 4 * FT(π)  # m
Lᶻ = 4 * FT(π)  # m

vtkpath = abspath(joinpath(ClimateMachine.Settings.output_dir, "vtk_bickley_3D"))

DOFs = [32]
Ns = [1, 2, 3, 4]
Novers = [0, 1]
fluxes = [RoeNumericalFlux(), RusanovNumericalFlux()]
periodicity = [false, true]

#=
DOFs = [32]
Ns = [4]
Novers = [1]
fluxes = [RoeNumericalFlux(), RusanovNumericalFlux()]
periodicity = [true]
=#
for DOF in DOFs, N in Ns, Nover in Novers, flux in fluxes, periodic in periodicity
    # simulation times
    timeend = FT(200) # s

    filename, Ne, dt = generate_name(DOF, N, Nover, flux, periodic)
    println("currently doing " * filename)
    # Domain Resolution
    Nˣ = Ne
    Nʸ = Ne
    Nᶻ = 1

    nout = round(Int, 2 / dt)
    # model params
    cₛ = sqrt(10) # m/s
    ρₒ = 1 # kg/m³
    μ = 0 # 1e-6,   # m²/s
    ν = 0 # 1e-6,   # m²/s
    κ = 0 # 1e-6,   # m²/s

    resolution = (; N, Nˣ, Nʸ, Nᶻ)
    domain = (; Lˣ, Lʸ, Lᶻ)
    timespan = (; dt, nout, timeend)
    params = (; cₛ, ρₒ, μ, ν, κ)
    if periodic
        config = Config(
            filename,
            resolution,
            domain,
            params;
            numerical_flux_first_order = flux,
            Nover = Nover,
            periodicity = (true, true, true),
            boundary = ((0, 0), (0, 0), (0, 0)),
            boundary_conditons = (ClimateMachine.Ocean.OceanBC(
                Impenetrable(FreeSlip()),
                Insulating(),
            ),),
        )
    else
        config = Config(
            filename,
            resolution,
            domain,
            params;
            numerical_flux_first_order = flux,
            Nover = Nover,
            periodicity = (true, periodic, true),
            boundary = ((0, 0), (1, 1), (0, 0)),
            boundary_conditons = (ClimateMachine.Ocean.OceanBC(
                Impenetrable(FreeSlip()),
                Insulating(),
            ),),
        )
    end
    tic = Base.time()

    run_CNSE(config, resolution, timespan; TimeStepper = SSPRK22Heuns)
    # LSRK54CarpenterKennedy
    toc = Base.time()
    time = toc - tic

    f = jldopen(filename * ".jld2", "a+")
    f["simulationtime"] = time
    f["threads"] = Threads.nthreads()
    f["arraytype"] = string(ClimateMachine.array_type())
    println(time, " seconds")
end