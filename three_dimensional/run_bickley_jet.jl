#!/usr/bin/env julia --project

include("bickley_jet.jl")

FT = Float64

#################
# RUN THE TESTS #
#################

vtkpath = abspath(joinpath(ClimateMachine.Settings.output_dir, "vtk_bickley_3D"))

let
    # simulation times
    timeend = FT(200) # s
    dt = FT(0.02) # s
    nout = Int(100)

    # Domain Resolution
    N = 3
    Nˣ = 8
    Nʸ = 8
    Nᶻ = 1

    # Domain size
    Lˣ = 4 * FT(π)  # m
    Lʸ = 4 * FT(π)  # m
    Lᶻ = 4 * FT(π)  # m

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

    config = Config(
        "rusanov_overintegration",
        resolution,
        domain,
        params;
        numerical_flux_first_order = RusanovNumericalFlux(),
        Nover = 0,
        periodicity = (true, true, true),
        boundary = ((0, 0), (0, 0), (0, 0)),
        boundary_conditons = (ClimateMachine.Ocean.OceanBC(
            Impenetrable(FreeSlip()),
            Insulating(),
        ),),
    )

    tic = Base.time()

    run_CNSE(config, resolution, timespan; TimeStepper = LSRK54CarpenterKennedy)

    toc = Base.time()
    time = toc - tic
    println(time)
end