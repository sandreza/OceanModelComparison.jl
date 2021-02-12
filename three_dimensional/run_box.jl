include("box.jl")
# ClimateMachine.inFT = Float64

#################
# Initial State #
#################
import ClimateMachine.Ocean: ocean_init_state!

function ocean_init_state!(
    model::ThreeDimensionalCompressibleNavierStokes.CNSE3D,
    state,
    aux,
    localgeo,
    t,
)

    x = aux.x
    y = aux.y
    z = aux.z

    ϵ = 0.01
    ρ = model.ρₒ
    state.ρ = ρ
    state.ρu = ρ * @SVector [-ϵ * sin(2π*y/100.0), -0, -0.00]
    state.ρθ = ρ * (20 + z / 100)

    return nothing
end

#################
# RUN THE TESTS #
#################

vtkpath = abspath(joinpath(ClimateMachine.Settings.output_dir, "vtk_box_3D"))

let
    # simulation times
    timeend = FT(200) # s
    dt = FT(10.0) # s
    nout = Int(20)

    # Domain Resolution
    N  = 1
    Nˣ = 2
    Nʸ = 2
    Nᶻ = 2

    # Domain size
    Lˣ = 100.0  # m
    Lʸ = 100.0  # m
    Lᶻ = 100.0  # m

    # model params
    cₛ = 5.0    # m/s
    ρₒ = 1.0    # kg/m³
    μ  = 0.0    # 1e-6,   # m²/s
    ν  = 1e-2   # m²/s
    κ  = 1e-2   # m²/s
    α  = 2e-4   # 1/K
    g  = 10.0   # m/s²

    resolution = (; N, Nˣ, Nʸ, Nᶻ)
    domain = (; Lˣ, Lʸ, Lᶻ)
    timespan = (; dt, nout, timeend)
    params = (; cₛ, ρₒ, μ, ν, κ, α, g)

    BC = (
        ClimateMachine.Ocean.OceanBC(
            Impenetrable(NoSlip()), 
            TemperatureFlux((state, aux, t) -> (0.0))
        ),
        ClimateMachine.Ocean.OceanBC(
            Impenetrable(KinematicStress(
                (state, aux, t) -> (@SVector [-0, -0, -0]),
            )),
            TemperatureFlux((state, aux, t) -> (5.0e-6)),
        ),
    )

    config = Config(
        "cool_the_box",
        resolution,
        domain,
        params;
        numerical_flux_first_order = RoeNumericalFlux(),
        Nover = 1,
        periodicity = (true, true, false),
        boundary = ((0, 0), (0, 0), (1, 2)),
        boundary_conditons = BC,
    )

    tic = Base.time()

    run_CNSE(config, resolution, timespan; TimeStepper = SSPRK22Heuns)

    toc = Base.time()
    time = toc - tic
    println(time)
end