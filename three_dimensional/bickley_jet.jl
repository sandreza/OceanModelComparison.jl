include("CNSE.jl")
include("ThreeDimensionalCompressibleNavierStokesEquations.jl")

function Config(
    name,
    resolution,
    domain,
    params;
    numerical_flux_first_order = RusanovNumericalFlux(),
    Nover = 0,
    periodicity = (true, true, true),
    boundary = ((0, 0), (0, 0), (0, 0)),
    boundary_conditons = (),
)
    mpicomm = MPI.COMM_WORLD
    ArrayType = ClimateMachine.array_type()

    xrange =
        range(-domain.Lˣ / 2; length = resolution.Nˣ + 1, stop = domain.Lˣ / 2)
    yrange =
        range(-domain.Lʸ / 2; length = resolution.Nʸ + 1, stop = domain.Lʸ / 2)
    zrange =
        range(-domain.Lᶻ / 2; length = resolution.Nᶻ + 1, stop = domain.Lᶻ / 2)

    brickrange = (xrange, yrange, zrange)

    topl = BrickTopology(
        mpicomm,
        brickrange,
        periodicity = periodicity,
        boundary = boundary,
    )

    grid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = ArrayType,
        polynomialorder = resolution.N + Nover,
    )

    model = ThreeDimensionalCompressibleNavierStokes.CNSE3D{FT}(
        (domain.Lˣ, domain.Lʸ, domain.Lᶻ),
        ClimateMachine.Ocean.NonLinearAdvectionTerm(),
        ThreeDimensionalCompressibleNavierStokes.ConstantViscosity{FT}(
            μ = params.μ,
            ν = params.ν,
            κ = params.κ,
        ),
        nothing,
        nothing,
        boundary_conditons;
        cₛ = params.cₛ,
        ρₒ = params.ρₒ,
    )

    dg = DGModel(
        model,
        grid,
        numerical_flux_first_order,
        CentralNumericalFluxSecondOrder(),
        CentralNumericalFluxGradient(),
    )

    return Config(name, dg, Nover, mpicomm, ArrayType)
end

import ClimateMachine.Ocean: ocean_init_state!, ocean_init_aux!

function ocean_init_state!(
    model::ThreeDimensionalCompressibleNavierStokes.CNSE3D,
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
    z = aux.z

    # The Bickley jet
    #=
    U = cosh(y)^(-2)

    # Slightly off-center vortical perturbations
    Ψ = exp(-(y + l / 10)^2 / (2 * (l^2))) * cos(k * x) * cos(k * y)

    # Vortical velocity fields (ũ, ṽ) = (-∂ʸ, +∂ˣ) ψ̃
    u = Ψ * (k * tan(k * y) + y / (l^2))
    v = -Ψ * k * tan(k * x)

    ρ = model.ρₒ
    state.ρ = ρ
    state.ρu = ρ * @SVector [U + ϵ * u, ϵ * v, -0]
    state.ρθ = ρ * sin(k * y)
    =#
    #=
    U = sech(y)^2
    V = 0
    W = 0
    # Slightly off-center vortical perturbations
    Ψ₁ = exp(-(y + l / 10)^2 / (2 * (l^2))) * cos(k * x) * cos(k * y)
    Ψ₂ = exp(-(z + l / 10)^2 / (2 * (l^2))) * cos(k * x) * cos(k * z)
    # Vortical velocity fields (u, v, w) = (-∂ʸ, +∂ˣ, 0) Ψ₁ + (0, -∂ᶻ, +∂ˣ)Ψ₂
    u =  Ψ₁ * (k * tan(k * y) + y / (l^2) + 1/(10 * l)) 
    v = -Ψ₁ * k * tan(k * x)  + Ψ₂ * (k * tan(k * z) + z / (l^2) + 1/(10 * l)) 
    w = -Ψ₂ * k * tan(k * x) 
    =#
    U = sech(y)^2
    V = 0
    W = 0
    # Slightly off-center vortical perturbations
    Ψ₁ = exp(-(y + l / 10)^2 / (2 * (l^2))) * cos(k * x) * cos(k * y)
    Ψ₂ = exp(-(z + l / 10)^2 / (2 * (l^2))) * cos(k * y) * cos(k * z)
    # Vortical velocity fields (u, v, w) = (-∂ʸ, +∂ˣ, 0) Ψ₁ + (0, -∂ᶻ, +∂ʸ)Ψ₂ 
    u =  Ψ₁ * (k * tan(k * y) + y / (l^2) + 1/(10 * l)) 
    v = -Ψ₁ * k * tan(k * x)  + Ψ₂ * (k * tan(k * z) + z / (l^2) + 1/(10 * l)) 
    w = -Ψ₂ * k * tan(k * y) 

    ρ = model.ρₒ
    state.ρ = ρ
    state.ρu = ρ * @SVector [U + ϵ * u, V + ϵ * v, W + ϵ * w]
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
        aux.z = geom.coord[3]
    end

    return nothing
end