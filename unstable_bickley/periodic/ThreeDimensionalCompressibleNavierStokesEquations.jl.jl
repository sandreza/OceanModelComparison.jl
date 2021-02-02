module ThreeDimensionalCompressibleNavierStokes

export ThreeDimensionalCompressibleNavierStokesEquations

using Test
using StaticArrays
using LinearAlgebra: dot, Diagonal, eigen

using ClimateMachine.Ocean
using ClimateMachine.VariableTemplates
using ClimateMachine.Mesh.Geometry
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.BalanceLaws
using ClimateMachine.Ocean: kinematic_stress, coriolis_parameter
using ClimateMachine.Mesh.Geometry: LocalGeometry
using ClimateMachine.MPIStateArrays: MPIStateArray

import ClimateMachine.BalanceLaws:
    vars_state,
    init_state_prognostic!,
    init_state_auxiliary!,
    compute_gradient_argument!,
    compute_gradient_flux!,
    flux_first_order!,
    flux_second_order!,
    source!,
    wavespeed,
    boundary_conditions,
    boundary_state!
import ClimateMachine.Ocean:
    ocean_init_state!,
    ocean_init_aux!,
    ocean_boundary_state!,
    _ocean_boundary_state!
import ClimateMachine.NumericalFluxes: numerical_flux_first_order!

×(a::SVector, b::SVector) = StaticArrays.cross(a, b)
⋅(a::SVector, b::SVector) = StaticArrays.dot(a, b)
⊗(a::SVector, b::SVector) = a * b'

abstract type TurbulenceClosure end
struct LinearDrag{T} <: TurbulenceClosure
    λ::T
end
struct ConstantViscosity{T} <: TurbulenceClosure
    ν::T
    κ::T
    function ConstantViscosity{T}(;
        ν = FT(1e-6),   # m²/s
        κ = FT(1e-6),   # m²/s
    ) where {T <: AbstractFloat}
        return new{T}(ν, κ)
    end
end

abstract type CoriolisForce end
struct fPlaneCoriolis{T} <: CoriolisForce
    fₒ::T
    β::T
    function fPlaneCoriolis{T}(;
        fₒ = T(1e-4), # Hz
        β = T(1e-11), # Hz/m
    ) where {T <: AbstractFloat}
        return new{T}(fₒ, β)
    end
end

abstract type Forcing end
struct KinematicStress{T} <: Forcing
    τₒ::T
    function KinematicStress{T}(; τₒ = T(1e-4)) where {T <: AbstractFloat}
        return new{T}(τₒ)
    end
end

"""
    ThreeDimensionalCompressibleNavierStokesEquations <: BalanceLaw
A `BalanceLaw` for shallow water modeling.
write out the equations here
# Usage
    ThreeDimensionalCompressibleNavierStokesEquations()
"""
struct ThreeDimensionalCompressibleNavierStokesEquations{D, A, T, C, F, BC, FT} <:
       BalanceLaw
    domain::D
    advection::A
    turbulence::T
    coriolis::C
    forcing::F
    boundary_conditions::BC
    g::FT
    c::FT
    function ThreeDimensionalCompressibleNavierStokesEquations{FT}(
        domain::D,
        advection::A,
        turbulence::T,
        coriolis::C,
        forcing::F,
        boundary_conditions::BC;
        g = FT(10), # m/s²
        c = FT(0),  #m/s
    ) where {FT <: AbstractFloat, D, A, T, C, F, BC}
        return new{D, A, T, C, F, BC, FT}(
            domain,
            advection,
            turbulence,
            coriolis,
            forcing,
            boundary_conditions,
            g,
            c,
        )
    end
end
CNSE3D = ThreeDimensionalCompressibleNavierStokesEquations

function vars_state(m::CNSE3D, ::Prognostic, T)
    @vars begin
        ρ::T
        ρu::SVector{3, T}
        ρθ::T
    end
end

function init_state_prognostic!(m::CNSE3D, state::Vars, aux::Vars, localgeo, t)
    ocean_init_state!(m, state, aux, localgeo, t)
end

function vars_state(m::CNSE3D, ::Auxiliary, T)
    @vars begin
        x::T
        y::T
    end
end

function init_state_auxiliary!(
    model::CNSE3D,
    state_auxiliary::MPIStateArray,
    grid,
    direction,
)
    init_state_auxiliary!(
        model,
        (model, aux, tmp, geom) -> ocean_init_aux!(model, aux, geom),
        state_auxiliary,
        grid,
        direction,
    )
end

function vars_state(m::CNSE3D, ::Gradient, T)
    @vars begin
        ∇u::SVector{3, T}
        ∇θ::T
    end
end

function compute_gradient_argument!(
    model::CNSE3D,
    grad::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    compute_gradient_argument!(model.turbulence, grad, state, aux, t)
end

compute_gradient_argument!(::LinearDrag, _...) = nothing

@inline function compute_gradient_argument!(
    ::ConstantViscosity,
    grad::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    ρ = state.ρ
    ρu = state.ρu
    ρθ = state.ρθ

    u = ρu / ρ
    θ = ρθ / ρ

    grad.∇u = u
    grad.∇θ = θ

    return nothing
end

function vars_state(m::CNSE3D, ::GradientFlux, T)
    @vars begin
        ν∇u::SMatrix{3, 3, T, 9}
        κ∇θ::SVector{3, T}
    end
end

function compute_gradient_flux!(
    model::CNSE3D,
    gradflux::Vars,
    grad::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    compute_gradient_flux!(
        model,
        model.turbulence,
        gradflux,
        grad,
        state,
        aux,
        t,
    )
end

compute_gradient_flux!(::CNSE3D, ::LinearDrag, _...) = nothing

@inline function compute_gradient_flux!(
    ::CNSE3D,
    turb::ConstantViscosity,
    gradflux::Vars,
    grad::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    ν = Diagonal(@SVector [turb.ν, turb.ν, -0])
    κ = Diagonal(@SVector [turb.κ, turb.κ, -0])

    gradflux.ν∇u = -ν * grad.∇u
    gradflux.κ∇θ = -κ * grad.∇θ

    return nothing
end

@inline function flux_first_order!(
    model::CNSE3D,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    ρ = state.ρ
    ρu = @SVector [state.ρu[1], state.ρu[2], state.ρu[3]]
    ρθ = state.ρθ
    g = model.g

    I = @SMatrix [
        1 -0 -0
        -0 1 -0
        -0 -0 1
    ]

    flux.ρ += ρu
    flux.ρu += g * ρ^2 * I / 2

    advective_flux!(model, model.advection, flux, state, aux, t)

    return nothing
end

advective_flux!(::CNSE3D, ::Nothing, _...) = nothing

@inline function advective_flux!(
    ::CNSE3D,
    ::NonLinearAdvectionTerm,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    ρ = state.ρ
    ρu = state.ρu
    ρv = @SVector [state.ρu[1], state.ρu[2], state.ρu[3]]
    ρθ = state.ρθ

    flux.ρu += ρv ⊗ ρu / ρ
    flux.ρθ += ρv * ρθ / ρ

    return nothing
end

function flux_second_order!(
    model::CNSE3D,
    flux::Grad,
    state::Vars,
    gradflux::Vars,
    ::Vars,
    aux::Vars,
    t::Real,
)
    flux_second_order!(model, model.turbulence, flux, state, gradflux, aux, t)
end

flux_second_order!(::CNSE3D, ::LinearDrag, _...) = nothing

@inline function flux_second_order!(
    ::CNSE3D,
    ::ConstantViscosity,
    flux::Grad,
    state::Vars,
    gradflux::Vars,
    aux::Vars,
    t::Real,
)
    flux.ρu += gradflux.ν∇u
    flux.ρθ += gradflux.κ∇θ

    return nothing
end

@inline function source!(
    model::CNSE3D,
    source::Vars,
    state::Vars,
    gradflux::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    coriolis_force!(model, model.coriolis, source, state, aux, t)
    forcing_term!(model, model.forcing, source, state, aux, t)
    linear_drag!(model, model.turbulence, source, state, aux, t)

    return nothing
end

coriolis_force!(::CNSE3D, ::Nothing, _...) = nothing

@inline function coriolis_force!(
    model::CNSE3D,
    coriolis::fPlaneCoriolis,
    source,
    state,
    aux,
    t,
)
    ρu = @SVector [state.ρu[1], state.ρu[2], state.ρu[3]]

    # f × u
    f = [-0, -0, coriolis_parameter(model, coriolis, aux.coords)]

    source.ρu -= f × ρu

    return nothing
end

forcing_term!(::CNSE3D, ::Nothing, _...) = nothing

@inline function forcing_term!(
    model::CNSE3D,
    forcing::KinematicStress,
    source,
    state,
    aux,
    t,
)
    source.ρu += kinematic_stress(model, forcing, aux.coords)

    return nothing
end

linear_drag!(::CNSE3D, ::ConstantViscosity, _...) = nothing

@inline function linear_drag!(::CNSE3D, turb::LinearDrag, source, state, aux, t)
    source.ρu -= turb.λ * state.ρu

    return nothing
end

@inline wavespeed(m::CNSE3D, _...) = m.c

roe_average(ρ⁻, ρ⁺, var⁻, var⁺) =
    (sqrt(ρ⁻) * var⁻ + sqrt(ρ⁺) * var⁺) / (sqrt(ρ⁻) + sqrt(ρ⁺))

function numerical_flux_first_order!(
    ::RoeNumericalFlux,
    model::CNSE3D,
    fluxᵀn::Vars{S},
    n⁻::SVector,
    state⁻::Vars{S},
    aux⁻::Vars{A},
    state⁺::Vars{S},
    aux⁺::Vars{A},
    t,
    direction,
) where {S, A}
    numerical_flux_first_order!(
        CentralNumericalFluxFirstOrder(),
        model,
        fluxᵀn,
        n⁻,
        state⁻,
        aux⁻,
        state⁺,
        aux⁺,
        t,
        direction,
    )

    FT = eltype(fluxᵀn)

    # constants and normal vectors
    g = model.g
    @inbounds nˣ = n⁻[1]
    @inbounds nʸ = n⁻[2]
    @inbounds nᶻ = n⁻[3]

    # get minus side states
    ρ⁻ = state⁻.ρ
    @inbounds ρu⁻ = state⁻.ρu[1]
    @inbounds ρv⁻ = state⁻.ρu[2]
    @inbounds ρw⁻ = state⁻.ρu[3]
    ρθ⁻ = state⁻.ρθ

    u⁻ = ρu⁻ / ρ⁻
    v⁻ = ρv⁻ / ρ⁻
    w⁻ = ρw⁻ / ρ⁻
    θ⁻ = ρθ⁻ / ρ⁻

    # get plus side states
    ρ⁺ = state⁺.ρ
    @inbounds ρu⁺ = state⁺.ρu[1]
    @inbounds ρv⁺ = state⁺.ρu[2]
    @inbounds ρw⁺ = state⁺.ρu[3]
    ρθ⁺ = state⁺.ρθ

    u⁺ = ρu⁺ / ρ⁺
    v⁺ = ρv⁺ / ρ⁺
    w⁺ = ρw⁺ / ρ⁺
    θ⁺ = ρθ⁺ / ρ⁺

    # averages for roe fluxes
    ρ  = (ρ⁺  + ρ⁻) / 2
    ρu = (ρu⁺ + ρu⁻) / 2
    ρv = (ρv⁺ + ρv⁻) / 2
    ρw = (ρw⁺ + ρw⁻) / 2
    ρθ = (ρθ⁺ + ρθ⁻) / 2

    u = roe_average(ρ⁻, ρ⁺, u⁻, u⁺)
    v = roe_average(ρ⁻, ρ⁺, v⁻, v⁺)
    w = roe_average(ρ⁻, ρ⁺, w⁻, w⁺)
    θ = roe_average(ρ⁻, ρ⁺, θ⁻, θ⁺)

    # normal and tangent velocities
    uₙ = nˣ * u + nʸ * v + nᶻ * v
    # uₚ = nˣ * v - nʸ * u

    # differences for difference vector
    Δρ  = ρ⁺ - ρ⁻
    Δρu = ρu⁺ - ρu⁻
    Δρv = ρv⁺ - ρv⁻
    Δρw = ρw⁺ - ρw⁻
    Δρθ = ρθ⁺ - ρθ⁻

    Δφ = @SVector [Δρ, Δρu, Δρv, Δρw, Δρθ]

    """
    # jacobian
    ∂F∂φ = [
        0 nˣ nʸ 0
        (nˣ * c^2 - u * uₙ) (uₙ + nˣ * u) (nʸ * u) 0
        (nʸ * c^2 - v * uₙ) (nˣ * v) (uₙ + nʸ * v) 0
        (-θ * uₙ) (nˣ * θ) (nʸ * θ) uₙ
    ]
    # eigen decomposition
    λ, R = eigen(∂F∂φ)
    """

    # eigen values matrix
    c = sqrt(g * ρ)
    λ = @SVector [abs(uₙ), abs(uₙ), abs(uₙ + c), abs(uₙ - c), abs(uₙ)]
    # Λ = Diagonal(abs.(λ))


    # eigenvector matrix (doesn't work in 3D)
    R = @SMatrix [
     0    0       1            1           0
     nʸ   0   (u + nˣ * c)  (u - nˣ * c)   0
    -nˣ  -nᶻ  (v + nʸ * c)  (v - nʸ * c)   0
     0    nʸ  (w + nᶻ * c)  (w - nᶻ * c)   0
     0    0       θ            θ           1
    ]

    # actually calculate flux
    parent(fluxᵀn) .-= R * (λ .* (R \ Δφ)) * 0.5
    # parent(fluxᵀn) .-= R * Λ * R⁻¹ * Δφ / 2

    return nothing
end

boundary_conditions(model::CNSE3D) = model.boundary_conditions

"""
    boundary_state!(nf, ::CNSE3D, args...)
applies boundary conditions for the hyperbolic fluxes
dispatches to a function in OceanBoundaryConditions
"""
@inline function boundary_state!(nf, bc, model::CNSE3D, args...)
    return _ocean_boundary_state!(nf, bc, model, args...)
end

"""
    ocean_boundary_state!(nf, bc::OceanBC, ::CNSE3D)
splits boundary condition application into velocity
"""
@inline function ocean_boundary_state!(nf, bc::OceanBC, m::CNSE3D, args...)
    return ocean_boundary_state!(nf, bc.velocity, m, m.turbulence, args...)
    return ocean_boundary_state!(nf, bc.temperature, m, args...)
end

include("bc_velocity.jl")

end