using ..Ocean: surface_flux

"""
    ocean_boundary_state!(::Union{NumericalFluxFirstOrder, NumericalFluxGradient}, ::Insulating, ::HBModel)
apply insulating boundary condition for temperature
sets transmissive ghost point
"""
function ocean_boundary_state!(
    ::Union{NumericalFluxFirstOrder, NumericalFluxGradient},
    ::Insulating,
    ::CNSE2D,
    state⁺,
    aux⁺,
    n⁻,
    state⁻,
    aux⁻,
    t,
)
    state⁺.ρθ = state⁻.ρθ

    return nothing
end

"""
    ocean_boundary_state!(::NumericalFluxSecondOrder, ::Insulating, ::HBModel)
apply insulating boundary condition for velocity
sets ghost point to have no numerical flux on the boundary for κ∇θ
"""
@inline function ocean_boundary_state!(
    ::NumericalFluxSecondOrder,
    ::Insulating,
    ::CNSE2D,
    state⁺,
    gradflux⁺,
    aux⁺,
    n⁻,
    state⁻,
    gradflux⁻,
    aux⁻,
    t,
)
    state⁺.ρθ = state⁻.ρθ
    gradflux⁺.κ∇θ = n⁻ * -0

    return nothing
end

"""
    ocean_boundary_state!(::Union{NumericalFluxFirstOrder, NumericalFluxGradient}, ::TemperatureFlux, ::HBModel)
apply temperature flux boundary condition for velocity
applies insulating conditions for first-order and gradient fluxes
"""
function ocean_boundary_state!(
    nf::Union{NumericalFluxFirstOrder, NumericalFluxGradient},
    ::TemperatureFlux,
    model::CNSE2D,
    args...,
)
    return ocean_boundary_state!(nf, Insulating(), model, args...)
end

"""
    ocean_boundary_state!(::NumericalFluxSecondOrder, ::TemperatureFlux, ::HBModel)
apply insulating boundary condition for velocity
sets ghost point to have specified flux on the boundary for κ∇θ
"""
@inline function ocean_boundary_state!(
    ::NumericalFluxSecondOrder,
    ::TemperatureFlux,
    model::CNSE2D,
    state⁺,
    gradflux⁺,
    aux⁺,
    n⁻,
    state⁻,
    gradflux⁻,
    aux⁻,
    t,
)
    state⁺.ρθ = state⁻.ρθ
    gradflux⁺.κ∇θ = n⁻ * surface_flux(model, aux⁻.coords, state⁻.ρθ)

    return nothing
end
