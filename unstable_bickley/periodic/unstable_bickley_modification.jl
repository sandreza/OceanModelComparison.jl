using Printf
using Plots
using Revise
using CUDA
using ClimateMachine

ClimateMachine.init()

using ClimateMachine.Ocean
using ClimateMachine.Ocean.Domains
using ClimateMachine.Ocean.Fields

using ClimateMachine.Mesh.Grids: DiscontinuousSpectralElementGrid
using ClimateMachine.GenericCallbacks: EveryXSimulationTime
using ClimateMachine.GenericCallbacks: EveryXSimulationSteps
using ClimateMachine.Ocean: current_step, Δt, current_time
using ClimateMachine.Ocean: JLD2Writer, OutputTimeSeries, write!
using CLIMAParameters: AbstractEarthParameterSet, Planet

struct NonDimensionalParameters <: AbstractEarthParameterSet end
Planet.grav(::NonDimensionalParameters) = 10
c = sqrt(Planet.grav(NonDimensionalParameters())) # gravity wave speed for unit depth

include("Bickley.jl")
using .Bickley

function cartesianify!(g)
    offdiagonals = [2,3,4,6,7,8]
    for i in offdiagonals
        g.vgeo[:, i, :] .= -0.0
    end
    return nothing
end

function exactz!(model)
    np, ne = size(model.fields.u.data)
    p = round(Int, np^(1/3))
    for field in eachindex(model.fields)
        newf = reshape(model.fields[field].data, (p,p,p,ne))
        for i in 1:ne
            newf[:,:,:,i] .= newf[:,:,1,i]
        end
    end
    for j in [1, 5, 9]
        newf = reshape(model.grid.vgeo[:,j,:], (p,p,p,ne))
        newf[:,:,:,:] .= newf[1,1,1,1]
        #=
        for i in 1:ne
            newf[:,:,:,i] .= newf[:,:,1,i]
        end
        =#
    end
    return nothing
end

# Low-p assumption:
effective_node_spacing(Ne, Np, Lx=4π) = Lx / (Ne * (Np + 1)^2)

function run(;
             Ne = 4,
             Np = 4,
             ν = 0,
             time_step = 0.1 * effective_node_spacing(Ne, Np) / c,
             array_type = Array,
             output_time_interval = 2,
             stabilizing_dissipation = nothing,
             stop_time = 200,
             augment_name = "")

    name = @sprintf("climate_machine_unstable_bickley_jet_Ne%d_Np%d_ν%.1e_no_rotation", Ne, Np, ν)
    name *= augment_name
    ClimateMachine.Settings.array_type = array_type

    # Domain

    domain = RectangularDomain(Ne = (Ne, Ne, 1), Np = Np,
                               x = (-2π, 2π), y = (-2π, 2π), z = (0, 1),
                               periodicity = (true, true, false))

    # Physical parameters:
    g = Planet.grav(NonDimensionalParameters())

    # Non-dimensional parameters
    ϵ = 0.1 # Perturbation amplitude
    ℓ = 0.5 # Perturbation width
    k = 0.5 # Perturbation wavenumber

    # Initial conditions: Jet/tracer + perturbations
    uᵢ(x, y, z) = Bickley.U(y) + ϵ * Bickley.ũ(x, y, ℓ, k)
    vᵢ(x, y, z) = ϵ * Bickley.ṽ(x, y, ℓ, k)
    θᵢ(x, y, z) = Bickley.C(y, domain.L.y)

    initial_conditions = InitialConditions(u=uᵢ, v=vᵢ, θ=θᵢ)

    model = Ocean.HydrostaticBoussinesqSuperModel(
        domain = domain,
        time_step = time_step,
        initial_conditions = initial_conditions,
        parameters = NonDimensionalParameters(),
        turbulence_closure = (νʰ = ν, κʰ = ν, νᶻ = ν, κᶻ = ν),
        rusanov_wave_speeds = (cʰ = sqrt(g * domain.L.z), cᶻ = 1e-2),
        stabilizing_dissipation = stabilizing_dissipation,
        coriolis = (f₀ = 0, β = 0),
        buoyancy = (αᵀ = 0,),
        boundary_tags = ((0, 0), (1, 1), (1, 2)),
        boundary_conditions = (OceanBC(Impenetrable(FreeSlip()), Insulating()),
                               OceanBC(Penetrable(FreeSlip()), Insulating()))
    )
    exactz!(model)
    # println("modifying the grid")
    # cartesianify!(model.grid)
    # We prepare a callback that periodically fetches the horizontal velocity and
    # tracer concentration for later animation,

    writer = JLD2Writer(model, filepath = name * ".jld2", overwrite_existing = true)
    cpu_grid = DiscontinuousSpectralElementGrid(domain, array_type=Array)
    # cartesianify!(cpu_grid)
    # println("modifying the grid again")
    start_time = time_ns()

    data_fetcher = EveryXSimulationTime(output_time_interval) do
        write!(writer)

        cpu_data = convert(Array, model.state.realdata)
        u = SpectralElementField(domain, cpu_grid, view(cpu_data, :, 1, :))

        # Print a helpful message
        step = @sprintf("Step: %d", current_step(model))
        time = @sprintf("time: %.2f", current_time(model))
        max_u = @sprintf("max|u|: %.6f", maximum(abs, u))

        elapsed = (time_ns() - start_time) * 1e-9
        wall_time = @sprintf("elapsed wall time: %.2f min", elapsed / 60)  

        isnan(maximum(abs, u)) && error("NaNs.") 

        @info "$step, $time, $max_u, $wall_time"
    end

    # and then run the simulation.

    model.solver_configuration.timeend = stop_time

    total_steps = ceil(Int, stop_time / time_step)
    @info @sprintf("Running a simulation of the instability of the Bickley jet (Δt=%.2e, steps=%d)", time_step, total_steps)

    try
        result = ClimateMachine.invoke!(model.solver_configuration;
                                        user_callbacks = [data_fetcher])
    catch err
        @warn "Simulation ended prematurely because $(sprint(showerror, err))"
    end

    return name
end

function visualize(name, contours=false)

    filepath = name * ".jld2"

    u_timeseries = OutputTimeSeries(:u, filepath)
    v_timeseries = OutputTimeSeries(:v, filepath)
    η_timeseries = OutputTimeSeries(:η, filepath)
    c_timeseries = OutputTimeSeries(:θ, filepath)

    u₀ = u_timeseries[1]
    domain = u₀.domain
    assembled_u₀ = assemble(u₀)
    x = assembled_u₀.x[:, 1, 1]
    y = assembled_u₀.y[1, :, 1]

    times = u_timeseries.times

    animation = @animate for i = 1:length(u_timeseries)

        @info "Plotting frame $i of $(length(u_timeseries))..."

        kwargs = (xlim = domain.x, ylim = domain.y, linewidth = 0, aspectratio = 1)

        u = assemble(u_timeseries[i]).data[:, :, 1]
        v = assemble(v_timeseries[i]).data[:, :, 1]
        η = assemble(η_timeseries[i]).data[:, :, 1]
        c = assemble(c_timeseries[i]).data[:, :, 1]

        if ~isnan(maximum(abs, u))

            s = @. sqrt(u^2 + v^2)

            ηmin = minimum(η)
            ηmax = maximum(η)
            ηlevels = range(ηmin, ηmax, length=31)

            kwargs = Dict(:xlabel => "x",
                          :ylabel => "y",
                          :aspectratio => 1,
                          :linewidth => 0,
                          :colorbar => true,
                          :xlims => (-domain.L.x/2, domain.L.x/2),
                          :ylims => (-domain.L.y/2, domain.L.y/2))

            plotter = contours ? contourf : heatmap

            kwargs[:clims] = (0, 1)
            contours && (kwargs[:levels] = range(0, 1, length=31))

            s_plot = plotter(x, y, clamp.(s, 0, 1)'; color = :thermal, kwargs...)

            kwargs[:clims] = (-1, 1)
            contours && (kwargs[:levels] = range(-1, 1, length=31))

            u_plot = plotter(x, y, clamp.(u, -1, 1)'; color = :balance, kwargs...)
            c_plot = plotter(x, y, clamp.(c, -1, 1)'; color = :thermal, kwargs...)
            #η_plot = plotter(x, y, clamp.(c, -1, 1)'; color = :thermal, kwargs...)

            u_title = @sprintf("u at t = %.2f", times[i])
            s_title = @sprintf("speed at t = %.2f", times[i])
            c_title = @sprintf("c at t = %.2f", times[i])

            plot(u_plot, s_plot, c_plot,
                 title = [u_title s_title c_title],
                 layout = (1, 3),
                 size = (1600, 400))

        end
    end

    gif(animation, name * ".gif", fps = 8)

    return nothing
end
##
include("StabilizingDissipations.jl")
using .StabilizingDissipations: StabilizingDissipation

# begin loop, Currently using vorticity criteria
#=
for DOF in (32)
for Np  in (3)# (2,3,4,5,6)
for diffusive_cfl in [1e-4] # [1e-1, 1e-4]
for reduction in [1e4] # [1e0, 1e4] 
for smoothness_exponent in [2] # [1, 10]
=#
DOF = 32
Np = 3
diffusive_cfl = 1e-4
reduction  = 1e8
smoothness_exponent = 2
ν = 0.0
array_type = Array
output_time_interval = 2
stabilizing_dissipation = nothing
stop_time = 200
Ne = round(Int, DOF / (Np+1))
println("Doing DOF=", DOF)
println("Polynomial Order=", Np)
time_step = 1.0 * effective_node_spacing(Ne, Np) / c 
κ = effective_node_spacing(Ne, Np)^2 / time_step
##
test_dissipation = StabilizingDissipation(minimum_node_spacing = effective_node_spacing(Ne, Np),
                                          time_step = time_step,
                                          Δu = 20000.0 / Np,
                                          Δθ = 10000.0 / Np,
                                          diffusive_cfl = diffusive_cfl,
                                          κʰ_min = diffusive_cfl * κ /reduction,
                                          νʰ_min = diffusive_cfl * κ /reduction,
                                          smoothness_exponent = smoothness_exponent)
# augment_name = @sprintf("_diffusive_cfl%0.1e_reduction%0.1e_smoothness_exponent%.1e", diffusive_cfl, reduction, smoothness_exponent)
augment_name = "_test80"
#=
name = run(Ne=Ne, Np=Np, 
           stabilizing_dissipation=test_dissipation,
           time_step = time_step,
           augment_name = augment_name)
=#
name = run(Ne=Ne, Np=Np, 
           time_step = time_step,
           augment_name = augment_name)
println("Done with " * name)
# visualize(name)
##
#=
end
end
end
end
end
=#

## Start here
augment_name = "_derivative_4"
name = @sprintf("climate_machine_unstable_bickley_jet_Ne%d_Np%d_ν%.1e_no_rotation", Ne, Np, ν)
name *= augment_name
ClimateMachine.Settings.array_type = array_type

# Domain
DOF = 32 * 2
Np = 4
Ne = round(Int, DOF / (Np+1))
println("Doing DOF=", DOF)
println("Polynomial Order=", Np)
time_step = 1.0 * effective_node_spacing(Ne, Np) / c
domain = RectangularDomain(Ne = (Ne, Ne, 1), Np = Np,
                            x = (-2π, 2π), y = (-2π, 2π), z = (0, 1),
                            periodicity = (true, true, false))

# Physical parameters:
g = Planet.grav(NonDimensionalParameters())

# Non-dimensional parameters
ϵ = 0.1 # Perturbation amplitude
ℓ = 0.5 # Perturbation width
k = 0.5 # Perturbation wavenumber

# Initial conditions: Jet/tracer + perturbations
uᵢ(x, y, z) = Bickley.U(y) + ϵ * Bickley.ũ(x, y, ℓ, k)
vᵢ(x, y, z) = ϵ * Bickley.ṽ(x, y, ℓ, k)
θᵢ(x, y, z) = Bickley.C(y, domain.L.y)

initial_conditions = InitialConditions(u=uᵢ, v=vᵢ, θ=θᵢ)

model = Ocean.HydrostaticBoussinesqSuperModel(
    domain = domain,
    time_step = time_step,
    initial_conditions = initial_conditions,
    parameters = NonDimensionalParameters(),
    turbulence_closure = (νʰ = ν, κʰ = ν, νᶻ = ν, κᶻ = ν),
    rusanov_wave_speeds = (cʰ = sqrt(g * domain.L.z), cᶻ = 1e-2),
    stabilizing_dissipation = stabilizing_dissipation,
    coriolis = (f₀ = 0, β = 0),
    buoyancy = (αᵀ = 0,),
    boundary_tags = ((0, 0), (1, 1), (1, 2)),
    boundary_conditions = (OceanBC(Impenetrable(FreeSlip()), Insulating()),
                            OceanBC(Penetrable(FreeSlip()), Insulating()))
)
println("modifying the grid")
#=
exactz!(model) # make the vertical coordinates all have the same value within an element
cartesianify!(model.grid) # remove crossterms
# remove vertical derivatives
model.grid.vgeo[:,9,:]  .= 0.0 
model.grid.sgeo[5, :, 5, :] .= 0.0
model.grid.sgeo[5, :, 6, :] .= 0.0
=#
# model.grid.vgeo[:,16,:] .= -0.0  # remove JcV terms which come in https://github.com/CliMA/ClimateMachine.jl/blob/433994334c547fb439917c84fb82548a85e81a5f/src/Numerics/DGMethods/DGModel_kernels.jl#L2090
# We prepare a callback that periodically fetches the horizontal velocity and
# tracer concentration for later animation,

writer = JLD2Writer(model, filepath = name * ".jld2", overwrite_existing = true)
cpu_grid = DiscontinuousSpectralElementGrid(domain, array_type=Array)
cartesianify!(cpu_grid)
println("modifying the grid again")
start_time = time_ns()

data_fetcher = EveryXSimulationTime(output_time_interval) do
    write!(writer)

    cpu_data = convert(Array, model.state.realdata)
    u = SpectralElementField(domain, cpu_grid, view(cpu_data, :, 1, :))

    # Print a helpful message
    step = @sprintf("Step: %d", current_step(model))
    time = @sprintf("time: %.2f", current_time(model))
    max_u = @sprintf("max|u|: %.6f", maximum(abs, u))

    elapsed = (time_ns() - start_time) * 1e-9
    wall_time = @sprintf("elapsed wall time: %.2f min", elapsed / 60)  

    isnan(maximum(abs, u)) && error("NaNs.") 

    @info "$step, $time, $max_u, $wall_time"
end

# and then run the simulation.

model.solver_configuration.timeend = stop_time

total_steps = ceil(Int, stop_time / time_step)
@info @sprintf("Running a simulation of the instability of the Bickley jet (Δt=%.2e, steps=%d)", time_step, total_steps)
##
field = model.fields.u
u1 = assemble(field).data[:,:,1]
u2 = assemble(field).data[:,:,2]
u3 = assemble(field).data[:,:,3]
##
try
    result = ClimateMachine.invoke!(model.solver_configuration;
                                    user_callbacks = [data_fetcher])
catch err
    @warn "Simulation ended prematurely because $(sprint(showerror, err))"
end
