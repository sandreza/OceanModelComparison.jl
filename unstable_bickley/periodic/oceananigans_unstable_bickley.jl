using JLD2
using Plots
using Printf
using Statistics

using Oceananigans
using Oceananigans.Advection
using Oceananigans.AbstractOperations
using Oceananigans.OutputWriters
using Oceananigans.Grids
using Oceananigans.Fields
using Oceananigans.Forcings

include("Bickley.jl")

using .Bickley

function run(; Nh = 128,
               output_time_interval = 2,
               stop_time = 200,
               arch = CPU(),
               ν = 0,
               advection = WENO5())

    name = "oceananigans_unstable_bickley_jet_Nh$(Nh)_$(typeof(advection).name.wrapper)"

    grid = RegularCartesianGrid(size=(Nh, Nh, 1),
                                x = (-2π, 2π), y=(-2π, 2π), z=(0, 1),
                                topology = (Periodic, Periodic, Bounded))

    model = IncompressibleModel(architecture = arch,
                                 timestepper = :RungeKutta3, 
                                   advection = advection,
                                        grid = grid,
                                     tracers = :c,
                                     closure = IsotropicDiffusivity(ν=ν, κ=ν),
                                    buoyancy = nothing)

    # ** Initial conditions **
    #
    # u, v: Large-scale jet + vortical perturbations
    #    c: Sinusoid
    
    # Parameters
    ϵ = 0.1 # perturbation magnitude
    ℓ = 0.5 # Gaussian width
    k = 0.5 # Sinusoidal wavenumber

    # Total initial conditions
    uᵢ(x, y, z) = Bickley.U(y) + ϵ * Bickley.ũ(x, y, ℓ, k)
    vᵢ(x, y, z) = ϵ * Bickley.ṽ(x, y, ℓ, k)
    cᵢ(x, y, z) = Bickley.C(y, grid.Ly)

    set!(model, u=uᵢ, v=vᵢ, c=cᵢ)
    
    progress(sim) = @info(@sprintf("Iter: %d, time: %.1f, Δt: %.3f, max|u|: %.2f",
                                   sim.model.clock.iteration, sim.model.clock.time,
                                   sim.Δt.Δt, maximum(abs, u.data.parent)))

    wizard = TimeStepWizard(cfl=1.0, Δt=1e-1, max_change=1.1, max_Δt=10.0)

    simulation = Simulation(model, Δt=wizard, stop_time=stop_time,
                            iteration_interval=10, progress=progress)

    # Output: primitive fields + computations
    u, v, w, c = primitives = merge(model.velocities, model.tracers)

   
    outputs = primitives

    save_grid = (file, model) -> file["serialized/grid"] = model.grid

    simulation.output_writers[:fields] =
        JLD2OutputWriter(model, merge(model.velocities, model.tracers),
                                schedule = TimeInterval(output_time_interval),
                                init = save_grid,
                                prefix = name * "_fields",
                                force = true)



    @info "Running a simulation of an unstable Bickley jet with $(Nh)² degrees of freedom..."

    run!(simulation)

    return name
end

function analyze(name)
    @info "Analyzing the results of an unstable Bickley jet simulation..."

    statistics_file = jldopen(name * "_statistics.jld2")

    iterations = parse.(Int, keys(statistics_file["timeseries/t"]))
    grid = statistics_file["serialized/grid"]

    ∇c² = [statistics_file["timeseries/∇c²/$iter"][1, 1, 1] for iter in iterations]
    c²  = [statistics_file["timeseries/c²/$iter"][1, 1, 1] for iter in iterations]
    t   = [statistics_file["timeseries/t/$iter"][1, 1, 1] for iter in iterations]

    close(statistics_file)

    plot(t, 2π / grid.Ly .* ∇c² ./ c²)
    savefig(name * ".png")

    return nothing
end

function visualize(name, contours=false)
    @info "Making a fun movie about an unstable Bickley jet..."

    fields_file = jldopen(name * "_fields.jld2")

    iterations = parse.(Int, keys(fields_file["timeseries/t"]))
    grid = fields_file["serialized/grid"]
    #=
    xu, yu, zu = nodes((Face, Cell, Cell), grid)
    xω, yω, zω = nodes((Face, Face, Cell), grid)
    xc, yc, zc = nodes((Cell, Cell, Cell), grid)
    =#
    anim = @animate for (i, iteration) in enumerate(iterations)

        @info "    Plotting frame $i from iteration $iteration..."
        
        t = fields_file["timeseries/t/$iteration"]
        u = fields_file["timeseries/u/$iteration"][:, :, 1]
        c = fields_file["timeseries/c/$iteration"][:, :, 1]

        kwargs = Dict(:xlabel => "x",
                      :ylabel => "y",
                      :aspectratio => 1,
                      :linewidth => 0,
                      :colorbar => true,
                      :clims => (-1, 1),
                      :xlims => (-grid.Lx/2, grid.Lx/2),
                      :ylims => (-grid.Ly/2, grid.Ly/2))

        contours && (kwargs[:levels] = range(-1, 1, length=31))
        plotter = contours ? contourf : heatmap

        u_plot = plotter(xu, yu, clamp.(u, -1, 1)'; color = :balance, kwargs...)
        c_plot = plotter(xc, yc, clamp.(c, -1, 1)'; color = :balance, kwargs...)

        u_title = @sprintf("u at t = %.1f", t)
        c_title = @sprintf("c at t = %.1f", t)

        plot(u_plot, c_plot,
             title = [u_title c_title],
             layout = (1, 2),
             size = (1600, 400))
    end

    gif(anim, name * ".gif", fps = 8)

    return nothing
end
##
Nh = 32
advection = WENO5()

tic = time()
name = run(Nh=Nh, advection=advection)
toc = time()
println(toc - tic)
##
fields_file = jldopen(name * "_fields.jld2")

iterations = parse.(Int, keys(fields_file["timeseries/t"]))
grid = fields_file["serialized/grid"]
ct = zeros(Nh, Nh, length(iterations))
for (i, iteration) in enumerate(iterations) 
        ct[:,:,i] = fields_file["timeseries/c/$iteration"][:, :, 1]
end

##
comparewith = false # need to run roeviz.jl first for comparison
using GLMakie
resolution = (2404, 1308)
interpolate = true
scene, layout = layoutscene(resolution = resolution )
lscene = layout[2:4,2:4] = LScene(scene)
if comparewith
lscene2 = layout[2:4,5:7] = LScene(scene)
end
layout[1,1] = LText(scene, "θ, DOF = $(Nh)^2", textsize = 50)

if comparewith
    layout[1,2:4] = LText(scene, "Oceananigans", textsize = 50)
layout[1,5:7] = LText(scene, "Climate Machine", textsize = 50)
end
time_slider = LSlider(scene, range = Int.(range(1, 100, length=100)), startvalue = 1)
time_node = time_slider.value

estate = @lift(ct[:,:,$time_node + 1])
clims = (-1,1)
heatmap1 = GLMakie.heatmap!(lscene, 0..1, 0..1, estate, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)

if comparewith
estate2 = @lift(oiroestates[4][:,:,$time_node])
clims = (-1,1)
heatmap1 = GLMakie.heatmap!(lscene2, 0..1, 0..1, estate2, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
end
estate2 = @lift(oiroestates[2][:,:,$time_node])
heatmap1 = GLMakie.heatmap!(lscene, 1..2, 0..1, estate2, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)

cbar = LColorbar(scene, heatmap1)
cbar.height = Relative(1/3)
cbar.width = Relative(1/3)
cbar.halign = :left
cbar.labelsize = 50

slidertext = @lift("Time t = " * string(2 *  $time_node))
layout[2:4, 1] = vgrid!(
    LText(scene, slidertext, width = nothing),
    time_slider,
    cbar,
)
display(scene)
##
seconds = 15
fps = 10
frames = round(Int, fps * seconds )
record(scene, pwd() * "/sideby.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end

##
#=
for Nh in (128)
    #for advection in (CenteredSecondOrder(), CenteredFourthOrder(), UpwindBiasedThirdOrder(), UpwindBiasedFifthOrder(), WENO5())
    for advection in (WENO5(),)
        name = run(Nh=Nh, advection=advection)
        # analyze(name)
        visualize(name)
    end
end
=#
