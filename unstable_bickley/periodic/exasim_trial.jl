using JLD2
using ClimateMachine
ClimateMachine.init()
using ClimateMachine.Mesh.Grids: DiscontinuousSpectralElementGrid
polynomialorders(::DiscontinuousSpectralElementGrid{T, dim, N}) where {T, dim, N} = Tuple([N for i in 1:dim])
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.MPIStateArrays
using MPI
using LinearAlgebra
include(pwd() * "/unstable_bickley/periodic/imperohooks.jl")
include(pwd() * "/unstable_bickley/periodic/vizinanigans2.jl")

DOF = 32

states = []
estates = []
istates = []
for N in 1:4
Ne = round(Int, DOF / (N+1))
filename = "compare_p" * string(N) * "_N" * string(Ne)
climatemachine = jldopen(filename * ".jld2", "r+")
filename = "exasim_p" * string(N) * "_N" * string(Ne) 
exasim  = jldopen(filename * ".jld2", "r+")

filename = "overint_p" * string(N) * "_N" * string(Ne)
oiclimatemachine  = jldopen(filename * ".jld2", "r+")

cmt = climatemachine["6threadsimulationtime"]
et = exasim["simtime"][2]
ethr = exasim["simtime"][1]
oicmt = oiclimatemachine["6threadsimulationtime"]
println("ClimateMachine on 6threads took $(cmt) seconds")
println("Exasim on $(ethr) thread took $(et) seconds")
println("oi ClimateMachine on 6threads took $(oicmt) seconds")
if et < cmt
    println("Exasim was faster by $(cmt-et) seconds")
else
    println("ClimateMachine was faster by $(et - cmt) seconds")
    println("OI ClimateMachine was faster by $(et - oicmt) seconds")
end

oi_grid = oiclimatemachine["grid"]
oigridhelper = GridHelper(oi_grid)  
dg_grid = climatemachine["grid"]
gridhelper = GridHelper(dg_grid)   
x, y, z = coordinates(dg_grid)
xC, yC, zC = cellcenters(dg_grid)
ϕ =  ScalarField(copy(x), gridhelper)

oix, oiy, oiz = coordinates(oi_grid)
oiϕ = ScalarField(copy(oix), oigridhelper)

newx = range(-2π, 2π, length = 64)
newy = range(-2π, 2π, length = 64)
ρθ = zeros(length(newx), length(newy), 100)
e_ρθ = zeros(length(newx), length(newy), 100)
oi_ρθ = zeros(length(newx), length(newy), 100)
tic = time()
for i in 1:100
    Q = climatemachine[string(i)]
    ϕ .= Q[:,4,:]
    ρθ[:,:,i] = ϕ(newx, newy, threads = true)
    ϕ.data[:, gridhelper.element.permutation] .= exasim[string(i)][:,4,:]
    e_ρθ[:,:,i] = ϕ(newx, newy, threads = true)
    # oi
    Q = oiclimatemachine[string(i)]
    oiϕ .= Q[:,4,:]
    oi_ρθ[:,:,i] = oiϕ(newx, newy, threads = true)
end
push!(states, ρθ)
push!(estates, e_ρθ)
push!(istates, oi_ρθ)
toc = time()
println("Interpolation took $(toc-tic) seconds")
println("-------------------")
close(climatemachine)
close(exasim)
close(oiclimatemachine)
end

##
resolution = (1956+400, 852)
interpolate = false
scene, layout = layoutscene(resolution = resolution )
lscene = layout[2:4,2:4] = LScene(scene)
lscene2 = layout[2:4, 5:7] = LScene(scene)
lscene3 = layout[2:4, 8:10] = LScene(scene)
layout[1, 2:4] = LText(scene, "Exasim ρθ", textsize = 50)
layout[1, 5:7] = LText(scene, "ClimateMachine ρθ", textsize = 50)
layout[1, 8:10] = LText(scene, "OI ClimateMachine ρθ", textsize = 50)
layout[1,1] = LText(scene, "DOF = $(DOF)^2", textsize = 50)
time_slider = LSlider(scene, range = Int.(range(1, 100, length=100)), startvalue = 1)
time_node = time_slider.value

# Exasim
estate1 = @lift(estates[1][:,:,$time_node])
estate2 = @lift(estates[2][:,:,$time_node])
estate3 = @lift(estates[3][:,:,$time_node])
estate4 = @lift(estates[4][:,:,$time_node])
clims = (-1,1)
heatmap1 = heatmap!(lscene, 0..1, 1..2, estate1, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene, 1..2, 1..2, estate2, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene, 0..1, 0..1, estate3, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene, 1..2, 0..1, estate4, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)

# ClimateMachine
state1 = @lift(states[1][:,:,$time_node])
state2 = @lift(states[2][:,:,$time_node])
state3 = @lift(states[3][:,:,$time_node])
state4 = @lift(states[4][:,:,$time_node])
clims = (-1,1)
heatmap1 = heatmap!(lscene2, 0..1, 1..2, state1, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene2, 1..2, 1..2, state2, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene2, 0..1, 0..1, state3, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene2, 1..2, 0..1, state4, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)

# Overintegration
istate1 = @lift(istates[1][:,:,$time_node])
istate2 = @lift(istates[2][:,:,$time_node])
istate3 = @lift(istates[3][:,:,$time_node])
istate4 = @lift(istates[4][:,:,$time_node])
clims = (-1,1)
heatmap1 = heatmap!(lscene3, 0..1, 1..2, istate1, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene3, 1..2, 1..2, istate2, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene3, 0..1, 0..1, istate3, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene3, 1..2, 0..1, istate4, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)

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
    LText(scene, "top left p=1, top right p=2"),
    LText(scene, "bottom left p=3, bottom right p=4"),
)
display(scene)

##
record_interaction = true
seconds = 20
fps = 10
frames = round(Int, fps * seconds )
if record_interaction
record(scene, pwd() * "/comparison.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end
end