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

DOF = 128

states = []
estates = []
istates = []
oiroestates = []
for N in [1,2,3, 4]
Ne = round(Int, DOF / (N+1))
filename = "compare_p" * string(N) * "_N" * string(Ne)
climatemachine = jldopen(filename * ".jld2", "r+")

filename = "roe_p" * string(N) * "_N" * string(Ne) 
roe  = jldopen(filename * ".jld2", "r+")

filename = "overint_p" * string(N) * "_N" * string(Ne)
oiclimatemachine  = jldopen(filename * ".jld2", "r+")

filename = "roe_overint_p" * string(N) * "_N" * string(Ne)
oiroe  = jldopen(filename * ".jld2", "r+")

cmt = climatemachine["6threadsimulationtime"]
et = roe["6threadsimulationtime"]
oicmt = oiclimatemachine["6threadsimulationtime"]
oiet = oiroe["6threadsimulationtime"]

println("Rusanov on 6 threads took $(cmt) seconds")
println("Roe on 6 threads took $(et) seconds")
println("OI Rusanov on 6threads took $(oicmt) seconds")
println("OI Roe on 6threads took $(oicmt) seconds")

oi_grid = oiclimatemachine["grid"]
oigridhelper = GridHelper(oi_grid)  
dg_grid = climatemachine["grid"]
gridhelper = GridHelper(dg_grid)   
x, y, z = coordinates(dg_grid)
xC, yC, zC = cellcenters(dg_grid)
ϕ =  ScalarField(copy(x), gridhelper)

oix, oiy, oiz = coordinates(oi_grid)
oiϕ = ScalarField(copy(oix), oigridhelper)

newx = range(-2π, 2π, length = DOF * 2)
newy = range(-2π, 2π, length = DOF * 2)
ρθ = zeros(length(newx), length(newy), 100)
e_ρθ = zeros(length(newx), length(newy), 100)
oi_ρθ = zeros(length(newx), length(newy), 100)
oi_e_ρθ = zeros(length(newx), length(newy), 100)
tic = time()
for i in 1:100
    Q = climatemachine[string(i)]
    ϕ .= Q[:,4,:]
    ρθ[:,:,i] = ϕ(newx, newy, threads = true)
    ϕ.data .= roe[string(i)][:,4,:]
    e_ρθ[:,:,i] = ϕ(newx, newy, threads = true)
    # oi
    Q = oiclimatemachine[string(i)]
    oiϕ .= Q[:,4,:]
    oi_ρθ[:,:,i] = oiϕ(newx, newy, threads = true)
    oiϕ.data .= oiroe[string(i)][:,4,:]
    oi_e_ρθ[:,:,i] = oiϕ(newx, newy, threads = true)
end
push!(states, ρθ)
push!(estates, e_ρθ)
push!(istates, oi_ρθ)
push!(oiroestates, oi_e_ρθ)
toc = time()
println("Interpolation took $(toc-tic) seconds")
println("-------------------")
close(climatemachine)
close(roe)
close(oiclimatemachine)
close(oiroe)
end

##
resolution = (2404, 1308)
interpolate = false
scene, layout = layoutscene(resolution = resolution )
lscene = layout[2:4,3:5] = LScene(scene)
lscene2 = layout[2:4, 6:8] = LScene(scene)
lscene3 = layout[5:7, 3:5] = LScene(scene)
lscene4 = layout[5:7, 6:8] = LScene(scene)

layout[1, 3:5] = LText(scene, "Rusanov", textsize = 50)
layout[1, 6:8] = LText(scene, "Roe", textsize = 50)
layout[3, 2] = LText(scene, "Underintegration", textsize = 50)
layout[6, 2] = LText(scene, "Overintegration", textsize = 50)
layout[1,1] = LText(scene, "ρθ, DOF = $(DOF)^2", textsize = 50)
time_slider = LSlider(scene, range = Int.(range(1, 100, length=100)), startvalue = 1)
time_node = time_slider.value

# roe
estate1 = @lift(estates[1][:,:,$time_node])
estate2 = @lift(estates[2][:,:,$time_node])
estate3 = @lift(estates[3][:,:,$time_node])
estate4 = @lift(estates[4][:,:,$time_node])
clims = (-1,1)
heatmap1 = heatmap!(lscene2, 0..1, 1..2, estate1, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene2, 1..2, 1..2, estate2, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene2, 0..1, 0..1, estate3, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene2, 1..2, 0..1, estate4, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)

# ClimateMachine
state1 = @lift(states[1][:,:,$time_node])
state2 = @lift(states[2][:,:,$time_node])
state3 = @lift(states[3][:,:,$time_node])
state4 = @lift(states[4][:,:,$time_node])
clims = (-1,1)
heatmap1 = heatmap!(lscene, 0..1, 1..2, state1, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene, 1..2, 1..2, state2, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene, 0..1, 0..1, state3, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene, 1..2, 0..1, state4, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)

# Overintegration Rusanov
istate1 = @lift(istates[1][:,:,$time_node])
istate2 = @lift(istates[2][:,:,$time_node])
istate3 = @lift(istates[3][:,:,$time_node])
istate4 = @lift(istates[4][:,:,$time_node])
clims = (-1,1)
heatmap1 = heatmap!(lscene3, 0..1, 1..2, istate1, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene3, 1..2, 1..2, istate2, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene3, 0..1, 0..1, istate3, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene3, 1..2, 0..1, istate4, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)

# Overintegration Roe
oiroestate1 = @lift(oiroestates[1][:,:,$time_node])
oiroestate2 = @lift(oiroestates[2][:,:,$time_node])
oiroestate3 = @lift(oiroestates[3][:,:,$time_node])
oiroestate4 = @lift(oiroestates[4][:,:,$time_node])
clims = (-1,1)
heatmap1 = heatmap!(lscene4, 0..1, 1..2, oiroestate1, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene4, 1..2, 1..2, oiroestate2, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene4, 0..1, 0..1, oiroestate3, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)
heatmap!(lscene4, 1..2, 0..1, oiroestate4, colorrange = clims, colormap = to_colormap(:balance), interpolate = interpolate)

cbar = LColorbar(scene, heatmap1)
cbar.height = Relative(1/3)
cbar.width = Relative(1/3)
cbar.halign = :left
cbar.labelsize = 50

slidertext = @lift("Time t = " * string(2 *  $time_node))
layout[2:6, 1] = vgrid!(
    LText(scene, slidertext, width = nothing),
    time_slider,
    cbar,
    LText(scene, "top left p=1, top right p=2"),
    LText(scene, "bottom left p=3, bottom right p=4"),
)
display(scene)

##
seconds = 15
fps = 10
frames = round(Int, fps * seconds )
record(scene, pwd() * "/roe_overint.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end