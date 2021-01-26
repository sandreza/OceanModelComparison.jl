

DOF = 32
N = 1
Ne = round(Int, DOF / (N+1))
filename = "compare_p" * string(N) * "_N" * string(Ne)
climatemachine = jldopen(filename * ".jld2", "r+")
filename = "exasim_p" * string(N) * "_N" * string(Ne)
exasim  = jldopen(filename * ".jld2", "r+")

cmt = climatemachine["6threadsimulationtime"]
et = exasim["simtime"][2]
ethr = exasim["simtime"][1]
println("ClimateMachine on 6threads took $(cmt) seconds")
println("Exasim on $(ethr) took $(et) seconds")
if et < cmt
    println("Exasim was faster by $(cmt-et) seconds")
else
    println("ClimateMachine was faster by $(et - cmt) seconds")
end
println("-------------------")
dg_grid = climatemachine["grid"]
gridhelper = GridHelper(dg_grid)   
x, y, z = coordinates(dg_grid)
xC, yC, zC = cellcenters(dg_grid)
ϕ =  ScalarField(copy(x), gridhelper)
##

DOF = 32

states = []
estates = []
for N in 1:4
Ne = round(Int, DOF / (N+1))
filename = "compare_p" * string(N) * "_N" * string(Ne)
climatemachine = jldopen(filename * ".jld2", "r+")
filename = "exasim_p" * string(N) * "_N" * string(Ne)
exasim  = jldopen(filename * ".jld2", "r+")

cmt = climatemachine["6threadsimulationtime"]
et = exasim["simtime"][2]
ethr = exasim["simtime"][1]
println("ClimateMachine on 6threads took $(cmt) seconds")
println("Exasim on $(ethr) thread took $(et) seconds")
if et < cmt
    println("Exasim was faster by $(cmt-et) seconds")
else
    println("ClimateMachine was faster by $(et - cmt) seconds")
end
dg_grid = climatemachine["grid"]
gridhelper = GridHelper(dg_grid)   
x, y, z = coordinates(dg_grid)
xC, yC, zC = cellcenters(dg_grid)
ϕ =  ScalarField(copy(x), gridhelper)

newx = range(-2π, 2π, length = 128)
newy = range(-2π, 2π, length = 128)
ρθ = zeros(length(newx), length(newy), 100)
e_ρθ = zeros(length(newx), length(newy), 100)
tic = time()
for i in 1:100
    Q = climatemachine[string(i)]
    ϕ .= Q[:,4,:]
    ρθ[:,:,i] = ϕ(newx, newy, threads = true)
    ϕ.data[:, gridhelper.element.permutation] .= exasim[string(i)][:,4,:]
    e_ρθ[:,:,i] = ϕ(newx, newy, threads = true)
end
push!(states, ρθ)
push!(estates, e_ρθ)
toc = time()
println("Interpolation took $(toc-tic) seconds")
println("-------------------")
close(climatemachine)
close(exasim)
end
#=
newx = range(-2π, 2π, length = 128)
newy = range(-2π, 2π, length = 128)
ρθ = zeros(length(newx), length(newy), 100)
tic = time()
for i in 1:100
    ϕ.data[:, gridhelper.element.permutation] .= exasim[string(i)][:,4,:]
    ρθ[:,:,i] .= ϕ(newx, newy)
end
toc = time()
println("The time for interpolation is $(toc - tic)")
=#
##
scene, layout = layoutscene()
lscene = layout[2:4,2:4] = LScene(scene)
layout[1, 2:4] = LText(scene, "Exasim")
time_slider = LSlider(scene, range = Int.(range(1, 100, length=100)), startvalue = 1)
time_node = time_slider.value

estate1 = @lift(estates[1][:,:,$time_node])
estate2 = @lift(estates[2][:,:,$time_node])
estate3 = @lift(estates[3][:,:,$time_node])
estate4 = @lift(estates[4][:,:,$time_node])
clims = (-1,1)
heatmap1 = heatmap!(lscene, 0..1, 1..2, estate1, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene, 1..2, 1..2, estate2, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene, 0..1, 0..1, estate3, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene, 1..2, 0..1, estate4, colorrange = clims, colormap = to_colormap(:balance))

cbar = LColorbar(scene, heatmap1)
cbar.height = Relative(1/3)
cbar.width = Relative(1/3)
cbar.halign = :left
cbar.labelsize = 50

layout[2:4, 1] = vgrid!(
    LText(scene, "Time", width = nothing),
    time_slider,
    cbar
)
display(scene)