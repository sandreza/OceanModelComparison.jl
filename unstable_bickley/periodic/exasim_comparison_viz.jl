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
for N in 1:4
Ne = round(Int, DOF / (N+1))
filename = "compare_p" * string(N) * "_N" * string(Ne)
f = jldopen(filename * ".jld2", "r+")
dg_grid = f["grid"]
gridhelper = GridHelper(dg_grid)   
x, y, z = coordinates(dg_grid)
xC, yC, zC = cellcenters(dg_grid)
ϕ =  ScalarField(copy(x), gridhelper)
newx = range(-2π, 2π, length = 128)
newy = range(-2π, 2π, length = 128)
ρθ = zeros(length(newx), length(newy), 101)
tic = time()
for i in 0:100
    Q = f[string(i)]
    ϕ .= Q[:,4,:]
    ρθ[:,:,i+1] = ϕ(newx, newy, threads = true)
end
push!(states, ρθ)
toc = time()
println("Interpolation took $(toc-tic) seconds")
close(f)
end
##
scene, layout = layoutscene()
lscene = layout[2:4,2:4] = LScene(scene)

time_slider = LSlider(scene, range = Int.(range(1, 101, length=101)), startvalue = 1)
time_node = time_slider.value
state1 = @lift(states[1][:,:,$time_node])
state2 = @lift(states[2][:,:,$time_node])
state3 = @lift(states[3][:,:,$time_node])
state4 = @lift(states[4][:,:,$time_node])
clims = (-1,1)
heatmap1 = heatmap!(lscene, 0..1, 1..2, state1, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene, 1..2, 1..2, state2, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene, 0..1, 0..1, state3, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene, 1..2, 0..1, state4, colorrange = clims, colormap = to_colormap(:balance))

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