

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

ϕ.data[:,gridhelper.element.permutation] .= exasim["1"][:,2,:]

newx = range(-2π, 2π, length = 128)
newy = range(-2π, 2π, length = 128)
ρθ = zeros(length(newx), length(newy), 101)
ρθ[:,:,1] .= ϕ(newx, newy)

heatmap(0..1, 0..1, ρθ[:,:,1])