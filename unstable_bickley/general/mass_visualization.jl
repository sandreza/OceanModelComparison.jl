using JLD2
using ClimateMachine
ClimateMachine.init()
using ClimateMachine.Mesh.Grids: DiscontinuousSpectralElementGrid
polynomialorders(::DiscontinuousSpectralElementGrid{T, dim, N}) where {T, dim, N} = N
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.MPIStateArrays
using ClimateMachine.DGMethods.NumericalFluxes
using MPI
using LinearAlgebra
include(pwd() * "/unstable_bickley/periodic/imperohooks.jl")
include(pwd() * "/unstable_bickley/periodic/vizinanigans.jl")
include(pwd() * "/unstable_bickley/general/convenience.jl")

DOFs = [32]
Ns = [1, 2, 3, 4]
Novers = [0, 1]
fluxes = [RusanovNumericalFlux(), RoeNumericalFlux()]
periodicity = [false] # [true, false]

states = []
namelist = []
for DOF in DOFs, Nover in Novers, flux in fluxes, periodic in periodicity, N in Ns
    name = just_generate_name(DOF, N, Nover, flux, periodic)
    push!(namelist, name)
    print(name)
    println(" ")
end
##
for name in namelist
    f = jldopen(name * ".jld2")
    simtime = f["simulationtime"]
    threadnum = f["threads"]
    array = f["arraytype"]
    if array == "Array" 
        archstring = " on the CPU with " * string(threadnum) * " threads"
    else
        archstring = " on the GPU"
    end
    prettyname = nameprettifier(name)
    println("The simulation time was " * @sprintf("%0.2f", simtime) * " seconds for " * prettyname * archstring)
    println("------------------------------")
    # get old grid
    newgrid = f["grid"]
    gridhelper = GridHelper(newgrid)  
    x, y, z = coordinates(newgrid)
    ϕ =  ScalarField(copy(x), gridhelper)
    # new grid
    newx = range(-2π, 2π, length = DOFs[1] * 2)
    newy = range(-2π, 2π, length = DOFs[1] * 2)
    ρθ = zeros(length(newx), length(newy), 100)
    # interpolate
    M = view(newgrid.vgeo, :, newgrid.Mid, :)
    for i in 1:100
        Q = f[string(i)]
        ϕ .= Q[:,4,:]
        ρθ[:,:,i] = ϕ(newx, newy, threads = true)
    end
    Q = f[string(100)]
    val = sum(M .* ϕ.data)/sum(M .* abs.(ϕ.data))
    if val ≥ eps(1e0)
        println(name * " not conserved to machine precision at " * string(100))
        println("the relative value is ", val)
        println("the value is ", sum(M .* ϕ.data))
    end
    push!(states, ρθ)
    close(f)
end

##
DOF = DOFs[1]
resolution = (2880, 1998)
interpolate = false
scene, layout = layoutscene(resolution = resolution )
lscene1 = layout[2:4, 3:5] = LScene(scene)
lscene2 = layout[2:4, 6:8] = LScene(scene)
lscene3 = layout[5:7, 3:5] = LScene(scene)
lscene4 = layout[5:7, 6:8] = LScene(scene)

layout[1, 3:5] = Label(scene, "Rusanov", textsize = 50)
layout[1, 6:8] = Label(scene, "Roe", textsize = 50)
layout[3, 2] = Label(scene, "Underintegration", textsize = 50)
layout[6, 2] = Label(scene, "Overintegration", textsize = 50)
layout[1,1] = Label(scene, "ρθ, DOF = $(DOF)^2", textsize = 50)

#=
layout[1,1] = Label(scene, "ρθ, Roe Fluxes", textsize = 50)
layout[1, 3:5-1] = Label(scene, "Underintegration", textsize = 50)
layout[1, 6:8-1] = Label(scene, "Overintegration", textsize = 50)
layout[3, 2] = Label(scene, "32²", textsize = 50)
layout[6, 2] = Label(scene, "128²", textsize = 50)
=#

time_slider = Slider(scene, range = Int.(range(1, 100, length=100)), startvalue = 1)
time_node = time_slider.value

# roe
state1 = @lift(states[1][:,:,$time_node])
state2 = @lift(states[2][:,:,$time_node])
state3 = @lift(states[3][:,:,$time_node])
state4 = @lift(states[4][:,:,$time_node])
clims = (-1,1)

heatmap1 = heatmap!(lscene1, 0..1, 1.1..2.1, state1, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene1, 1..2, 1.1..2.1, state2, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene1, 0..1, 0..1, state3, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene1, 1..2, 0..1, state4, colorrange = clims, colormap = :balance, interpolate = interpolate)

state5 = @lift(states[5][:,:,$time_node])
state6 = @lift(states[6][:,:,$time_node])
state7 = @lift(states[7][:,:,$time_node])
state8 = @lift(states[8][:,:,$time_node])

heatmap1 = heatmap!(lscene2, 0..1, 1.1..2.1, state5, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene2, 1..2, 1.1..2.1, state6, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene2, 0..1, 0..1, state7, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene2, 1..2, 0..1, state8, colorrange = clims, colormap = :balance, interpolate = interpolate)

state9  = @lift(states[9][:,:,$time_node])
state10 = @lift(states[10][:,:,$time_node])
state11 = @lift(states[11][:,:,$time_node])
state12 = @lift(states[12][:,:,$time_node])

heatmap1 = heatmap!(lscene3, 0..1, 1.1..2.1, state9, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene3, 1..2, 1.1..2.1, state10, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene3, 0..1, 0..1, state11, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene3, 1..2, 0..1, state12, colorrange = clims, colormap = :balance, interpolate = interpolate)

state13  = @lift(states[13][:,:,$time_node])
state14 = @lift(states[14][:,:,$time_node])
state15 = @lift(states[15][:,:,$time_node])
state16 = @lift(states[16][:,:,$time_node])

heatmap1 = heatmap!(lscene4, 0..1, 1.1..2.1, state13, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene4, 1..2, 1.1..2.1, state14, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene4, 0..1, 0..1, state15, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene4, 1..2, 0..1, state16, colorrange = clims, colormap = :balance, interpolate = interpolate)

#=
heatmap1 = heatmap!(lscene1, 0..1, 0..1, state1, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene1, 1..2, 0..1, state2, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene1, 2..3, 0..1, state3, colorrange = clims, colormap = :balance, interpolate = interpolate)
heatmap!(lscene1, 3..4, 0..1, state4, colorrange = clims, colormap = :balance, interpolate = interpolate)
=#


cbar = Colorbar(scene, heatmap1)
cbar.height = Relative(1/3)
cbar.width = Relative(1/3)
cbar.halign = :left
cbar.labelsize = 50

slidertext = @lift("Time t = " * string(2 *  $time_node))
layout[2:6, 1] = vgrid!(
    Label(scene, slidertext, width = nothing),
    time_slider,
    cbar,
    Label(scene, "top left p=1, top right p=2"),
    Label(scene, "bottom left p=3, bottom right p=4"),
)
display(scene)

##
seconds = 15
fps = 10
frames = round(Int, fps * seconds )
record(scene, pwd() * "/512dof.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end

##
record(scene, "512dof.mp4", 1:100, framerate=10) do n
    time_node[] = n
end