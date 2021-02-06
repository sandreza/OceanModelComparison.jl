using JLD2
using ClimateMachine
ClimateMachine.init()
using ClimateMachine.Mesh.Grids: DiscontinuousSpectralElementGrid
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.MPIStateArrays
using ClimateMachine.DGMethods.NumericalFluxes
using MPI, Printf
using LinearAlgebra

include("convenience.jl")

DOFs = [32]
Ns = [1, 2, 3, 4]
Novers = [0, 1]
fluxes = [RusanovNumericalFlux(), RoeNumericalFlux()]
periodicity = [true] # [true, false]

DOFs = [32]
Ns = [1, 2, 3, 4]
Novers = [1]
fluxes = [RoeNumericalFlux(), RusanovNumericalFlux()]
periodicity = [true]

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
    f = jldopen(name * ".jld2", "r+")
    simtime = f["simulationtime"]
    threadnum = f["threads"]
    array = f["arraytype"]
    if array == "Array" 
        archstring = " on the CPU with " * string(threadnum) * " threads"
    else
        archstring = " on the GPU"
    end
    prettyname = nameprettifier(name)
    println("For simulation " * name)
    println("The simulation time was " * @sprintf("%0.2f", simtime) * " seconds for " * prettyname * archstring)
    # get old grid
    newgrid = f["grid"]
    # interpolate
    M = view(newgrid.vgeo, :, newgrid.Mid, :)
    for timeindex in 100:100
        println("time " , timeindex)
    Q = f[string(timeindex)]
    for stateindex in 1:5
        field = Q[:,stateindex,:] # this is the tracer
        μQ = sum(M .* field) / sum(M) 
        l1Q = sum(M .* abs.(field)) / sum(M)
        val = μQ / l1Q
        println("For stateindex " * string(stateindex))
        println("the mean value is ", μQ)
        println("the l1 norm is ", l1Q)
        println("the relative mean value is ", val)
    end
        println("---------")
    end

    println("------------------------------")
    close(f)
end