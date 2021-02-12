using Test
using JLD2
using ClimateMachine
ClimateMachine.init(disable_gpu = true)
using ClimateMachine.ODESolvers
using ClimateMachine.Mesh.Filters
using ClimateMachine.VariableTemplates
using ClimateMachine.Mesh.Grids: polynomialorders
using ClimateMachine.Ocean

using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.BalanceLaws:
    vars_state, Prognostic, Auxiliary, number_states
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.MPIStateArrays
using ClimateMachine.VTK

using MPI
using LinearAlgebra
using StaticArrays
using Logging, Printf, Dates

function generate_name(DOF, N, Nover, flux, periodic; L = 4π, mpicomm = MPI.COMM_WORLD, endtime = 200, FT = Float64)
    Ne = round(Int, DOF / (N+1))
    Nˣ = Ne
    Nʸ = Ne
    Nᶻ = Ne
    Lˣ = L 
    Lʸ = L
    Lᶻ = L
    # grid
    xrange = range(-Lˣ / 2; length = Nˣ + 1, stop = Lˣ / 2)
    yrange = range(-Lʸ / 2; length = Nʸ + 1, stop = Lʸ / 2)
    zrange = range(-Lᶻ / 2; length = Nᶻ + 1, stop = Lᶻ / 2)
    mpicomm = mpicomm
    brickrange = (xrange, yrange, zrange)
    topl = BrickTopology(
        mpicomm,
        brickrange,
        periodicity = (true, periodic, true),
        boundary = ((0, 0), (0, 0), (0,0)),
    )
    grid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = Array,
        polynomialorder = N + Nover,
    )
    cflgrid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = Array,
        polynomialorder = N,
    )
    
    Δx =  min_node_distance(cflgrid)
    cfl = 0.1 
    dt = cfl * Δx / √10

    # run
    timeend = FT(endtime) # s
    nout = round(Int, 2 / dt)
    dt = 2 / nout
    filename = "3Dflux_" * strip(string(flux), ['(', ')']) * "_p" * string(N) * "_N" * string(Ne) * "_Nover" * string(Nover) * "_periodicity_" * string(periodic)
    f = jldopen(filename * ".jld2", "a+")
    f["grid"] = grid
    close(f)
    return filename, Ne, dt
end

function just_generate_name(DOF, N, Nover, flux, periodic; L = 4π, mpicomm = MPI.COMM_WORLD, endtime = 200)
    Ne = round(Int, DOF / (N+1))
    Nˣ = Ne
    Nʸ = Ne
    Lˣ = L 
    Lʸ = L
    filename = "3Dflux_" * strip(string(flux), ['(', ')']) * "_p" * string(N) * "_N" * string(Ne) * "_Nover" * string(Nover) * "_periodicity_" * string(periodic)
    return filename
end

function nameprettifier(name)
    tmp = split(name, "_")
    dof = "polynomial order " * tmp[3][2] * " and " * tmp[4][2:end] * " elements"
    stringflux = ", a " * tmp[2] 
    intorder = ", and integration order 2(" * tmp[3][2] * "+" * tmp[5][end] * ")- 1," 
    return dof * stringflux * intorder
end
