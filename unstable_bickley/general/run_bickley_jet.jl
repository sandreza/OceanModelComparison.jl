#!/usr/bin/env julia --project
# path = "/home/sandre/Desktop/Julia/OceanModelComparison.jl/unstable_bickley/general/"
# include(path * "run_bickley_jet.jl")
include("bickley_jet.jl")
include("convenience.jl")
ClimateMachine.init()

const FT = Float64

#################
# RUN THE TESTS #
#################
DOFs = [512]
Ns = [1, 2, 3, 4]
Novers = [0, 1]
fluxes = [RoeNumericalFlux(), RusanovNumericalFlux()]
periodicity = [false, true]

for DOF in DOFs, N in Ns, Nover in Novers, flux in fluxes, periodic in periodicity
    
L = 4π
endtime = 200.0

filename, Ne, dt = generate_name_2(DOF, N, Nover, flux, periodic, L = L, endtime = endtime)
println("Doing " * filename)
# simulation times
timeend = FT(endtime) # s
dt = FT(dt) # s
nout = round(Int, 2 / dt)
# Domain size
Lˣ = L  # m
Lʸ = L  # m
Nˣ = Ne
Nʸ = Ne
params = (; N, Nˣ, Nʸ, Lˣ, Lʸ, dt, nout, timeend)
if periodic
    config = Config(
        filename,
        params;
        numerical_flux_first_order = flux,
        Nover = Nover,
        periodicity = (true, true),
        boundary = ((0, 0), (0, 0)),
        boundary_conditons = (ClimateMachine.Ocean.OceanBC(
            Impenetrable(FreeSlip()),
            Insulating(),
        ),),
    )
else
    config = Config(
        filename,
        params;
        numerical_flux_first_order = flux,
        Nover = Nover,
        periodicity = (true, false),
        boundary = ((0, 0), (1,1)),
        boundary_conditons = (ClimateMachine.Ocean.OceanBC(
            Impenetrable(FreeSlip()),
            Insulating(),
        ),),
    )
end

tic = Base.time()

run_bickley_jet(config, params; TimeStepper = SSPRK22Heuns, vtkpath = "")
#  LSRK54CarpenterKennedy
toc = Base.time()
time = toc - tic
println(time)

f = jldopen(filename * ".jld2", "a+")
f["simulationtime"] = toc - tic
f["threads"] = Threads.nthreads()
f["arraytype"] = string(ClimateMachine.array_type())
close(f)
end
##
for  i in [1], j in [2,3]
    println(i,j)
end