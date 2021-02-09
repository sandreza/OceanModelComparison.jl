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

# for Rusanov flux:
# cfl of 0.333 for SPPRK22Heuns and no overintegration 0.5 failed, 56
# cfl of 0.25 for SPPRK22Heuns and  overintegration 0.5 failed, 133


DOFs = [256]
Ns = [1]
Novers = [1]
fluxes = [RusanovNumericalFlux()]
periodicity = [true]
cfls = [6,5]
timesteppers = [SSPRK22Heuns]
#  SSPRK22Heuns, LSRK54CarpenterKennedy, LSRK144NiegemannDiehlBusch, LSRKEulerMethod

for DOF in DOFs, N in Ns, Nover in Novers, flux in fluxes, periodic in periodicity, cfl in cfls, timestepper in timesteppers
    
L = 4π
endtime = 200.0

filename, Ne, dt = generate_name_3(DOF, N, Nover, flux, periodic, L = L, endtime = endtime, cflind = cfl)
println("Doing " * filename)
# simulation times
timeend = FT(endtime) # s
dt = FT(dt) # s
nout = round(Int, 2 / dt)
println("dt is ", dt)
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

run_bickley_jet(config, params; TimeStepper = timestepper, vtkpath = "")
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
