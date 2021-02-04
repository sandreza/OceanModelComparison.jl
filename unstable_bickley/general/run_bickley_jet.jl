#!/usr/bin/env julia --project

include("bickley_jet.jl")
include("convenience.jl")
ClimateMachine.init()

const FT = Float64

#################
# RUN THE TESTS #
#################

DOF = 32 # 32, 128, 512
N = 3 # 1,2,3,4
Nover = 1 # 1
flux = RoeNumericalFlux() # RusanovNumericalFlux()
periodic = true #false
L = 4π
endtime = 200.0

filename, Ne, dt = generate_name(DOF, N, Nover, flux, periodic, L = L, endtime = endtime)

# simulation times
timeend = FT(endtime) # s
dt = FT(dt) # s
nout = Int(100)
# Domain size
Lˣ = L  # m
Lʸ = L  # m
Nˣ = Ne
Nʸ = Ne
params = (; N, Nˣ, Nʸ, Lˣ, Lʸ, dt, nout, timeend)
##
if periodic
    config = Config(
        filename,
        params;
        numerical_flux_first_order = flux,
        Nover = 0,
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
        Nover = 0,
        periodicity = (true, false),
        boundary = ((0, 0), (1,1)),
        boundary_conditons = (ClimateMachine.Ocean.OceanBC(
            Impenetrable(FreeSlip()),
            Insulating(),
        ),),
    )
end

tic = Base.time()

run_bickley_jet(config, params; TimeStepper = SSPRK22Heuns)
#  LSRK54CarpenterKennedy
toc = Base.time()
time = toc - tic
println(time)

f = jldopen(filename * ".jld2", "a+")
f["simulationtime"] = toc - tic
f["threads"] = Threads.nthreads()
f["arraytype"] = string(ClimateMachine.array_type())
close(f)
