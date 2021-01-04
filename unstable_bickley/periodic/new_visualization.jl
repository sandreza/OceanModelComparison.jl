# Unstable Bickley jet

using Printf
using ClimateMachine
ClimateMachine.init()
using ClimateMachine.Ocean
using ClimateMachine.Ocean.Domains
using ClimateMachine.Ocean.Fields
using ClimateMachine.Mesh.Grids: DiscontinuousSpectralElementGrid
using ClimateMachine.GenericCallbacks: EveryXSimulationTime
using ClimateMachine.GenericCallbacks: EveryXSimulationSteps
using ClimateMachine.Ocean: current_step, Δt, current_time
using ClimateMachine.Ocean: JLD2Writer, OutputTimeSeries, write!
using CLIMAParameters: AbstractEarthParameterSet, Planet
##
# hack for current case

polynomialorders(::DiscontinuousSpectralElementGrid{T, dim, N}) where {T, dim, N} = Tuple([N for i in 1:dim])
include(pwd() * "/unstable_bickley/periodic/imperohooks.jl")

name = "climate_machine_unstable_bickley_jet_Ne43_Np2_ν0.0e+00_no_rotation_diffusive_cfl1.0e-02_reduction1.0e+04_smoothness_exponent2.0e+00"
filepath = name * ".jld2"

u_timeseries = OutputTimeSeries(:u, filepath);
v_timeseries = OutputTimeSeries(:v, filepath);
η_timeseries = OutputTimeSeries(:η, filepath);
c_timeseries = OutputTimeSeries(:θ, filepath);


## Note that simulation was actually 3D
grid = u_timeseries.grid
gridhelper = GridHelper(grid)
eh = gridhelper.element 
ih = gridhelper.interpolation     
x, y, z = coordinates(grid)
xC,yC,zC = cellcenters(grid)
ϕ =  ScalarField(copy(x), gridhelper)
ϕ((0,0,0))
xnew = range(-2π, 2π, length = 3*43)
ynew = range(-2π, 2π, length = 3*43)
znew = range(0,0, length = 1 )
## comment, not a fair comparison, needs to be divided by polynomial order
## u = assemble(tmp).data[:, :, 1] about 10x slower

scene = Scene()
for i in 1:100
    tmp = u_timeseries[i]
    ϕ .= tmp.data
    ϕu = ϕ(xnew, ynew, znew)
    contourf!(scene, xnew, ynew, ϕu[:,:,1], colormap = :balance, levels = 20)
end
i = 90
tmp = c_timeseries[i]
ϕ .= tmp.data
ϕu = ϕ(xnew, ynew, znew)
contourf!(scene, xnew, ynew, ϕu[:,:,1], colormap = :balance, levels = 20)