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
xC, yC, zC = cellcenters(grid)
ϕ =  ScalarField(copy(x), gridhelper)
ϕ((0,0,0))
xnew = range(-2π, 2π, length = 3*43)
ynew = range(-2π, 2π, length = 3*43)
znew = range(0,0, length = 1 )
## comment, not a fair comparison, needs to be divided by polynomial order
## u = assemble(tmp).data[:, :, 1] about 10x slower


#=
nt = 100
ut = zeros(length(xnew), length(ynew), nt)
vt = zeros(length(xnew), length(ynew), nt)
ηt = zeros(length(xnew), length(ynew), nt)
ct = zeros(length(xnew), length(ynew), nt)
tic = time()
trange = collect(1:nt)
for i in 1:nt
    ϕ .= u_timeseries[i].data
    ut[:, :, i] .= view(ϕ(xnew, ynew, znew), :, :, 1)
    ϕ .= v_timeseries[i].data
    vt[:, :, i] .= view(ϕ(xnew, ynew, znew), :, :, 1)
    ϕ .= η_timeseries[i].data
    ηt[:, :, i] .= view(ϕ(xnew, ynew, znew), :, :, 1)
    ϕ .= c_timeseries[i].data
    ct[:, :, i] .= view(ϕ(xnew, ynew, znew), :, :, 1)
end
toc = time()
=#
##

resolution = (2358, 1364)
width = round(Int, resolution[1] / 5)
scene, layout = layoutscene(resolution = resolution)
lscene = layout[2, 2] = LScene(scene)
lscene2 = layout[2, 3] = LScene(scene)
lscene3 = layout[3, 2] = LScene(scene)
lscene4 = layout[3, 3] = LScene(scene)
layout[1,1] = LText(scene, "Menu", width = width, textsize = 50)

# Clim sliders
timerange = Int.(collect(range(1, 100, length = 100)))
time_slider = LSlider(scene, range = timerange, startvalue = 1)
time_node = time_slider.value

xrange = Int.(collect(range(8, 128, length = 128-8+1)))
x_slider = LSlider(scene, range = xrange, startvalue = 32)
x_node = x_slider.value



# u
xx = @lift(range(-2π, 2π, length = $x_node))
yy = @lift(range(-2π, 2π, length = $x_node))
ϕtmp =  @lift(ScalarField(u_timeseries[$time_node].data, gridhelper))
ϕu = @lift(($ϕtmp($xx, $yy, znew))[:,:,1])
urange = extrema(ut)
# ϕu = @lift(ut[:, :, $time_node])
heatmap!(lscene, xx, yy, ϕu, colormap = :balance, levels = 20, colorrange = urange, interpolate = true)
axis = scene.children[1][Axis]

# v 
vrange = extrema(vt)
#ϕv = @lift(vt[:, :, $time_node])
ϕtmp =  @lift(ScalarField(v_timeseries[$time_node].data, gridhelper))
ϕv = @lift(($ϕtmp($xx, $yy, znew))[:,:,1])
vrange = extrema(vt)
heatmap!(lscene2, xx, yy, ϕv, colormap = :balance, levels = 20, colorrange = vrange, interpolate = true)

# η 
ηrange = extrema(ηt)
# ϕη = @lift(ηt[:, :, $time_node])

ϕtmp =  @lift(ScalarField(η_timeseries[$time_node].data, gridhelper))
ϕη = @lift(($ϕtmp($xx, $yy, znew))[:,:,1])
# wierd bug
heatmap!(lscene3, xx, yy,  ϕη, colormap = :balance, levels = 20, colorrange = ηrange, interpolate = true)

# c 
crange = extrema(ct)
# ϕc = @lift(ct[:, :, $time_node])
ϕtmp =  @lift(ScalarField(c_timeseries[$time_node].data, gridhelper))
ϕc = @lift(($ϕtmp($xx, $yy, znew))[:,:,1])
heatmap!(lscene4, xx, yy, ϕc, colormap = :balance, levels = 20, colorrange = crange, interpolate = true)

# slider
layout[2:3, 1] = vgrid!(
    LText(scene, "time", width = nothing),
    time_slider,
    LText(scene, "coarseness", width = nothing),
    x_slider,
)
display(scene)