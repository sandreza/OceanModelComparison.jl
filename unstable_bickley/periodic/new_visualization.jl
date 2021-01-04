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
# scene.px_area to get current resolution
resolution = (2062, 1430)
width = round(Int, resolution[1] / 5)
scene, layout = layoutscene(resolution = resolution)
lscene = layout[2:3, 2:3] = LScene(scene, title = "u")
lscene2 = layout[2:3, 4:5] = LScene(scene)
lscene3 = layout[5:6, 2:3] = LScene(scene)
lscene4 = layout[5:6, 4:5] = LScene(scene)
layout[1,1] = LText(scene, "Menu", width = width, textsize = 50)

padding = (width/2, 0 , 0, 0)
layout[1,2:3] = LText(scene, "u", width = width, textsize = 50, padding = padding)
layout[4,2:3] = LText(scene, "η", width = width, textsize = 50, padding = padding)
layout[1,4:5] = LText(scene, "v", width = width, textsize = 50, padding = padding)
layout[4,4:5] = LText(scene, "c", width = width, textsize = 50, padding = padding)

# Clim sliders
timerange = Int.(collect(range(1, 100, length = 100)))
time_slider = LSlider(scene, range = timerange, startvalue = 1)
time_node = time_slider.value

xrange = Int.(collect(range(8, 128, length = 128-8+1)))
yrange = xrange
x_slider = LSlider(scene, range = xrange, startvalue = 32)
x_node = x_slider.value

y_slider = LSlider(scene, range = yrange, startvalue = 32)
y_node = y_slider.value

xx = @lift(range(-2π, 2π, length = $x_node))
yy = @lift(range(-2π, 2π, length = $y_node))

# Menu
interpolationlabels = ["contour", "heatmap"]
interpolationchoices = [true, false]
interpolationnode = Node(interpolationchoices[1])


# u
ϕtmp =  @lift(ScalarField(u_timeseries[$time_node].data, gridhelper))
ϕu = @lift(($ϕtmp($xx, $yy, znew))[:,:,1])
urange = extrema(ut)
# ϕu = @lift(ut[:, :, $time_node])
heatmap1 = heatmap!(lscene, xx, yy, ϕu, colormap = :balance, levels = 20, colorrange = urange, interpolate = interpolationnode)

# color bar
cbar = LColorbar(scene, heatmap1)
cbar.width = Relative(1/3)
cbar.height = Relative(2/3)
cbar.ticksvisible = false
cbar.ticklabelsize = 0
cbar.halign = :left
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

#
interpolationmenu = LMenu(scene, options = zip(interpolationlabels, interpolationchoices))
on(interpolationmenu.selection) do s
    interpolationnode[] = s
    # hack
    heatmap!(lscene, xx, yy, ϕu, colormap = :balance, colorrange = urange, interpolate = s)
    heatmap!(lscene2, xx, yy, ϕv, colormap = :balance, colorrange = vrange, interpolate = s)
    heatmap!(lscene3, xx, yy, ϕη, colormap = :balance, colorrange = ηrange, interpolate = s)
    heatmap!(lscene4, xx, yy, ϕc, colormap = :balance, colorrange = crange, interpolate = s)
end

timetext = @lift("time, t = " * @sprintf("%0.1f", 2 * ($time_node - 1)))
xtext = @lift("x interpolation = " * string($x_node) * " points")
ytext = @lift("y interpolation = " * string($y_node) * " points")
# slider padding = (left, right, bottom, top)
# last entry moves title up and down
# second to last pushes things away
layout[2:6, 1] = vgrid!(
    LText(scene, "plotting options", width = width, textsize = 30, padding = (0,0, 10, 0)),
    interpolationmenu,
    LText(scene, timetext, width = width, textsize = 30),
    time_slider,
    LText(scene, xtext, width = width, textsize = 30),
    x_slider,
    LText(scene, ytext, width = width, textsize = 30),
    y_slider,
    LText(scene, "ColorBar", width = width, textsize = 50, padding = (0, 0, 0, 100)),
    cbar,
)
display(scene)