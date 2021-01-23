# Unstable Bickley jet

using Printf, JLD2
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

names = [
"climate_machine_unstable_bickley_jet_Ne43_Np2_ν0.0e+00_no_rotation_diffusive_cfl1.0e-02_reduction1.0e+04_smoothness_exponent2.0e+00",
"climate_machine_unstable_bickley_jet_Ne8_Np3_ν0.0e+00_no_rotation_test43",
"climate_machine_unstable_bickley_jet_Ne73_Np6_ν0.0e+00_no_rotation_diffusive_cfl1.0e+00_reduction1.0e+08_smoothness_exponent1.0e+00",
"climate_machine_unstable_bickley_jet_Ne8_Np3_ν0.0e+00_no_rotation_no_derivative",
"climate_machine_unstable_bickley_jet_Ne8_Np3_ν0.0e+00_no_rotation_derivative",
"climate_machine_unstable_bickley_jet_Ne8_Np3_ν0.0e+00_no_rotation_more_levels",
"climate_machine_unstable_bickley_jet_Ne8_Np3_ν0.0e+00_no_rotation_no_derivative_3",
"climate_machine_unstable_bickley_jet_Ne8_Np3_ν0.0e+00_no_rotation_no_derivative_4",
"climate_machine_unstable_bickley_jet_Ne21_Np2_ν0.0e+00_no_rotation_derivative_2",
"climate_machine_unstable_bickley_jet_Ne8_Np3_ν0.0e+00_no_rotation_exasim_comparison",
]
name = names[end]
filepath = name * ".jld2"
jldfile = jldopen(filepath)
u_timeseries = OutputTimeSeries(:u, filepath);
v_timeseries = OutputTimeSeries(:v, filepath);
η_timeseries = OutputTimeSeries(:η, filepath);
c_timeseries = OutputTimeSeries(:θ, filepath);
extrema(c_timeseries[100].data)

# reexport grid data
dg_grid = jldfile["grid"]
x, y, z = coordinates(dg_grid)
xC, yC, zC = cellcenters(dg_grid)
Ω = extrema.((x,y,z))
topologystring = string(typeof(dg_grid.topology))
polynomialorder = polynomialorders(dg_grid)
addup(xC, tol) = sum(abs.(xC[1] .- xC) .≤ tol)
ne = length(xC)
ex = round(Int64, ne / addup(xC, 10^4 * eps(maximum(abs.(x)))))
ey = round(Int64, ne / addup(yC, 10^4 * eps(maximum(abs.(y)))))
ez = round(Int64, ne / addup(zC, 10^4 * eps(maximum(abs.(z)))))
elements = (ex,ey,ez)
velocity_u = [Array(u_timeseries[i].data) for i in 1:2:length(u_timeseries)]
velocity_v = [Array(v_timeseries[i].data) for i in 1:2:length(u_timeseries)]
scalar_η = [Array(η_timeseries[i].data) for i in 1:2:length(u_timeseries)]
scalar_c = [Array(c_timeseries[i].data) for i in 1:2:length(u_timeseries)]


@save "example.jld2" Ω polynomialorder elements topologystring velocity_u velocity_v scalar_η scalar_c

##
## Note that simulation was actually 3D
gridu = u_timeseries.grid
gridhelper = GridHelper(gridu)
eh = gridhelper.element 
ih = gridhelper.interpolation     
x, y, z = coordinates(gridu)
xC, yC, zC = cellcenters(gridu)
ϕ =  ScalarField(copy(x), gridhelper)
ϕ((0,0,0))
xnew = range(-2π, 2π, length = 3*43)
ynew = range(-2π, 2π, length = 3*43)
znew = range(1,1, length = 1 )
ϕ(xnew, ynew, znew)
## comment, not a fair comparison, needs to be divided by polynomial order
## u = assemble(tmp).data[:, :, 1] about 10x slower
ti = 50
field = u_timeseries[ti]
u1 = assemble(field).data[:,:,1]
u2 = assemble(field).data[:,:,2]
u3 = assemble(field).data[:,:,3]

nt = length(u_timeseries) - 1
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
println(toc - tic)

##
i = nt - 1
states = [ut, vt, ηt, ct]
statenames = ["u", "v", "η", "c"]

scene = volumeslice(states, statenames = statenames, bins = 30)

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
timerange = Int.(collect(range(1, nt, length = nt)))
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
heatmap1 = AbstractPlotting.heatmap!(lscene, xx, yy, ϕu, colormap = :balance, levels = 20, colorrange = urange, interpolate = interpolationnode)

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
AbstractPlotting.heatmap!(lscene2, xx, yy, ϕv, colormap = :balance, levels = 20, colorrange = vrange, interpolate = true)


# η 
# ηrange = extrema(ut) ./ 100000
ϕη = @lift(ηt[:, :, $time_node])

# ϕtmp =  @lift(ScalarField(u_timeseries[$time_node].data, gridhelper))

# ϕη = @lift(($ϕtmp($xx, $yy, znew))[:,:,1] - ($ϕtmp($xx, $yy, range(0.5,0.5, length = 1 )))[:,:,1])
ηrange = @lift(extrema($ϕη))
# wierd bug
AbstractPlotting.heatmap!(lscene3, xx, yy,  ϕη, colormap = :balance, levels = 20, colorrange = ηrange, interpolate = true)

# c 
crange = extrema(ct)
# ϕc = @lift(ct[:, :, $time_node])
ϕtmp =  @lift(ScalarField(c_timeseries[$time_node].data, gridhelper))
ϕc = @lift(($ϕtmp($xx, $yy, znew))[:,:,1])
AbstractPlotting.heatmap!(lscene4, xx, yy, ϕc, colormap = :balance, levels = 20, colorrange = crange, interpolate = true)

#
interpolationmenu = LMenu(scene, options = zip(interpolationlabels, interpolationchoices))
on(interpolationmenu.selection) do s
    interpolationnode[] = s
    # hack
    AbstractPlotting.heatmap!(lscene, xx, yy, ϕu, colormap = :balance, colorrange = urange, interpolate = s)
    AbstractPlotting.heatmap!(lscene2, xx, yy, ϕv, colormap = :balance, colorrange = vrange, interpolate = s)
    AbstractPlotting.heatmap!(lscene3, xx, yy, ϕη, colormap = :balance, colorrange = ηrange, interpolate = s)
    AbstractPlotting.heatmap!(lscene4, xx, yy, ϕc, colormap = :balance, colorrange = crange, interpolate = s)
end

timetext = @lift("time, t = " * @sprintf("%0.1f", 2 * ($time_node - 1)/99*100))
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
    LText(scene, "state min/max over simulation", width = width, textsize = 20, padding = (0, 0, 0, 00)),
    cbar,
)
display(scene)
##
record_interaction = false
seconds = 20
fps = 10
frames = round(Int, fps * seconds )
if record_interaction
record(scene, pwd() * "/test.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end
end