states = vizstates
statenames = string.(1:length(states))
units = ["" for i in eachindex(states)]
aspect = (1, 1, 1)
resolution = (2678, 1030)
statistics = true
title = "Volume plot of "
bins = 300
statlabelsize = (20,20)

scene, layout = layoutscene(resolution = resolution)
volumescene = layout[2:4, 2:4] = LScene(scene)
menuwidth = round(Int, 350)
layout[1,1] = Label(scene, "Menu", width = menuwidth, textsize = 50)

slice_slider = Slider(scene, range = range(0, 1, length = 101), startvalue = 0.0)
slice_node = slice_slider.value

directionindex = [1, 2, 3]
directionnames = ["x-slice", "y-slice", "z-slice"]
directionnode = Node(directionindex[1])

stateindex = collect(1:length(states))
stateslider = Slider(scene, range = Int.(collect(1:length(states))), startvalue = 1)
statenode = stateslider.value

layout[1, 2:4] = Label(scene, @lift(title * statenames[$statenode]), textsize = 50)

colorchoices = [:balance, :thermal, :dense, :deep, :curl, :thermometer]
colornode = Node(colorchoices[1])

state = @lift(states[$statenode])
statename = @lift(statenames[$statenode])
unit = @lift(units[$statenode])
nx = @lift(size($state)[1])
ny = @lift(size($state)[2])
nz = @lift(size($state)[3])
nr = @lift([$nx, $ny, $nz])

nslider = 100
xrange = range(0.00, aspect[1], length = nslider)
yrange = range(0.00, aspect[2], length = nslider)
zrange = range(0.00, aspect[3], length = nslider)
constx = collect(reshape(xrange, (nslider,1,1)))
consty = collect(reshape(yrange, (1,nslider,1)))
constz = collect(reshape(zrange, (1,1,nslider)))
matx = zeros(nslider,nslider,nslider)
maty = zeros(nslider,nslider,nslider)
matz = zeros(nslider,nslider,nslider)
matx .= constx
maty .= consty
matz .= constz
sliceconst = [matx, maty, matz]
planeslice = @lift(sliceconst[$directionnode])

upperclim_slider = Slider(scene, range = range(0, 1, length = 101), startvalue = 0.99)
upperclim_node = upperclim_slider.value
lowerclim_slider = Slider(scene, range = range(0, 1, length = 101), startvalue = 0.01)
lowerclim_node = lowerclim_slider.value

clims = @lift((quantile($state[:], $lowerclim_node) , quantile($state[:], $upperclim_node)))

volume!(volumescene, 0..aspect[1], 0..aspect[2], 0..aspect[3], state, overdraw = false, colorrange = clims, colormap = @lift(to_colormap($colornode)))


alpha_slider = Slider(scene, range = range(0, 1, length = 101), startvalue = 0.5)
alphanode = alpha_slider.value

slicecolormap = @lift(cgrad(:viridis, alpha = $alphanode))
v = volume!(volumescene, 0..aspect[1], 0..aspect[2], 0..aspect[3],
        planeslice, algorithm = :iso, isorange = 0.005, 
        isovalue = @lift($slice_node * aspect[$directionnode]),
        transparency = true, overdraw = false, visible = true, 
        colormap = slicecolormap, colorrange = [-1,0] )

# Volume histogram

layout[3, 1] = Label(scene, "Statistics", textsize = 50)
hscene = layout[4, 1] = Axis(scene, xlabel = @lift(statenames[$statenode] * " " * units[$statenode]), 
                    xlabelcolor = :black, ylabel = "pdf", 
                    ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
                    xticklabelsize = statlabelsize[1], yticklabelsize = statlabelsize[2],
                    xtickcolor = :black, ytickcolor = :black,
                    xticklabelcolor  = :black, yticklabelcolor = :black)

histogram_node = @lift(histogram($state, bins = bins))
vxs = @lift($histogram_node[1])
vys = @lift($histogram_node[2])
pdf = GLMakie.AbstractPlotting.barplot!(hscene, vxs, vys, color = :red, 
                strokecolor = :red, 
                strokewidth = 1)

@lift(GLMakie.AbstractPlotting.xlims!(hscene, extrema($vxs))) 
@lift(GLMakie.AbstractPlotting.ylims!(hscene, extrema($vys)))
vlines!(hscene, @lift($clims[1]), color = :black, linewidth = menuwidth / 100)
vlines!(hscene, @lift($clims[2]), color = :black, linewidth = menuwidth / 100)


# Slice
sliceupperclim_slider = Slider(scene, range = range(0, 1, length = 101), startvalue = 0.99)
sliceupperclim_node = sliceupperclim_slider.value
slicelowerclim_slider = Slider(scene, range = range(0, 1, length = 101), startvalue = 0.01)
slicelowerclim_node = slicelowerclim_slider.value


slicexaxislabel = @lift(["y", "x", "x"][$directionnode])
sliceyaxislabel = @lift(["z", "z", "y"][$directionnode])

slicexaxis = @lift([[1,$ny], [1, $nx], [1,$nx]][$directionnode])
sliceyaxis = @lift([[1,$nz], [1, $nz], [1,$ny]][$directionnode])

slicescene = layout[2:4, 5:6] = Axis(scene, xlabel = slicexaxislabel, ylabel = sliceyaxislabel)

sliced_state1 = @lift( $state[round(Int, 1 + $slice_node * (size($state)[1]-1)), 1:size($state)[2], 1:size($state)[3]])
sliced_state2 = @lift( $state[1:size($state)[1], round(Int, 1 + $slice_node * (size($state)[2]-1)), 1:size($state)[3]])
sliced_state3 = @lift( $state[1:size($state)[1], 1:size($state)[2], round(Int, 1 + $slice_node * (size($state)[3]-1))]) 
sliced_states = @lift([$sliced_state1, $sliced_state2, $sliced_state3])
sliced_state = @lift($sliced_states[$directionnode]) 

oclims = @lift((quantile($sliced_state[:], $slicelowerclim_node) , quantile($sliced_state[:], $sliceupperclim_node)))
slicecolormapnode = @lift($oclims[1] < $oclims[2] ? $colornode : $colornode )
sliceclims = @lift($oclims[1] != $oclims[2] ? (minimum($oclims), maximum($oclims)) : (minimum($oclims)-1, maximum($oclims)+1))

heatmap1 = heatmap!(slicescene, slicexaxis, sliceyaxis, sliced_state, interpolate = true, colormap = slicecolormapnode, colorrange = sliceclims)

# Colorbar
newlabel = @lift($statename * " " * $unit)
cbar = Colorbar(scene, heatmap1, label = newlabel)
cbar.width = Relative(1/3)
# cbar.height = Relative(5/6)
cbar.halign = :left
# cbar.flipaxisposition = true
# cbar.labelpadding = -250
cbar.labelsize = 50

@lift(GLMakie.AbstractPlotting.xlims!(slicescene, extrema($slicexaxis))) 
@lift(GLMakie.AbstractPlotting.ylims!(slicescene, extrema($sliceyaxis)))

sliceindex = @lift([round(Int, 1 + $slice_node * ($nx-1)), round(Int, 1 + $slice_node * ($ny-1)), round(Int, 1 + $slice_node * ($nz-1))][$directionnode])
slicestring = @lift(directionnames[$directionnode] * " of " * statenames[$statenode] ) 
layout[1, 5:6] = Label(scene, slicestring, textsize = 50)

#=
axis = scene.children[1][OldAxis] 
axis[:names][:axisnames] = ("↓", "↓ ", "↓ ")
axis[:names][:align] = ((:left, :center), (:right, :center), (:right, :center))
axis[:names][:textsize] = (50.0, 50.0, 50.0)
axis[:ticks][:textsize] = (00.0, 00.0, 00.0)
=#

# Menus
colormenu = Menu(scene, options = zip(colorchoices, colorchoices))
on(colormenu.selection) do s
    colornode[] = s
end


# Slice Statistics
layout[1, 7] = Label(scene, "Slice Menu", width = menuwidth, textsize = 50)
layout[3, 7] = Label(scene, "Slice Statistics", textsize = 50)
hslicescene = layout[4, 7] = Axis(scene, xlabel = @lift(statenames[$statenode] * " " * units[$statenode]), 
                    xlabelcolor = :black, ylabel = "pdf", 
                    ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
                    xticklabelsize = statlabelsize[1], yticklabelsize = statlabelsize[2],
                    xtickcolor = :black, ytickcolor = :black,
                    xticklabelcolor  = :black, yticklabelcolor = :black)

slicehistogram_node = @lift(histogram($sliced_state, bins = bins))
xs = @lift($slicehistogram_node[1])
ys = @lift($slicehistogram_node[2])
pdf = GLMakie.AbstractPlotting.barplot!(hslicescene, xs, ys, color = :blue, 
                strokecolor = :blue, 
                strokewidth = 1)

@lift(GLMakie.AbstractPlotting.xlims!(hslicescene, extrema($xs))) 
@lift(GLMakie.AbstractPlotting.ylims!(hslicescene, extrema($ys)))
vlines!(hslicescene, @lift($sliceclims[1]), color = :black, linewidth = menuwidth / 100)
vlines!(hslicescene, @lift($sliceclims[2]), color = :black, linewidth = menuwidth / 100)

interpolationnames = ["contour", "heatmap"]
interpolationchoices = [true, false]
interpolationnode = Node(interpolationchoices[1])
interpolationmenu = Menu(scene, options = zip(interpolationnames, interpolationchoices))

on(interpolationmenu.selection) do s
    interpolationnode[] = s
    # hack
    heatmap!(slicescene, slicexaxis, sliceyaxis, sliced_state, interpolate = s, colormap = slicecolormapnode, colorrange = sliceclims)
end

directionmenu = Menu(scene, options = zip(directionnames, directionindex))

on(directionmenu.selection) do s
    directionnode[] = s
end

statesliderstring = @lift("State at t = " * string($statenode))
slicemenustring = @lift(directionnames[$directionnode] * " at index "  * string(round(Int, 1 + $slice_node * ($nr[$directionnode]-1)))) 
lowerclim_string = @lift("quantile = " *  @sprintf("%0.2f", $lowerclim_node) * ", value = " * @sprintf("%0.1e", $clims[1]))
upperclim_string = @lift("quantile = " *  @sprintf("%0.2f", $upperclim_node) * ", value = " * @sprintf("%0.1e", $clims[2]))
alphastring = @lift("Slice alpha = " * @sprintf("%0.2f", $alphanode))
layout[2, 1] = vgrid!(
    Label(scene, statesliderstring, width = nothing),
    stateslider,
    Label(scene, "Color", width = nothing),
    colormenu,
    Label(scene, "Slice Direction", width = nothing),
    directionmenu,
    Label(scene, alphastring, width = nothing),
    alpha_slider,
    Label(scene, slicemenustring, width = nothing),
    slice_slider,
    Label(scene, lowerclim_string, width = nothing),
    lowerclim_slider,
    Label(scene, upperclim_string, width = nothing),
    upperclim_slider,
)

slicelowerclim_string = @lift("quantile = " *  @sprintf("%0.2f", $slicelowerclim_node) * ", value = " * @sprintf("%0.1e", $sliceclims[1]))
sliceupperclim_string = @lift("quantile = " *  @sprintf("%0.2f", $sliceupperclim_node) * ", value = " * @sprintf("%0.1e", $sliceclims[2]))

layout[2,7] = vgrid!(
    Label(scene, "Contour Plot Type", width = nothing), 
    interpolationmenu,
    Label(scene, slicelowerclim_string, width = nothing),
    slicelowerclim_slider,
    Label(scene, sliceupperclim_string, width = nothing), 
    sliceupperclim_slider,
    cbar,
)
display(scene)
##
record(scene, "64dof.mp4", 1:100, framerate=10) do n
    statenode[] = n
    θ = -0.005 * 2π
    alphanode[] = max(0.5 - n/100 ,0)
    rotate_cam!(scene.children[1], (θ, 0, 0))
end