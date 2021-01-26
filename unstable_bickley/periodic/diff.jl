resolution = (1956+400, 852)
scene, layout = layoutscene(resolution = resolution )
lscene = layout[2:4,2:4] = LScene(scene)
lscene2 = layout[2:4, 5:7] = LScene(scene)
lscene3 = layout[2:4, 8:10] = LScene(scene)
layout[1, 2:4] = LText(scene, "Exasim ρθ", textsize = 50)
layout[1, 5:7] = LText(scene, "OI ClimateMachine ρθ", textsize = 50)
layout[1, 8:10] = LText(scene, "Difference ρθ", textsize = 50)
layout[1,1] = LText(scene, "DOF = $(DOF)^2", textsize = 50)
time_slider = LSlider(scene, range = Int.(range(1, 100, length=100)), startvalue = 1)
time_node = time_slider.value

# Exasim
estate1 = @lift(estates[1][:,:,$time_node])
estate2 = @lift(estates[2][:,:,$time_node])
estate3 = @lift(estates[3][:,:,$time_node])
estate4 = @lift(estates[4][:,:,$time_node])
clims = (-1,1)
heatmap1 = heatmap!(lscene, 0..1, 1..2, estate1, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene, 1..2, 1..2, estate2, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene, 0..1, 0..1, estate3, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene, 1..2, 0..1, estate4, colorrange = clims, colormap = to_colormap(:balance))

# Overintegration
istate1 = @lift(istates[1][:,:,$time_node])
istate2 = @lift(istates[2][:,:,$time_node])
istate3 = @lift(istates[3][:,:,$time_node])
istate4 = @lift(istates[4][:,:,$time_node])
clims = (-1,1)
heatmap1 = heatmap!(lscene2, 0..1, 1..2, istate1, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene2, 1..2, 1..2, istate2, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene2, 0..1, 0..1, istate3, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene2, 1..2, 0..1, istate4, colorrange = clims, colormap = to_colormap(:balance))


# Difference
diffstate1 = @lift($estate1 - $istate1)
diffstate2 = @lift($estate2 - $istate2)
diffstate3 = @lift($estate3 - $istate3)
diffstate4 = @lift($estate4 - $istate4)

clims = (-1,1)
heatmap1 = heatmap!(lscene3, 0..1, 1..2, diffstate1, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene3, 1..2, 1..2, diffstate2, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene3, 0..1, 0..1, diffstate3, colorrange = clims, colormap = to_colormap(:balance))
heatmap!(lscene3, 1..2, 0..1, diffstate4, colorrange = clims, colormap = to_colormap(:balance))

cbar = LColorbar(scene, heatmap1)
cbar.height = Relative(1/3)
cbar.width = Relative(1/3)
cbar.halign = :left
cbar.labelsize = 50

slidertext = @lift("Time t = " * string(2 *  $time_node))
p1norm = @lift(norm($diffstate1)/norm($istate1))
p2norm = @lift(norm($diffstate2)/norm($istate2))
p3norm = @lift(norm($diffstate3)/norm($istate3))
p4norm = @lift(norm($diffstate4)/norm($istate4))

relerror1 = @lift("relative difference p1 = " *  @sprintf("%0.2e", $p1norm) )
relerror2 = @lift("relative difference p2 = " *  @sprintf("%0.2e", $p2norm) )
relerror3 = @lift("relative difference p3 = " *  @sprintf("%0.2e", $p3norm) )
relerror4 = @lift("relative difference p4 = " *  @sprintf("%0.2e", $p4norm) )
layout[2:4, 1] = vgrid!(
    LText(scene, slidertext, width = nothing),
    time_slider,
    cbar,
    LText(scene, "top left p=1, top right p=2"),
    LText(scene, "bottom left p=3, bottom right p=4"),
    LText(scene, relerror1),
    LText(scene, relerror2),
    LText(scene, relerror3),
    LText(scene, relerror4)
)
display(scene)

##
record_interaction = true
seconds = 20
fps = 10
frames = round(Int, fps * seconds )
if record_interaction
record(scene, pwd() * "/comparison.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end
end