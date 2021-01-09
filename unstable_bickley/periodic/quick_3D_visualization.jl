# Assumes new visualization.jl has been round
# and vizinanigans
using Random, Statistics
x, y, z = coordinates(gridu)
xC, yC, zC = cellcenters(gridu)
ϕ =  ScalarField(copy(x), gridhelper)

xnew = range(-2π, 2π, length = 3*43)
ynew = range(-2π, 2π, length = 3*43)
znew = range(0,1, length = 10 )

i = 90
ϕ .=  u_timeseries[i].data
ϕu = ϕ(xnew, ynew, znew)

ϕ .=  v_timeseries[i].data
ϕv = ϕ(xnew, ynew, znew)

ϕ .=  η_timeseries[i].data
ϕη = ϕ(xnew, ynew, znew)

ϕ .=  c_timeseries[i].data
ϕc = ϕ(xnew, ynew, znew)
states = [ϕu, ϕv, ϕη, ϕc]
statenames = ["u", "v", "η", "c"]
statistics = true
title = "Field at t = " * string(2*i) * " :"
scene = visualize(states, statenames = statenames, aspect = (1,1,0.1), statistics = statistics, title = title);
##
scene = volumeslice(states, statenames = statenames, aspect = (1,1,0.1), statistics = statistics)
##
record_interaction = false
seconds = 20
fps = 10
frames = round(Int, fps * seconds )
if record_interaction
record(scene, pwd() * "/3dbickley.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end
end