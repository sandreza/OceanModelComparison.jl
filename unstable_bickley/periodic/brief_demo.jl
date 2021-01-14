# include dependencies
include(pwd() * "/unstable_bickley/periodic/dependencies.jl")

# filenames / load
names = [
"climate_machine_unstable_bickley_jet_Ne43_Np2_ν0.0e+00_no_rotation_diffusive_cfl1.0e-02_reduction1.0e+04_smoothness_exponent2.0e+00",
"climate_machine_unstable_bickley_jet_Ne8_Np3_ν0.0e+00_no_rotation_exasim_comparison",
"climate_machine_unstable_bickley_jet_Ne11_Np2_ν0.0e+00_no_rotation_diffusive_cfl1.0e-03_reduction1.0e+02_smoothness_exponent1.0e+00",
"climate_machine_unstable_bickley_jet_Ne73_Np6_ν0.0e+00_no_rotation_diffusive_cfl1.0e+00_reduction1.0e+08_smoothness_exponent1.0e+00",
]
name = names[2]
filepath = name * ".jld2"

u_timeseries = OutputTimeSeries(:u, filepath);
v_timeseries = OutputTimeSeries(:v, filepath);
η_timeseries = OutputTimeSeries(:η, filepath);
c_timeseries = OutputTimeSeries(:θ, filepath);

## Add ImperoHooks utils
gridu = u_timeseries.grid
gridhelper = GridHelper(gridu)
eh = gridhelper.element 
ih = gridhelper.interpolation     
x, y, z = coordinates(gridu)
xC, yC, zC = cellcenters(gridu)
# Now we can evaluate an MPIStateArray at an arbitrary point
ϕ =  ScalarField(copy(x), gridhelper)
ϕ(0,0,0)
xnew = range(-2π, 2π, length = 3*43)
ynew = range(-2π, 2π, length = 3*43)
znew = range(0, 1, length = 10 )
# Or even a range of points
newϕ =  ϕ(xnew, ynew, znew)

# Futhermore we can visualize it, for example
# set data 
nt = length(u_timeseries)
ϕ .= c_timeseries[nt].data .* u_timeseries[nt].data
states = [ϕ(xnew,ynew,znew)]
statenames = ["uη"]
volumeslice(states, statenames = statenames, statistics = true, aspect = (1,1,0.1))
##
# interpolate the timeseries
xnew = range(-2π, 2π, length = 3*43)
ynew = range(-2π, 2π, length = 3*43)
znew = range(1, 1, length = 1 )

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
states = [ut, vt, ηt, ct]
statenames = ["u", "v", "η", "c"]

scene = volumeslice(states, statenames = statenames, bins = 30)
