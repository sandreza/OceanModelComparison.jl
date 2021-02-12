include(pwd() * "/three_dimensional/convenience.jl")
include(pwd() * "/unstable_bickley/periodic/imperohooks.jl")
include(pwd() * "/unstable_bickley/periodic/vizinanigans.jl")

DOFs = [64]
Ns = [2]
Novers = [1]
fluxes = [RoeNumericalFlux()]
periodicity = [true] # [true, false]

states = []
namelist = []
for DOF in DOFs, Nover in Novers, flux in fluxes, periodic in periodicity, N in Ns
    name = just_generate_name(DOF, N, Nover, flux, periodic)
    push!(namelist, name)
    print(name)
    println(" ")
end
##
for name in namelist
    f = jldopen(name * ".jld2")
    simtime = f["simulationtime"]
    threadnum = f["threads"]
    array = f["arraytype"]
    if array == "Array" 
        archstring = " on the CPU with " * string(threadnum) * " threads"
    else
        archstring = " on the GPU"
    end
    prettyname = nameprettifier(name)
    println("The simulation time was " * @sprintf("%0.2f", simtime) * " seconds for " * prettyname * archstring)
    println("------------------------------")
    # get old grid
    newgrid = f["grid"]
    gridhelper = GridHelper(newgrid)  
    x, y, z = coordinates(newgrid)
    ϕ =  ScalarField(copy(x), gridhelper)
    # new grid
    newx = range(-2π, 2π, length = DOFs[1]* 2)
    newy = range(-2π, 2π, length = DOFs[1]* 2)
    newz = range(-2π, 2π, length = DOFs[1]* 2)
    ρθ = zeros(length(newx), length(newy), length(newz), 100)
    # interpolate
    M = view(newgrid.vgeo, :, newgrid.Mid, :)
    for i in 1:100
        Q = f[string(i)]
        ϕ .= Q[:,5,:]
        ρθ[:,:,:,i] = ϕ(newx, newy, newz, threads = true)
    end
    println("the mean value is ", sum(M .* ϕ.data))
    val = sum(M .* ϕ.data)/sum(M .* abs.(ϕ.data))
    if val ≥ eps(1e0)
        println(name * " not conserved to machine precision at " * string(100))
        println("the relative value is ", val)
        println("the value is ", sum(M .* ϕ.data))
    end
    push!(states, ρθ)
    close(f)
end
##
tmpstates = copy(states)
vizstates = [states[1][:,:,:,i] for i in 1:100]
scene = volumeslicetime(vizstates, aspect = (1,1,1))


##* 2
seconds = 10
fps = 10
frames = round(Int, fps * seconds )
record(scene, pwd() * "/3Dbickley.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end
##
record(scene, "64dof.mp4", 1:100, framerate=10) do n
    time_node[] = n
    θ = -0.005 * 2π
    rotate_cam!(scene.children[1], (θ, 0, 0))
end