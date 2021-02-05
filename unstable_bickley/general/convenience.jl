function generate_name(DOF, N, Nover, flux, periodic; L = 4π, mpicomm = MPI.COMM_WORLD, endtime = 200)
    Ne = round(Int, DOF / (N+1))
    Nˣ = Ne
    Nʸ = Ne
    Lˣ = L 
    Lʸ = L
    # grid
    xrange = range(-Lˣ / 2; length = Nˣ + 1, stop = Lˣ / 2)
    yrange = range(-Lʸ / 2; length = Nʸ + 1, stop = Lʸ / 2)
    mpicomm = mpicomm
    brickrange = (xrange, yrange)
    topl = BrickTopology(
        mpicomm,
        brickrange,
        periodicity = (true, periodic),
        boundary = ((0, 0), (0, 0)),
    )
    grid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = Array,
        polynomialorder = N + Nover,
    )
    
    Δx =  min_node_distance(grid)
    cfl = 0.15
    dt = cfl * Δx / √10
    # run
    timeend = FT(endtime) # s
    nout = round(Int, 2 / dt)
    dt = 2 / nout
    filename = "flux_" * strip(string(flux), ['(', ')']) * "_p" * string(N) * "_N" * string(Ne) * "_Nover" * string(Nover) * "_periodicity_" * string(periodic)
    f = jldopen(filename * ".jld2", "a+")
    f["grid"] = grid
    close(f)
    return filename, Ne, dt
end

function just_generate_name(DOF, N, Nover, flux, periodic)
    Ne = round(Int, DOF / (N+1))
    Nˣ = Ne
    Nʸ = Ne
    filename = "flux_" * strip(string(flux), ['(', ')']) * "_p" * string(N) * "_N" * string(Ne) * "_Nover" * string(Nover) * "_periodicity_" * string(periodic)
    return filename
end

function nameprettifier(name)
    tmp = split(name, "_")
    dof = "polynomial order " * tmp[3][2] * " and " * tmp[4][2:end] * " elements"
    stringflux = ", a " * tmp[2] 
    intorder = ", and integration order 2(" * tmp[3][2] * "+" * tmp[5][end] * ")- 1," 
    return dof * stringflux * intorder
end

