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

function generate_name_2(DOF, N, Nover, flux, periodic; L = 4π, mpicomm = MPI.COMM_WORLD, endtime = 200)
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
    cflgrid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = Array,
        polynomialorder = N + 0,
    )
    Δx =  min_node_distance(cflgrid)
    cfl = 0.15
    dt = cfl * Δx / √10

    grid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = Array,
        polynomialorder = N + Nover,
    )
    # run
    timeend = FT(endtime) # s
    nout = round(Int, 2 / dt)
    dt = 2 / nout
    filename = "new_flux_" * strip(string(flux), ['(', ')']) * "_p" * string(N) * "_N" * string(Ne) * "_Nover" * string(Nover) * "_periodicity_" * string(periodic)
    f = jldopen(filename * ".jld2", "a+")
    f["grid"] = grid
    close(f)
    return filename, Ne, dt
end

function just_generate_name_2(DOF, N, Nover, flux, periodic; L = 4π, mpicomm = MPI.COMM_WORLD, endtime = 200)
    Ne = round(Int, DOF / (N+1))
    Nˣ = Ne
    Nʸ = Ne
    Lˣ = L 
    Lʸ = L
    filename = "new_flux_" * strip(string(flux), ['(', ')']) * "_p" * string(N) * "_N" * string(Ne) * "_Nover" * string(Nover) * "_periodicity_" * string(periodic)
    return filename
end

function just_generate_name(DOF, N, Nover, flux, periodic; L = 4π, mpicomm = MPI.COMM_WORLD, endtime = 200)
    Ne = round(Int, DOF / (N+1))
    Nˣ = Ne
    Nʸ = Ne
    Lˣ = L 
    Lʸ = L
    filename = "flux_" * strip(string(flux), ['(', ')']) * "_p" * string(N) * "_N" * string(Ne) * "_Nover" * string(Nover) * "_periodicity_" * string(periodic)
    return filename
end

function nameprettifier(name)
    tmp = split(name, "_")
    dof = "polynomial order " * tmp[3+1][2] * " and " * tmp[4+1][2:end] * " elements"
    stringflux = ", a " * tmp[1+2] 
    intorder = ", and integration order 2(" * tmp[1+3][2] * "+" * tmp[1+5][end] * ")- 1," 
    return dof * stringflux * intorder
end
