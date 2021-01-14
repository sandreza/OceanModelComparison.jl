using Test
using ClimateMachine
using ClimateMachine.Mesh.Grids
using ClimateMachine.Mesh.Topologies
using ClimateMachine.MPIStateArrays
using Impero, Printf, MPI, LinearAlgebra, Statistics, GaussQuadrature

ClimateMachine.init()
const ArrayType = ClimateMachine.array_type()
const mpicomm = MPI.COMM_WORLD
const FT = Float64

N = 8
Ne_vert  = 1
Ne_horz  = 3
Rrange = range(0.5; length = Ne_vert + 1, stop = FT(1))

topl = StackedCubedSphereTopology(mpicomm, Ne_horz, Rrange)
spheregrid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = ArrayType,
        polynomialorder = N,
        meshwarp = Topologies.cubedshellwarp,
)

scene = visualize(spheregrid)