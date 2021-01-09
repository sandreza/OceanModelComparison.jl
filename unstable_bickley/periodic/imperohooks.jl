using Impero
using ClimateMachine

using ClimateMachine.Mesh.Grids
using ClimateMachine.Mesh.Elements
import ClimateMachine.Mesh.Elements: baryweights
using GaussQuadrature
using Base.Threads
using Makie, GLMakie, AbstractPlotting
using ImageTransformations, Colors
using AbstractPlotting.MakieLayout, LinearAlgebra


include(pwd() * "/imperohooks/permutations.jl")
include(pwd() * "/imperohooks/find_element.jl")
include(pwd() * "/imperohooks/lagrange_interpolation.jl")
include(pwd() * "/imperohooks/fields.jl")
include(pwd() * "/imperohooks/gridhelper.jl")



"""
function cellaverage(Q; M = nothing)
# Description
Compute the cell-average of Q given the mass matrix M.
Assumes that Q and M are the same size
# Arguments
`Q`: MPIStateArrays (array)
# Keyword Arguments
`M`: Mass Matrix (array)
# Return
The cell-average of Q
"""
function cellaverage(Q; M = nothing)
    if M==nothing
        return nothing
    end
    return (sum(M .* Q, dims = 1) ./ sum(M , dims = 1))[:]
end

"""
function coordinates(grid::DiscontinuousSpectralElementGrid)
# Description
Gets the (x,y,z) coordinates corresponding to the grid
# Arguments
- `grid`: DiscontinuousSpectralElementGrid
# Return
- `x, y, z`: views of x, y, z coordinates
"""
function coordinates(grid::DiscontinuousSpectralElementGrid)
    x = view(grid.vgeo, :, grid.x1id, :)   # x-direction	
    y = view(grid.vgeo, :, grid.x2id, :)   # y-direction	
    z = view(grid.vgeo, :, grid.x3id, :)   # z-direction
    return x, y, z
end

"""
function cellcenters(Q; M = nothing)
# Description
Get the cell-centers of every element in the grid
# Arguments
- `grid`: DiscontinuousSpectralElementGrid
# Return
- Tuple of cell-centers
"""
function cellcenters(grid::DiscontinuousSpectralElementGrid)
    x, y, z = coordinates(grid)
    M = view(grid.vgeo, :, grid.Mid, :)  # mass matrix
    xC = cellaverage(x, M = M)
    yC = cellaverage(y, M = M)
    zC = cellaverage(z, M = M)
    return xC[:], yC[:], zC[:]
end
