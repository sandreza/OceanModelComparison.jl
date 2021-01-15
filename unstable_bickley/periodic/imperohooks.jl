# using Impero
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
include(pwd() * "/imperohooks/utils.jl")

