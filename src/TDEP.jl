module TDEP

using TDEP_jll
using MPI
using StaticArrays
using DelimitedFiles
using ProgressMeter
using Printf
using FastGaussQuadrature

# Define interfaces for systems and force calculation
using AtomsBase
using AtomsCalculators

# Useful for end-user to have access to automatically
using Reexport 
@reexport using Unitful


const DefaultFloat = Float64

for file in readdir(joinpath(@__DIR__, "cmds"))
    include(joinpath("./cmds", file))
end

include("io.jl")
include("stubs.jl")
include("workflows/sTDEP.jl")

# Include the simplified TDEP wrappers
include("libolle/crystal_structure.jl")

end
