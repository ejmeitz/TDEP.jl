module TDEP

using TDEP_jll
using MPI

using DelimitedFiles
using ProgressMeter
using Printf

# Define interfaces for systems and force calculation
using AtomsBase
using AtomsCalculators

# Useful for end-user to have access to automatically
using Reexport 
@reexport using Unitful
@reexport using ASEconvert
@reexport using PythonCall

const DefaultFloat = Float64

for file in readdir(joinpath(@__DIR__, "cmds"))
    include(joinpath("./cmds", file))
end

include("io.jl")
include("stubs.jl")


end
