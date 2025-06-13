module TDEP

using TDEP_jll
using MPI
using DelimitedFiles
using Printf

const DefaultFloat = Float64

for file in readdir(joinpath(@__DIR__, "cmds"))
    include(joinpath("./cmds", file))
end

end
