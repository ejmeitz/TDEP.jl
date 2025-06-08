module TDEP

using TDEP_jll

const DefaultFloat = Float64

for file in readdir(@__DIR__, "cmds")
    include(joinpath("./cmds", file))
end

end
