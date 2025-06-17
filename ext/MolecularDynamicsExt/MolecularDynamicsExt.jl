module MolecularDynamicsExt

using TDEP
using Molly
using Unitful
using Printf
using ASEconvert

include("MDTDEP.jl")

# Needed for sTDEP with MD as force calculator
function TDEP.single_point_forces(sys::System{3})
    return Molly.forces(sys)
end


end