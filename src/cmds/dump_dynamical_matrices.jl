export DumpDynamicalMatrices

"""
    DumpDynamicalMatrices(; qpoint_grid = (26, 26, 26), meshtype = 1, readqpoints = false)

# Arguments
- `qpoint_grid::NTuple{3,Int} = (26, 26, 26)`: Density of the q-point mesh for Brillouin zone integrations.
- `meshtype::Int = 1`: Q-point mesh type:  
  1. Monkhorst–Pack  
  2. FFT  
  3. Wedge-based mesh (approximate density)  
  4. Commensurate mesh of an approximately cubic supercell
- `readqpoints::Bool = false`: Read q-points from `infile.dynmatqpoints` (fractional coordinates) instead of generating a mesh.

# Description
Utility to extract TDEP dynamical matrices in formats readable by other codes (currently supports ABINIT DDB format).

TDEP Documentation: https://tdep-developers.github.io/tdep/program/dump_dynamical_matrices/
"""
Base.@kwdef struct DumpDynamicalMatrices <: TDEP_Command{Nothing}
    qpoint_grid::NTuple{3,Int} = (26, 26, 26)
    meshtype::Int = 1
    readqpoints::Bool = false
end

cmd_name(::DumpDynamicalMatrices) = "dump_dynamical_matrices"