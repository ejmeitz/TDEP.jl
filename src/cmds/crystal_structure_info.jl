export CrystalStructureInfo

"""
    CrystalStructureInfo(; printsymmetry = false)

# Arguments
- `printsymmetry::Bool = false`: Also print the symmetry operations.

# Description
Diagnostic tool to verify symmetry heuristics: prints the identified Bravais lattice, inequivalent high-symmetry points in the Brillouin zone, and outputs the Brillouin zone and its irreducible wedges as polyhedra for inspection. Optionally prints all symmetry operations.

TDEP Documentation: https://tdep-developers.github.io/tdep/program/crystal_structure_info/
"""
Base.@kwdef struct CrystalStructureInfo <: TDEP_Command
    printsymmetry::Bool = false
end


cmd_name(::CrystalStructureInfo) = "crystal_structure_info"