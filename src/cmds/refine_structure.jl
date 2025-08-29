export RefineStructure

"""
    RefineStructure(; tolerance_lattice = 1e-5, tolerance_internal = 1e-5,
                    unitcell = "infile.ucposcar", supercell = "",
                    prototype = "")

# Arguments
- `tolerance_lattice::Float64 = 1e-5`: Tolerance for the lattice symmetry.
- `tolerance_internal::Float64 = 1e-5`: Tolerance for internal atomic degrees of freedom.
- `unitcell::String = "infile.ucposcar"`: Filename of the unit cell to refine.
- `supercell::String = ""`: Filename of the supercell to refine (if any).
- `prototype::String = ""`: Prototype unit cell with known correct symmetry.

# Description
Utility to ensure the input structure satisfies all symmetry constraints to high precision by refining lattice parameters and atomic positions.

TDEP Documentation: https://tdep-developers.github.io/tdep/program/refine_structure/
"""
Base.@kwdef struct RefineStructure <: TDEP_Command
    tolerance_lattice::Float64 = 1e-5
    tolerance_internal::Float64 = 1e-5
    unitcell::String = "infile.ucposcar"
    supercell::String = ""
    prototype::String = ""
end


cmd_name(::RefineStructure) = "refine_structure"