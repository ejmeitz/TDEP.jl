export GenerateStructure

"""
    GenerateStructure(; dimensions, nondiagonal_dimensions,
                      output_format = 1, desired_na = -1)

# Arguments
- `dimensions::NTuple{3,Int} = (5, 5, 5)`: Dimensions of the supercell.
- `nondiagonal_dimensions::NTuple{9,Int} = (0, 0, 0, 0, 0, 0, 0, 0, 0)`: Non-diagonal dimensions of the supercell.
- `output_format::Int = 1`: Output format.  
  - `1` = VASP  
  - `2` = Abinit  
  - `4` = FHI-Aims  
  - `5` = Siesta
- `desired_na::Int = -1`: Desired number of atoms in the supercell; tries to choose a cell as cubic as possible.

# Description
Builds supercells (both diagonal and non-diagonal) and can find the optimal supercell for a given lattice—handy for complicated structures.

TDEP Documentation: https://tdep-developers.github.io/tdep/program/generate_structure/
"""
Base.@kwdef struct GenerateStructure <: TDEP_Command{Nothing}
    dimensions::NTuple{3,Int} = (5, 5, 5)
    nondiagonal_dimensions::NTuple{9,Int} = (0, 0, 0, 0, 0, 0, 0, 0, 0)
    output_format::Int = 1
    desired_na::Int = -1
end

cmd_name(::GenerateStructure) = "generate_structure"