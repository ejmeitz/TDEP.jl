export ExtractForceConstants

"""
ExtractForceConstants(; secondorder_cutoff, thirdorder_cutoff = -1.0, fourthorder_cutoff = -1.0,
                        polar = false, stride = 1, firstorder = false, temperature = -1.0, 
                        norotational = false, nohuang = false, nohermitian = false, 
                        potential_energy_differences = true)

# Arguments
- `secondorder_cutoff::T`: Cutoff for the second order force constants.
- `thirdorder_cutoff::T = -1.0`: Cutoff for the third order force constants.
- `fourthorder_cutoff::T = -1.0`: Cutoff for the fourth order force constants.
- `polar::Bool = false`: Add dipole–dipole corrections for polar materials.
- `stride::Int = 1`: Use every N-th configuration instead of all (useful for long MD runs with linearly dependent frames).
- `firstorder::Bool = false`: Include the first order force constants (for finite-temperature equilibrium structure).
- `temperature::T = -1.0`: Temperature for the self-consistent solver.
- `norotational::Bool = false`: Turn off imposing rotational invariance (needed for 2D systems).
- `nohuang::Bool = false`: Turn off imposing Huang invariances (useful for 2D systems).
- `nohermitian::Bool = false`: Turn off imposing Hermitian invariance.

# Description
The main algorithm of the TDEP method. Starting from a symmetry analysis, this 
command finds the irreducible representation of the interatomic force constants and 
extracts them from position and force data.

TDEP Documentation: https://tdep-developers.github.io/tdep/program/extract_forceconstants/
"""

Base.@kwdef struct ExtractForceConstants{T} <: TDEP_Command{T}
    secondorder_cutoff::T
    thirdorder_cutoff::T = -1.0
    fourthorder_cutoff::T = -1.0
    polar::Bool = false
    stride::Int = 1
    firstorder::Bool = false
    temperature::T = -1.0
    norotational::Bool = false
    nohuang::Bool = false
    nohermitian::Bool = false
    potential_energy_differences::Bool = true
end

cmd_name(::ExtractForceConstants) = "extract_forceconstants"

function required_files(efc::ExtractForceConstants)

    required_files = ["infile.ucposcar", "infile.ssposcar", "infile.stat",
                      "infile.positions", "infile.forces", "infile.meta"]

    if efc.polar
        required_files = [required_files; "infile.lotosplitting"]
    end

    return required_files
end