export CanonicalConfiguration

"""
    CanonicalConfiguration(; temperature, nconf, quantum = false, outputformat = 1,
                            mindist = -1, debye_temperature = -1.0, maximum_frequency = -1.0)

# Arguments 
- `temperature::T`: Temperature to emulate.
- `nconf::Int`: Number of configurations to generate.
- `quantum::Bool = false` : Use Bose-Einstein statistics instead of Maxwell-Boltzmann.
- `output_format::Int = 1`: (1) VASP (2) Abinit (4) FHI-Aims (5) Siesta
- `mindist::T = -1`: Smallest distance between two atoms allowed, in units of the nearest neighbour distance.
- `debye_temperature::T = -1.0`: Generate forceconstants that match a Debye temperature,
                                 and build displacements according to these. 
- `maximum_frequency::T = -1.0`: Generate forceconstants that match a maximum frequency (in THz),
                                 and build displacements according to these.

# Description
This code takes a second order forceconstant and generates a supercell
with the atoms displaced corresponding to a harmonic canonical ensemble.
These configurations can be used directly, to generate force-displacement
statistics to determine force constants and other thermal properties, 
or as the starting point for molecular dynamics.

TDEP Documentatoin: https://tdep-developers.github.io/tdep/program/canonical_configuration/
"""
Base.@kwdef struct CanonicalConfiguration{T} <: TDEP_Command{T}
    temperature::T
    nconf::Int
    quantum::Bool = false
    output_format::Int = 1
    mindist::T = -1.0
    debye_temperature::T = -1.0
    maximum_frequency::T = -1.0
end

cmd_name(::CanonicalConfiguration) = "canonical_configuration"