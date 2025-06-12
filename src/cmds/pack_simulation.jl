export PackSimulation

"""
    PackSimulation(; stride = 1, nrand = -1, temperature = -1.0,
                    variable_cell = false, magnetic_moments = false,
                    dielectric = false, molecular_dynamics = false,
                    notidy = false)

# Arguments
- `stride::Int = 1`: Pack every Nth configuration instead of all.
- `nrand::Int = -1`: Pack N random configurations instead of all.
- `temperature::T = -1.0`: Override the simulation temperature.
- `variable_cell::Bool = false`: Store variable cell information.
- `magnetic_moments::Bool = false`: Store projected magnetic moments.
- `dielectric::Bool = false`: Store dielectric constant and Born charges.
- `molecular_dynamics::Bool = false`: Mark dataset as real molecular dynamics.
- `notidy::Bool = false`: Skip default cleanup (drift removal, etc.).

# Description
Utility to pack a simulation into HDF5 format with optional subsampling, random selection, and storage of additional metadata (cell variability, magnetic moments, dielectric properties, MD flag).

TDEP Documentation: https://tdep-developers.github.io/tdep/program/pack_simulation/
"""
Base.@kwdef struct PackSimulation{T} <: TDEP_Command{T}
    stride::Int = 1
    nrand::Int = -1
    temperature::T = -1.0
    variable_cell::Bool = false
    magnetic_moments::Bool = false
    dielectric::Bool = false
    molecular_dynamics::Bool = false
    notidy::Bool = false
end

cmd_name(::PackSimulation) = "pack_simulation"