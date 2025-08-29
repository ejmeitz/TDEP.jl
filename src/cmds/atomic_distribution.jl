export AtomicDistribution

"""
    AtomicDistribution(; cutoff = 8.0, stride = 1)

# Arguments
- `cutoff::Float64 = 8.0`: Consider atomic pairs up to this distance (Å).
- `stride::Int = 1`: Use every Nth configuration instead of all.

# Description
Calculates properties of the atomic distribution from molecular dynamics—mean-square displacement, pair distribution function, vector distribution functions, and probability densities. Useful for analyzing simulations near instabilities or phase transitions to understand atomic positions.

TDEP Documentation: https://tdep-developers.github.io/tdep/program/atomic_distribution/
"""
Base.@kwdef struct AtomicDistribution <: TDEP_Command
    cutoff::Float64 = 8.0
    stride::Int = 1
end

cmd_name(::AtomicDistribution) = "atomic_distribution"