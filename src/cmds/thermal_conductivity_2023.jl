export ThermalConductivity2023

"""
    ThermalConductivity2023(; readiso = false, qpoint_grid, integrationtype = 2,
                           sigma = 1.0, threshold = 4.0, readqmesh = false,
                           temperature = -1.0, temperature_range,
                           logtempaxis = false, max_mfp = -1.0,
                           dumpgrid = false, noisotope = false)

# Arguments
- `readiso::Bool = false`: Read the isotope distribution from `infile.isotopes` (format specified in docs).
- `qpoint_grid::NTuple{3,Int} = (26, 26, 26)`: Density of q-point mesh for Brillouin zone integrations.
- `integrationtype::Int = 2`: Integration type for phonon DOS integration:  
  1. Gaussian  
  2. Adaptive Gaussian  
  3. Tetrahedron
- `sigma::T = 1.0`: Global scaling factor for adaptive Gaussian smearing.
- `threshold::T = 4.0`: Number of standard deviations beyond which the Gaussian is considered zero.
- `readqmesh::Bool = false`: Read the q-point mesh from file (see `genkpoints` utility).
- `temperature::T = -1.0`: Evaluate thermal conductivity at a single temperature.
- `temperature_range::NTuple{3,T} = (100, 300, 5)`: Series of temperatures (min, max, number of points).
- `logtempaxis::Bool = false`: Space temperature points logarithmically instead of linearly.
- `max_mfp::T = -1.0`: Limit the phonon mean free path (approximate domain size).
- `dumpgrid::Bool = false`: Write files with q-vectors, frequencies, eigenvectors, and group velocities for a grid.
- `noisotope::Bool = false`: Do not consider isotope scattering.

# Description
Calculates lattice thermal conductivity via the iterative solution of the phonon Boltzmann equation. Cumulative plots and raw data dumps of intermediate values are available.

**Note:** Legacy implementation; use the improved `thermal_conductivity` program when possible.

TDEP Documentation: https://tdep-developers.github.io/tdep/program/thermal_conductivity_2023/
"""
Base.@kwdef struct ThermalConductivity2023{T} <: TDEP_Command{T}
    readiso::Bool = false
    qpoint_grid::NTuple{3,Int} = (26, 26, 26)
    integrationtype::Int = 2
    sigma::T = 1.0
    threshold::T = 4.0
    readqmesh::Bool = false
    temperature::T = -1.0
    temperature_range::NTuple{3,T} = (100, 300, 5)
    logtempaxis::Bool = false
    max_mfp::T = -1.0
    dumpgrid::Bool = false
    noisotope::Bool = false
end

cmd_name(::ThermalConductivity2023) = "thermal_conductivity_2023"