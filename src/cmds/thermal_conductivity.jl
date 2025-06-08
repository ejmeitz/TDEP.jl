export ThermalConductivity

"""
    ThermalConductivity(; readiso = false, qpoint_grid, qpoint_grid3ph,
                       qpoint_grid4ph, integrationtype = 2, sigma = 1.0,
                       readqmesh = false, fourthorder = false,
                       classical = false, temperature = 300.0,
                       max_mfp = -1.0, noisotope = false,
                       iterative_tolerance = 1e-5, iterative_maxsteps = 200,
                       seed = -1, mfppts = 1000, freqpts = 1000,
                       dosintegrationtype = 3, dossigma = 1.0)

# Arguments
- `readiso::Bool = false`: Read the isotope distribution from `infile.isotopes`.
- `qpoint_grid::NTuple{3,Int} = (26, 26, 26)`: Density of q-point mesh for Brillouin zone integrations.
- `qpoint_grid3ph::NTuple{3,Int} = (-1, -1, -1)`: Dimensions of the grid for three-phonon integration.
- `qpoint_grid4ph::NTuple{3,Int} = (-1, -1, -1)`: Dimensions of the grid for four-phonon integration.
- `integrationtype::Int = 2`: Integration type for phonon DOS:  
  1. Gaussian  
  2. Adaptive Gaussian
- `sigma::T = 1.0`: Global scaling factor for Gaussian/adaptive Gaussian smearing.
- `readqmesh::Bool = false`: Read the q-point mesh from file (see `genkpoints` utility).
- `fourthorder::Bool = false`: Include four-phonon contributions to scattering.
- `classical::Bool = false`: Use the classical limit for phonon occupation and heat capacity.
- `temperature::T = 300.0`: Evaluate thermal conductivity at a single temperature.
- `max_mfp::T = -1.0`: Limit phonon mean free path (approximate domain size in m).
- `noisotope::Bool = false`: Do not consider isotope scattering.
- `iterative_tolerance::T = 1e-5`: Tolerance for the iterative solution.
- `iterative_maxsteps::Int = 200`: Maximum iterations for the solver.
- `seed::Int = -1`: Seed for Monte-Carlo grids (positive integer).
- `mfppts::Int = 1000`: Number of mean-free-path points for cumulative conductivity.
- `freqpts::Int = 1000`: Number of frequency points for spectral conductivity.
- `dosintegrationtype::Int = 3`: Integration type for spectral conductivity:  
  1. Gaussian  
  2. Adaptive Gaussian  
  3. Tetrahedron
- `dossigma::T = 1.0`: Scaling factor for Gaussian integration in spectral conductivity.

# Description
Calculates lattice thermal conductivity using the mode-coupling formalism, including collective and off-diagonal contributions up to fourth order.

TDEP Documentation: https://tdep-developers.github.io/tdep/program/thermal_conductivity/
"""
Base.@kwdef struct ThermalConductivity{T} <: TDEP_Command{T}
    readiso::Bool = false
    qpoint_grid::NTuple{3,Int} = (26, 26, 26)
    qpoint_grid3ph::NTuple{3,Int} = (-1, -1, -1)
    qpoint_grid4ph::NTuple{3,Int} = (-1, -1, -1)
    integrationtype::Int = 2
    sigma::T = 1.0
    readqmesh::Bool = false
    fourthorder::Bool = false
    classical::Bool = false
    temperature::T = 300.0
    max_mfp::T = -1.0
    noisotope::Bool = false
    iterative_tolerance::T = 1e-5
    iterative_maxsteps::Int = 200
    seed::Int = -1
    mfppts::Int = 1000
    freqpts::Int = 1000
    dosintegrationtype::Int = 3
    dossigma::T = 1.0
end

cmd_name(::ThermalConductivity) = "thermal_conductivity"