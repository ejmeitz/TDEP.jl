export Lineshape

"""
    Lineshape(; unit = "thz", temperature = 300.0, n_energies = 1200,
              qpoint_grid = (26, 26, 26), meshtype = 2,
              integrationtype = 2, sigma = 1.0, path = false,
              readpath = false, nq_on_path = 100,
              qpoint = (0.0, 0.0, 0.0), highsymmetrypoint = "",
              max_energy = 1.4, no_isotope_scattering = false,
              no_thirdorder_scattering = false, qdirin = (1.0, 0.0, 0.0),
              grid = false, readiso = false, readqmesh = false)

# Arguments
- `unit::String = "thz"`  
  Output unit:  
  - `"thz"` = terahertz (frequency)  
  - `"icm"` = inverse cm  
  - `"mev"` = meV  
- `temperature::T = 300.0`  
  Temperature for occupation numbers (should match force-constant determination temperature).  
- `n_energies::Int = 1200`  
  Number of energy points for the self-energy calculation.  
- `qpoint_grid::NTuple{3,Int} = (26, 26, 26)`  
  Density of q-point mesh for Brillouin zone integrations.  
- `meshtype::Int = 2`  
  Q-point mesh type:  
  1. Monkhorst–Pack  
  2. FFT  
  3. Wedge-based mesh (approximate density)  
- `integrationtype::Int = 2`  
  Integration type for phase-space integrals:  
  1. Gaussian  
  2. Adaptive Gaussian  
  3. Tetrahedron  
  4. …  
  5. …  
- `sigma::T = 1.0`  
  Global scaling factor for Gaussian/adaptive Gaussian smearing.  
- `path::Bool = false`  
  Calculate self-energy and spectral function along a high-symmetry path.  
- `readpath::Bool = false`  
  Read the path from `infile.qpoints_dispersion`.  
- `nq_on_path::Int = 100`  
  Number of q-points between each high-symmetry point when `path = true`.  
- `qpoint::NTuple{3,T} = (0.0, 0.0, 0.0)`  
  Single q-point (fractional coordinates) for self-energy calculation.  
- `highsymmetrypoint::String = ""`  
  Label of a high-symmetry point (e.g. `"X"`, `"L"`) instead of explicit `qpoint`.  
- `max_energy::T = 1.4`  
  Maximum energy cutoff (in multiples of max harmonic frequency).  
- `no_isotope_scattering::Bool = false`  
  Disable isotope (mass-disorder) scattering.  
- `no_thirdorder_scattering::Bool = false`  
  Disable three-phonon scattering.  
- `qdirin::NTuple{3,T} = (1.0, 0.0, 0.0)`  
  Incident wavevector (Cartesian) for non-analytical zone-center behavior.  
- `grid::Bool = false`  
  Calculate spectral functions on a q-point grid.  
- `readiso::Bool = false`  
  Read isotope distribution from `infile.isotopes`.  
- `readqmesh::Bool = false`  
  Read q-point mesh from file (see `genkpoints` utility).

# Description
Calculate frequency-dependent phonon self-energy and spectral function, either along a path or on a grid, with options to disable specific scattering channels and control mesh/integration parameters.

TDEP Documentation: https://tdep-developers.github.io/tdep/program/lineshape/
"""
Base.@kwdef struct Lineshape{T} <: TDEP_Command{T}
    unit::String = "thz"
    temperature::T = 300.0
    n_energies::Int = 1200
    qpoint_grid::NTuple{3,Int} = (26, 26, 26)
    meshtype::Int = 2
    integrationtype::Int = 2
    sigma::T = 1.0
    path::Bool = false
    readpath::Bool = false
    nq_on_path::Int = 100
    qpoint::NTuple{3,T} = (0.0, 0.0, 0.0)
    highsymmetrypoint::String = ""
    max_energy::T = 1.4
    no_isotope_scattering::Bool = false
    no_thirdorder_scattering::Bool = false
    qdirin::NTuple{3,T} = (1.0, 0.0, 0.0)
    grid::Bool = false
    readiso::Bool = false
    readqmesh::Bool = false
end

cmd_name(::Lineshape) = "lineshape"