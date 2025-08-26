export PhononDispersionRelations

"""
    PhononDispersionRelations(; unit = "thz", nq_on_path = 100, readpath = false,
                              dos = false, qpoint_grid, meshtype = 1,
                              sigma = 1.0, readqmesh = false,
                              integrationtype = 2, dospoints = 400,
                              gruneisen = false, dumpgrid = false,
                              temperature = -1.0)

# Arguments
- `unit::String = "thz"`  
  Output unit:  
  - `"thz"` = terahertz (frequency)  
  - `"icm"` = inverse cm  
  - `"mev"` = meV  
- `nq_on_path::Int = 100`  
  Number of q-points between each high-symmetry point on the dispersion path.  
- `readpath::Bool = false`  
  Read the q-point path from `infile.qpoints_dispersion` (generates an example from the crystal structure).  
- `dos::Bool = false`  
  Calculate the phonon density of states (DOS).  
- `qpoint_grid::NTuple{3,Int} = (26, 26, 26)`  
  Density of q-point mesh for Brillouin zone integrations.  
- `meshtype::Int = 1`  
  Q-point mesh type:  
  1. Monkhorst–Pack  
  2. FFT  
  3. Wedge-based mesh (approximate density)  
  4. Commensurate mesh of an approximately cubic supercell  
- `sigma::Float64 = 1.0`  
  Global scaling factor for Gaussian/adaptive Gaussian smearing (baseline set procedurally).  
- `readqmesh::Bool = false`  
  Read the q-point mesh from file (see `genkpoints` utility).  
- `integrationtype::Int = 2`  
  Integration type for the phonon DOS:  
  1. Gaussian  
  2. Adaptive Gaussian  
  3. Tetrahedron  
- `dospoints::Int = 400`  
  Number of points on the frequency axis of the DOS.  
- `gruneisen::Bool = false`  
  Use third-order force constants to calculate mode Grüneisen parameters.  
- `dumpgrid::Bool = false`  
  Write files with q-vectors, frequencies, eigenvectors, and group velocities for a grid.
- `temperature::Float64 = -1`
  If positive, write an `outfile.free_energy` file which has the vibrational free energy, 
  entropy and heat capacity (always quantum) 

# Description
Calculate phonon dispersion relations and related quantities. By default, dispersion curves along a standard high-symmetry path are computed. Options allow calculation of mode Grüneisen parameters, projected phonon DOS, thermodynamic properties, or full data dumps.

TDEP Documentation: https://tdep-developers.github.io/tdep/program/phonon_dispersion_relations/
"""
Base.@kwdef struct PhononDispersionRelations <: TDEP_Command{Float64}
    unit::String = "thz"
    nq_on_path::Int = 100
    readpath::Bool = false
    dos::Bool = false
    qpoint_grid::NTuple{3,Int} = (26, 26, 26)
    meshtype::Int = 1
    sigma::Float64 = 1.0
    readqmesh::Bool = false
    integrationtype::Int = 2
    dospoints::Int = 400
    gruneisen::Bool = false
    dumpgrid::Bool = false
    temperature::Float64 = -1.0
end

cmd_name(::PhononDispersionRelations) = "phonon_dispersion_relations"