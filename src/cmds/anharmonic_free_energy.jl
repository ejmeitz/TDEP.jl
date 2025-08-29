export AnharmonicFreeEnergy

"""
    AnharmonicFreeEnergy(; temperature = 300.0, qpoint_grid = (26, 26, 26),
                         integrationtype = 2, sigma = 1.0,
                         temperature_range = (-1, -1, -1),
                         readiso = false, meshtype = 1,
                         readqmesh = false)

# Arguments
- `temperature::T = 300.0`: Temperature for occupation numbers (should match the force-constant determination temperature).
- `qpoint_grid::NTuple{3,Int} = (26, 26, 26)`: Density of the q-point mesh for Brillouin zone integrations.
- `integrationtype::Int = 2`: Integration type:  
  1. Gaussian  
  2. Adaptive Gaussian  
  3. Tetrahedron
- `sigma::T = 1.0`: Global scaling factor for Gaussian/adaptive Gaussian smearing.
- `temperature_range::NTuple{3,T} = (-1, -1, -1)`: Series of temperatures (min, max, number of points) for thermodynamic phonon properties.
- `readiso::Bool = false`: Read the isotope distribution from file.
- `meshtype::Int = 1`: Q-point mesh type:  
  1. Monkhorst–Pack  
  2. FFT  
  3. Wedge-based mesh  
  4. Commensurate cubic supercell mesh
- `readqmesh::Bool = false`: Read the q-point mesh from file (see `genkpoints` utility).

# Description
Compute phonon free energy including anharmonic contributions over a specified temperature range or at a single temperature, with options for mesh and integration control.

TDEP Documentation: https://tdep-developers.github.io/tdep/program/anharmonic_free_energy/
"""
Base.@kwdef struct AnharmonicFreeEnergy <: TDEP_Command
    temperature::Float64 = 300.0
    qpoint_grid::NTuple{3,Int} = (26, 26, 26)
    integrationtype::Int = 2
    sigma::Float64 = 1.0
    temperature_range::NTuple{3,Float64} = (-1.0, -1.0, -1.0)
    readiso::Bool = false
    meshtype::Int = 1
    readqmesh::Bool = false
end

cmd_name(::AnharmonicFreeEnergy) = "anharmonic_free_energy"