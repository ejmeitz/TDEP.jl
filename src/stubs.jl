export generate_MDTDEP_dataset, NVT, TI, LAMMPSCalculator


# This file contains function handles which are implemented
# as part of package extensions. These are separated as they
# require a heavy dependency (Molly.jl) which not all
# users may want to download.

function generate_MDTDEP_dataset end

struct NVT
    thermostat
    n_steps_warmup::Integer
    n_steps::Integer
    sample_every::Integer
end

function TI end


"""
    LAMMPSCalculator(
        sys::System{3, AT, T},
        potential_definition::Union{String, Array{String}};
        logfile_path::String = "none",
    )

Defines a general interaction that will call LAMMPS to calculate forces and energies. Forces
and energies are calculated on a single thread. You must call LAMMPS.MPI.Init() for LAMMPS.jl
to load the LAMMPS executable on systems where MPI is available.

The LAMMPS potential files can be found at:
`abspath(dirname(LAMMPS.locate()), "..", "share", "lammps", "potentials")`

Restrictions:
-------------
- CPU only
- Floats promote to Float64
- No triclinic boundary
- 3D systems only
- Fully periodic systems only

Arguments:
----------
- `sys::System{3}`: The system object this interaction will be applied to. You still have to add this
    interaction to the system object yourself after constructing the LAMMPSCalculator.
- `potential_definition::Union{String, Array{String}}` : Commands passed to lammps which define your interaction.
    For example, to define LJ you pass:
    `lj_cmds = ["pair_style lj/cut 8.5", "pair_coeff * * 0.0104 3.4", "pair_modify shift yes"]`
- `label_type_map::Dict{Symbol, Int} = Dict{Symbol, Int}()` : By default atom types are assigned in the 
    order they appear in the system. This can make defining the potential for multi-atomic systems
    difficult. By providing this dictionary you can overide the type label assigned to each unique species. 
- `logfile_path::String = "none"` : Path where LAMMPS logfile is written. Defaults to no log file. 
"""
mutable struct LAMMPSCalculator{T, E, L}
    lmp::T # T will be LMP but that is not available here
    last_updated::Int
    energy_unit::E
    length_unit::L
end