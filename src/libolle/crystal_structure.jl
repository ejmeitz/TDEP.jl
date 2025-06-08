# Types from the libolle/type_crystalstructure


# https://github.com/tdep-developers/tdep/blob/898f7a1ef2c764b471e344f45ea3217e49c3fef7/src/libolle/type_crystalstructure.f90#L167
Base.@kwdef struct CrystalStructure
    isotope::Vector{IsotopeDistribution{DefaultFloat}}
    alloyspecies::Vector{AlloyAtom{DefaultFloat}}
    bz::BrillouinZone
    ibz_wedge::IBZWedge
    space_group::SymSet
    classification::CrystalStructureClassification
    mag_info::MagInfo
    natoms::Int
    lavec::Vector{Vector{DefaultFloat}}
    inv_lavec::Vector{Vector{DefaultFloat}}
    recip_lavec::Vector{Vector{DefaultFloat}}
    inv_recip_lavec::Vector{Vector{DefaultFloat}}
    volume::DefaultFloat
    n_elements::Cint
    element_counts::Vector{Cint}
    atomic_symbols::Vector{Symbol}
    atomic_numbers::Vector{Cint}
    species::Vector{Cint}
    positions::Vector{Vector{DefaultFloat}}
    positions_frac::Vector{Vector{DefaultFloat}}
    velocities::Vector{Vector{DefaultFloat}} = Vector{Vector{DefaultFloat}}[]
    forces::Vector{Vector{DefaultFloat}} = Vector{Vector{DefaultFloat}}[]
    displacements::Vector{Vector{DefaultFloat}} = Vector{Vector{DefaultFloat}}[]
    masses::Vector{DefaultFloat}
    inv_sqrt_mass::Vector{DefaultFloat}
    inelastic_neutron_cross_section::Vector{DefaultFloat} = DefaultFloat[]
    flavors::Vector{Cint} = Cint[]
end


function CrystalStructure(filename::String; verbosity::Cint = typemin(Cint))
    return @ccall "TDEP".readfromfile(
        filename::Cstring, verbosity::Ref{Cint}
    )::CrystalStructure
end


# struct BZWedge
#     n_points::Cint = typemin(Cint)
#     nodes::Vector{Cint}
#     plane::Plane
# end

# Base.@kwdef struct IBZWedge
#     n_nodes::Cint = typemin(Cint)
#     n_faces::Cint = typemin(Cint)
#     bz_wedges::Vector{BZWedge}
#     high_symmetry_points::Matrix{DefaultFloat}
#     labels::Vector{string}
# end


# #https://github.com/tdep-developers/tdep/blob/898f7a1ef2c764b471e344f45ea3217e49c3fef7/src/libolle/type_crystalstructure.f90#L110
# Base.@kwdef struct CrystalStructureClassification
#     title::String = "empty header"
#     alloy::Bool = false
#     collinear_mag_mom::Bool = false
#     noncollinear_mag_mom::Bool = false
#     have_space_group::Bool = false
#     have_bz::Bool = false
#     have_bravais::Bool = false
#     points_labeled::Bool = false
#     decided_time_reversal::Bool = false
#     have_character_table::Bool = false
#     unique_primitive_basis::Matrix{DefaultFloat} = zeros(DefaultFloat, 3, 3)
#     unique_conventional_basis::Matrix{DefaultFloat} = zeros(DefaultFloat, 3, 3)
#     permutation_to_unique::Matrix{DefaultFloat} = zeros(DefaultFloat, 3, 3)
#     transformation_to_unique::Matrix{DefaultFloat} = zeros(DefaultFloat, 3, 3)
#     unit_cell_lattice_parameter::DefaultFloat = -DefaultFloat(Inf)
#     bravais_lattice::String = "nothing"
#     bravais_lattice_long_name::String = "nothing"
#     supercell::Bool = false
#     supercell_matrix::Matrix{DefaultFloat} = zeros(DefaultFloat, 3, 3)
#     cell_index::Matrix{Cint}
#     index_in_unitcell::Matrix{Cint}
#     verbosity::Cint = typemin(Cint)
# end