module LAMMPSExt

using LAMMPS
using TDEP
using Unitful
import AtomsBase
import LinearAlgebra

export single_point_potential_energy

# Modified from the Molly LAMMPS Calculator which I implemented
# this lets me remove some things and also not have to depend on Molly

function TDEP.LAMMPSCalculator(
        sys::TDEPSystem,
        potential_definition::Union{String, Array{String}};
        label_type_map::Dict{Symbol, Int} = Dict{Symbol, Int}(),
        logfile_path::String = "none",
    )

    cell_matrix = reduce(hcat, AtomsBase.cell_vectors(sys))

    if !LinearAlgebra.isdiag(cell_matrix)
        error(ArgumentError("LAMMPSCalculator does not support triclinic systems yet. Got non-diagonal cell."))
    end

    if length(potential_definition) == 0
        error(ArgumentError("Cannot pass emptry string as potential definition for LAMMPSCalculator"))
    end

    lmp = LMP(["-screen","none"], LAMMPS.MPI.COMM_WORLD)

    all_syms = AtomsBase.atomic_symbol(sys, :)
    unique_syms = unique(all_syms)
    unique_sym_idxs = Dict(sym => findfirst(x -> x == sym, all_syms) for sym in unique_syms)

    if any(unique_syms .== :unknown)
        error(ArgumentError("All atoms must have atomic symbols to use LAMMPSCalculator"))
    end
    
    ids = collect(Int32, 1:length(sys))
    xhi, yhi, zhi = ustrip.(u"angstrom", LinearAlgebra.diag(cell_matrix))

    if length(label_type_map) == 0
        label_type_map = Dict(sym => Int32(i) for (i, sym) in enumerate(unique_syms))
        types = [label_type_map[sym] for sym in all_syms]
    else 
        unique_sym_user = keys(label_type_map)
        if Set(unique_sym_user) != Set(unique_syms)
            error(ArgumentError("You provided a label_type_map with $(unique_sym_user) symbols, but" *
                " the system has $(unique_syms). They must match exactly if you pass label_type_map."))
        end
        types = [Int32(label_type_map[sym]) for sym in all_syms]
    end

    m_lmp = Dict(label_type_map[sym] => sys.masses[i] for (sym, i) in unique_sym_idxs)

    label_map_cmd = "labelmap atom " * join(["$(i) $(sym)" for (sym,i) in label_type_map], " ") 

    setup_cmd = """
            log $(logfile_path)
            units $(lammps_unit_system)
            atom_style atomic
            atom_modify map array sort 0 0
        """
    
    cell_cmd = """
            boundary p p p
            region cell block 0 $(xhi) 0 $(yhi) 0 $(zhi) units box
            create_box $(length(unique_syms)) cell
            $(label_map_cmd)
        """
    
    mass_cmd = join(["mass $(type) $(m)" for (type,m) in  m_lmp], "\n")

    command(lmp, setup_cmd)
    command(lmp, cell_cmd)
    command(lmp, mass_cmd)


    LAMMPS.create_atoms(
        lmp,
        reinterpret(reshape, Float64, sys.position),
        ids,
        types
    )   

    try
        command(lmp, potential_definition)
    catch e
        if startswith(e.msg, "Number of element to type mappings does")
            @info "Ensure path to potential definition is wrapped in quotes if there are spaces in path."
        end
        rethrow(e)
    end

    command(lmp, "compute pot_e all pe")

    # This allows LAMMPS to register the computes/fixes
    # and build the neighbor list. 
    command(lmp, "run 0 post no")

    return LAMMPSCalculator{typeof(lmp)}(lmp, -1)
end

# Expect Vector of Vectors or 3 x N Matrix
function single_point_potential_energy(r::AbstractVecOrMat, inter::LAMMPSCalculator)
    scatter!(lammps_calc.lmp, "x", reinterpret(reshape, Float64, r))
    command(lammps_calc.lmp, "run 0 pre no post no")
    return extract_compute(inter.lmp, "pot_e", STYLE_GLOBAL, TYPE_SCALAR)[1] * u"eV"
end

end