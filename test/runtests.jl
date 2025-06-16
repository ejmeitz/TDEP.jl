using Molly
using TDEP
using Test

data_dir = joinpath(@__DIR__, "data")


#TODO Test run all commands on dummy data to ensure they are wrapped properly
#TODO Try with only manditory args as well as ALL kwargs


@testset "MD Extension" begin
    
    rundir = "/mnt/merged/emeitz/tdep_jl_test"
    a = 5.2468u"Å" # Lattice parameter for FCC Argon at 10 K
    temp = 10.0u"K"
    fcc_crystal = SimpleCrystals.FCC(a, :Ar, SVector(4, 4, 4))
    
    n_atoms = length(fcc_crystal)
    
    r_cut = 8.5u"Å"

    pot = LennardJones(cutoff=ShiftedForceCutoff(r_cut), )

    sys = System(
        fcc_crystal;
        pairwise_inters=(pot,),
        energy_units=u"eV",
        force_units=u"eV * Å^-1",
    )

    # Set sigma and epsilon parameters properly
    σ = 3.4u"Å"
    ϵ = 0.010423u"eV"
    updated_atoms = []

    for i in eachindex(sys)
        push!(updated_atoms, Molly.Atom(index=sys.atoms[i].index, atom_type=sys.atoms[i].atom_type,
                                mass=sys.atoms[i].mass, charge=sys.atoms[i].charge,
                                σ=σ, ϵ=ϵ))
    end

    sys = System(sys; atoms=[updated_atoms...])
    sim = Langevin(dt = 1.0u"fs", temperature = temp, friction = 1.0u"ps^-1")

    generate_MDTDEP_dataset(sys, sim, 25_000, 50_000, 1_000, rundir; n_seeds = 2)

    to_ucposcar(fcc_crystal, joinpath(rundir, "infile.ucposcar"))
    efc = ExtractForceConstants(secondorder_cutoff = ustrip(r_cut))
    execute(efc, rundir)

end