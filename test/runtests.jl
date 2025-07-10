
using SimpleCrystals
using Molly
using TDEP
using Test
using PythonCall

data_dir = joinpath(@__DIR__, "data")


#TODO Test run all commands on dummy data to ensure they are wrapped properly
#TODO Try with only manditory args as well as ALL kwargs


function make_Ar_Molly_Sys()
    a = 5.2468u"Å" # Lattice parameter for FCC Argon at 10 K
    fcc_crystal = SimpleCrystals.FCC(a, :Ar, SVector(4, 4, 4))
    
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

    return System(sys; atoms=[updated_atoms...])
end

@testset "sTDEP SimpleCrystals" begin

    using ASEconvert
    using PythonCall
    ase_lj = pyimport("ase.calculators.lj")


    basedir = "/mnt/merged/emeitz/tdep_jl_test/cc_test"

    # Settings
    a = 5.2468u"Å"
    r_cut = 8.5
    temperature = 10.0 #K
    max_freq = 2.5 #THz
    nconf = 20 # number of configurations per iteration
    niter = 2 # number of self-consistent iterations

    # Create Structure
    fcc_crystal = SimpleCrystals.FCC(a, :Ar, SVector(4, 4, 4))
    to_ucposcar(fcc_crystal, joinpath(basedir, "infile.ucposcar"))

    # Create Force Calculator
    ase_calc = ase_lj.LennardJones(sigma = 3.4, epsilon = 0.010423, rc = r_cut)
    calc = ASEcalculator(ase_calc)

    # Run sTDEP
    sTDEP(
            fcc_crystal,
            calc,
            basedir,
            niter,
            nconf,
            r_cut,
            temperature,
            max_freq
        )

end

@testset "MD Extension" begin
    
    outdir = "/mnt/merged/emeitz/tdep_jl_test"

    # Settings
    a = 5.2468u"Å" 
    r_cut = 8.5u"Å"
    temp = 10.0u"K"
    damping = 1.0u"ps^-1"
    dt = 1.0u"fs"
    n_steps_warmup = 2000
    n_steps = 25_000
    sample_every = 1_000
    outdir = "/mnt/merged/emeitz/tdep_jl_test/group_meeting"

    fcc_crystal = SimpleCrystals.FCC(a, :Ar, SVector(4, 4, 4))
    to_ucposcar(fcc_crystal, joinpath(outdir, "infile.ucposcar"))


    # Configure MD (Molly.jl)
    # Details of this are unimportant, read their docs

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


    # Generate MD Dataset (TDEP.jl)
    sim = NVT(temp, damping, dt, n_steps_warmup, n_steps, sample_every)
    generate_MDTDEP_dataset(sys, sim, outdir; n_seeds = 2)

    # Extract IFCs (TDEP.jl)
    ncores = 10
    efc = ExtractForceConstants(secondorder_cutoff = 4.0)
    execute(efc, outdir, ncores)

end


@testset "TDEP Al MD Extension" begin

    # Settings
    a = 4.05u"Å" 
    temp = 10.0u"K"
    damping = 1.0u"ps^-1"
    dt = 1.0u"fs"
    n_steps_warmup = 2000
    n_steps = 25_000
    sample_every = 1_000
    outdir = "/mnt/merged/emeitz/tdep_jl_test/group_meeting"

    # Create Structure (SimpleCrystals.jl)
    al_crystal = SimpleCrystals.FCC(a, 26.981539u"u", SVector(4, 4, 4))
    to_ucposcar(al_crystal, joinpath(outdir, "infile.ucposcar"))

    # Configure MD (Molly.jl)
    sys = System(
        al_crystal;
        energy_units=u"eV",
        force_units=u"eV * Å^-1",
    )

    ase_eam = pyimport("ase.calculators.eam")
    pot = ase_eam.EAM(potential = "/home/emeitz/software/lammps/potentials/Al_zhou.eam.alloy")

    calc = ASECalculator(
        ase_calc=pot,
        atoms=sys.atoms,
        coords=sys.coords,
        boundary=sys.boundary,
        elements = fill("Al", length(sys.atoms))
    )

    sys = System(
        sys;
        general_inters=(calc,)
    )

    # Generate MD Dataset (TDEP.jl)
    sim = NVT(temperature, damping, dt, n_steps_warmup, n_steps, sample_every)
    generate_MDTDEP_dataset(sys, sim, outdir; n_seeds = 2)

end

@testset "TI Pair Potential" begin

    ifc_path = ""
    n_uc = 4
    a = 5.2468u"Å" # Lattice parameter for FCC Argon at 10 K
    temp = 10.0u"K"
    damping = 1.0u"ps^-1"
    dt = 1.0u"fs"
    n_steps_warmup = 2000
    n_steps = 25_000
    sample_every = 1_000
    n_lambda = 13

    fcc_crystal = SimpleCrystals.FCC(a, :Ar, SVector(n_uc, n_uc, n_uc))
    m = mass(fcc_crystal, 1)

    energy_unit = u"eV"
    length_unit = u"Å"
    ifc2_unit = energy_unit / length_unit^2

    ifc2 = (jldopen(ifc_path, "dynmat") .* ustrip(m)) * ifc2_unit

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

    sim = NVT(temperature, damping, dt, n_steps_warmup, n_steps, sample_every)

    TI(
        sys,
        sim,
        ifc2,
        n_lambda,
    )
end