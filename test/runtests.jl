
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

    using Molly
    using TDEP
    using SimpleCrystals
    using JLD2

    # ifc_path = "/mnt/mntsdb/emeitz/ForceConstants/LJ_ALM/LJ_10K_residual.jld2"
    # ifc_path = raw"C:/Users/ejmei/Desktop/LJ_TDEP_80K_5UC.jld2"
    ifc_path = raw"C:\Users\ejmei\Desktop\LJ_10K_residual.jld2"
    n_uc = 4
    a = 5.2468u"Å" # Lattice parameter for FCC Argon at 10 K
    temp = 10.0u"K"
    damping = 0.5u"ps^-1"
    dt = 2.0u"fs"
    r_cut = 8.5u"Å"
    n_steps_warmup = 50_000
    n_steps = 250_000
    sample_every = 50
    n_lambda = 13
    # U0 = -20.06716u"eV" #! DEPENDS ON NUMBER OF ATOMS and TEMP!
    # U0 = -19.903232u"eV" # 10K
    U0 = 0.0u"eV"

    ifc_conv = 23.060541945 # converts kcal/mol/Ang^2 to eV/Ang^2
    f_conv = sqrt(418.4) * 1e12 # converts freqs from real units --> rad /s 
    # f_conv = 9.82269474855602e13 # converts freqs from metal units --> rad / s
    
    fcc_crystal = SimpleCrystals.FCC(a, :Ar, SVector(n_uc, n_uc, n_uc))
    m = AtomsBase.mass(fcc_crystal, 1)

    energy_unit = u"eV"
    length_unit = u"Å"
    ifc2_unit = energy_unit / length_unit^2

    freqs_sq, dynmat = load(ifc_path, "freqs_sq", "dynmat");
    freqs = (f_conv .* sqrt.(freqs_sq)) * u"rad / s";
    ifc2 = (dynmat .* ustrip(m)  ./ ifc_conv) * ifc2_unit

    pot = LennardJones(; cutoff=ShiftedForceCutoff(r_cut), )

    sys = System(
        fcc_crystal;
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

    sim = NVT(temp, damping, dt, n_steps_warmup, n_steps, sample_every)

    ΔF, F0, mean_ΔUs = TI(
                            sys,
                            pot,
                            sim,
                            ifc2,
                            freqs,
                            U0,
                            n_lambda,
                        )
end


# T = 10.0u"K"
# pot = LennardJones(cutoff=ShiftedForceCutoff(r_cut), )

# sys = System(
#         fcc_crystal;
#         pairwise_inters=(pot,),
#         energy_units=u"eV",
#         force_units=u"eV * Å^-1",
#     )


# # Set sigma and epsilon parameters properly
# σ = 3.4u"Å"
# ϵ = 0.010423u"eV"
# updated_atoms = []

# for i in eachindex(sys)
#     push!(updated_atoms, Molly.Atom(index=sys.atoms[i].index, atom_type=sys.atoms[i].atom_type,
#                             mass=sys.atoms[i].mass, charge=sys.atoms[i].charge,
#                             σ=σ, ϵ=ϵ))
# end

# sys = System(sys; 
#             atoms=[updated_atoms...], 
#             loggers=(disp = DisplacementsLogger(50, sys.coords),
#             traj = XYZWriter(50, "C:/Users/ejmei/Desktop/LJ_traj.xyz"))
# )

# sim = NVT(temp, damping, dt, n_steps_warmup, n_steps, sample_every)

# random_velocities!(sys, T)

# simulate!(sys, sim.thermostat, sim.n_steps_warmup; run_loggers=true)
# simulate!(sys, sim.thermostat, sim.n_steps)


# # True: -0.0820714
# # 4UC : -0.06997523423449314 
# # 5UC : -0.0701189 


# struct XYZWriter
#     n_steps::Int
#     filepath::String
# end


# function Molly.log_property!(logger::XYZWriter, sys::System, neighbors=nothing,
#                        step_n::Integer=0; kwargs...)
#     if step_n % logger.n_steps == 0
#         open(logger.filepath, "a") do file
#             write(file, "$(length(sys.atoms))\n")
#             write(file, "Step: $step_n\n")
#             for atom in sys.atoms
#                 pos = ustrip.(sys.coords[atom.index])
#                 write(file, "$(atom.atom_type) $(pos[1]) $(pos[2]) $(pos[3])\n")
#             end
#         end
#     end
# end