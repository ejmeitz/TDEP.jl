using ASEconvert
using PythonCall
using Molly
using TDEP
using SimpleCrystals
using OhMyThreads
ase_eam = pyimport("ase.calculators.eam")

const THREADS_PER_SIM = 4
const NUM_THREADS_TOTAL = 40
const THREADS_PER_CHUNK = Int(NUM_THREADS_TOTAL / THREADS_PER_SIM)

ifc_r_cut = 11.0 # potential has 6 Ang cutoff, this should be safe
out_dir = (T) -> "/mnt/mntsdb/emeitz/ForceConstants/AlEAM/$(T)K"
a = 3.9860u"angstrom" # defined in EAM potential file
m = 26.982u"u"
temps = [100, 200, 300, 400, 500, 600]u"K"
dt = 1.0u"fs"
damping = 0.5u"ps^-1"
n_steps_warmup = 100_000
n_steps = 2_500_000
sample_every = 5000
n_seeds = 5

diamond_crystal = Diamond(a, :Al, SVector(4,4,4))

ase_calc = ase_eam.EAM("/home/emeitz/software/lammps/potentials/Al_jnp.eam")

sys = System(
            diamond_crystal;
            general_inters = (ASEcalculator(ase_calc), ),
            energy_units = u"eV",
            force_units = u"eV / angstrom"
        )

@tasks for T in temps
    @info "Starting $T"
    @set ntasks = THREADS_PER_CHUNK 
    mkpath(out_dir(T))
    to_ucposcar(diamond_crystal, joinpath(out_dir(T), "infile.ucposcar"))
    sim = NVT(T, damping, dt, n_steps_warmup, n_steps, sample_every)
    generate_MDTDEP_dataset(sys, sim, out_dir(T); n_seeds = n_seeds, n_threads = THREADS_PER_SIM)
    efc = ExtractForceConstants(secondorder_cutoff = ifc_r_cut)
    execute(efc, outdir(T), ncores)
end