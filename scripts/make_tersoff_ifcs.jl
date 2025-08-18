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

ifc_r_cut = 16.0u"angstrom" # potential has 10.1 Ang cutoff
out_dir = (T) -> "/mnt/mntsdb/emeitz/ForceConstants/AlEAM/$(ustrip(T))K"
a = 4.041u"angstrom" # defined in EAM potential file
m = 26.982u"u"
temps = [100, 200, 300, 400, 500, 600]u"K"
dt = 1.0u"fs"
damping = 0.5u"ps^-1"
n_steps_warmup = 100_000
n_steps = 500_000
sample_every = 10_000
n_seeds = 5

fcc_crystal = FCC(a, :Al, SVector(5,5,5))

sys = System(
        fcc_crystal,
        energy_units = u"eV",
        force_units = u"eV / angstrom"
    )


ase_calc = ase_eam.EAM(potential="/home/emeitz/software/lammps/potentials/Al_zhou.eam.alloy")
inter = ASECalculator(
        ase_calc = ase_calc,
        atoms = sys.atoms,
        coords = sys.coords,
        boundary = sys.boundary,
        atoms_data = sys.atoms_data
    )

sys = System(sys; general_inters = (inter, ))

for T in temps
    @info "Starting $T"
    # @set ntasks = THREADS_PER_CHUNK 
    mkpath(out_dir(T))
    to_ucposcar(fcc_crystal, joinpath(out_dir(T), "infile.ucposcar"))
    sim = NVT(T, damping, dt, n_steps_warmup, n_steps, sample_every)
    generate_MDTDEP_dataset(sys, sim, out_dir(T); n_seeds = n_seeds, n_threads = THREADS_PER_SIM)
    efc = ExtractForceConstants(secondorder_cutoff = ustrip(u"angstrom", ifc_r_cut))
    execute(efc, outdir(T), ncores)
end