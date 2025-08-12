using Molly
using TDEP
using SimpleCrystals
using JLD2
using Plots
using DelimitedFiles
using ProgressMeter
using OhMyThreads

const N_THREADS_TOTAL = 40
const N_THREADS_PER_SIM = 8
const N_THREADS_PER_TASK = Int(N_THREADS_TOTAL / N_THREADS_PER_SIM)

function lambda_study(sys::System, sim::NVT, pair_pot, ifc2, freqs, U0, n_lambdas::Vector{Int}, outpath::String)
    F = []
    delF_arr = []
    F0_arr = []
    delU_data = Dict()
    λs = Dict()
    # p = Progress(length(n_lambdas); dt = 5, desc = "Lambda Study:")
    @tasks for N in n_lambdas
        @info "Starting $N"
        @set ntasks = N_THREADS_PER_TASK
        delF, F0, delU_data[N], λs[N] = TI(sys, pair_pot, sim, ifc2, freqs, U0, N; nthreads = N_THREADS_PER_SIM)
        push!(F, delF + F0 + U0)
        push!(delF_arr, delF)
        push!(F0_arr, F0 + U0) # F_TDEP = F_harmonic + <V - V2>
        # next!(p)
    end
    # finish!(p)

    for N in n_lambdas
        s1 = scatter(λs[N], ustrip.(delU_data[N]) ./ length(sys); xlabel = "Lambda", ylabel = "<ΔU> [eV / atom]")
        savefig(s1, joinpath(outpath, "lambda$(N).png"))
    end

    header = ["n_lambda" "F0 + U0" "ΔF" "F"]
    open(joinpath(outpath, "F.txt"), "w") do f
        writedlm(f, [header; n_lambdas ustrip.(F0_arr) ustrip.(delF_arr) ustrip.(F)], ",")
    end
    s2 = scatter(n_lambdas, F ./ length(sys), xlabel = "Number of Quadrature Points", ylabel = "F [eV / atom]")
    savefig(s2, joinpath(outpath, "F_vs_Nlambda.png"))
end

function N_steps_study(sys::System, sims::Vector{NVT}, pair_pot, ifc2, freqs, U0, n_lambda, outpath::String)

    F = []
    delF_arr = []
    F0_arr = []
    delU_data = Dict()
    λs = Dict()
    # p = Progress(length(n_lambdas); dt = 5, desc = "Lambda Study:")
    @tasks for sim in sims
        NSTEPS = sim.n_steps
        @info "Starting $NSTEPS"
        @set ntasks = N_THREADS_PER_TASK
        delF, F0, delU_data[NSTEPS], λs[NSTEPS] = TI(sys, pair_pot, sim, ifc2, freqs, U0, n_lambda; nthreads = N_THREADS_PER_SIM)
        push!(F, delF + F0 + U0)
        push!(delF_arr, delF)
        push!(F0_arr, F0 + U0) # F_TDEP = F_harmonic + <V - V2>
        # next!(p)
    end

    all_n_steps = [sim.n_steps for sim in sims]
    SAMPLE_EVERY = first(sims).sample_every

    println("THE FREE ENERGIES ARE: $(F)")

    header = ["n_steps" "F0 + U0" "ΔF" "F"]
    open(joinpath(outpath, "F.txt"), "w") do f
        writedlm(f, [header; all_n_steps ustrip.(F0_arr) ustrip.(delF_arr) ustrip.(F)], ",")
    end
    s1 = scatter(all_n_steps, F ./ length(sys), xlabel = "Number of Steps (Sampled every $(SAMPLE_EVERY))", ylabel = "F [eV / atom]")
    savefig(s1, joinpath(outpath, "F_vs_nsteps.png"))

end

function size_effects_study(systems::Vector{<:System}, n_uc::Vector{Int}, sim::NVT, pair_pots, ifc2s, freqs, U0_per_atom, n_lambda, outpath::String)
    F = []
    delF_arr = []
    F0_arr = []
    # p = Progress(length(n_lambdas); dt = 5, desc = "Lambda Study:")
    @tasks for (i,sys) in collect(enumerate(systems))
        @info "Starting sys with $(length(sys)) atoms"
        @set ntasks = N_THREADS_PER_TASK
        U0 = U0_per_atom*length(sys)
        delF, F0, _, _ = TI(sys, pair_pots[i], sim, ifc2s[i], freqs[i], U0, n_lambda; nthreads = N_THREADS_PER_SIM)
        push!(F, delF + F0 + U0)
        push!(delF_arr, delF)
        push!(F0_arr, F0 + U0) # F_TDEP = F_harmonic + <V - V2>
        # next!(p)
    end

    println("THE FREE ENERGIES ARE: $(F)")

    n_atoms = length.(systems)

    println("THE FREE ENERGIES PER ATOM ARE: $(F ./ n_atoms)")

    header = ["n_atoms" "F0 + U0" "ΔF" "F"]
    open(joinpath(outpath, "F.txt"), "w") do f
        writedlm(f, [header; n_atoms ustrip.(F0_arr) ustrip.(delF_arr) ustrip.(F)], ",")
    end
    s1 = scatter(n_uc, F ./ n_atoms, xlabel = "Number of Unitcells (N x N x N)", ylabel = "F [eV / atom]")
    savefig(s1, joinpath(outpath, "F_vs_nuc.png"))
end

function init_LJ_system(n_uc; a = 5.2468u"Å", r_cut = 8.5u"Å")
    fcc_crystal = SimpleCrystals.FCC(a, :Ar, SVector(n_uc, n_uc, n_uc))
    pot = LennardJones(; cutoff=ShiftedPotentialCutoff(r_cut))

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

    return sys, pot

end

function setup_ifcs(path::String, m, f_conv, ifc_conv, ifc2_unit)
    freqs_sq, dynmat = load(path, "freqs_sq", "dynmat");
    freqs = (f_conv .* sqrt.(freqs_sq)) * u"rad / s";
    ifc2 = (dynmat .* ustrip(m)  ./ ifc_conv) * ifc2_unit
    return freqs, ifc2
end

# Global Settings
damping = 0.5u"ps^-1"
n_steps_warmup = 100_000
n_steps = 1_000_000
n_steps_arr = [50_000, 100_000, 250_000, 500_000, 1_000_000]
n_steps_optimal = 750_000
sample_every = 50
n_lambdas = collect(range(5, 15, step = 2))
n_lambda_optimal = 9 # probably can go lower, this feels safe though
energy_unit = u"eV"
length_unit = u"Å"
ifc2_unit = energy_unit / length_unit^2
# outpath = (T) -> "/mnt/merged/emeitz/TDEP_TI_Test/$(ustrip(T))K/lambda_study"
# outpath = (T) -> "/mnt/merged/emeitz/TDEP_TI_Test/$(ustrip(T))K/nsteps_study"
outpath = (T) -> "/mnt/merged/emeitz/TDEP_TI_Test/$(ustrip(T))K/size_study"

#LJ Settings
ifc_path = (n,T) -> "/mnt/mntsdb/emeitz/ForceConstants/LJ_ALM/LJ_TDEP_$(n)UC_$(ustrip(T))K.jld2"
temps = [10, 80]u"K"
n_uc = [4,5,6]
dt = 2.0u"fs"
U0s = [-19.903232, -20.06716]u"eV" # FOR 4 UNIT CELLS
U0s_per_atom = U0s ./ 256
ifc_conv = 23.060541945 # converts kcal/mol/Ang^2 to eV/Ang^2
f_conv = sqrt(418.4) * 1e12 # converts freqs from real units --> rad /s 
m = 39.95u"u"

# LAMBDA STUDY
# for (i,T) in enumerate(temps)
#     @info "T = $(T)"
#     sys, pot = init_LJ_system(n_uc[1])
#     sim = NVT(T, damping, dt, n_steps_warmup, n_steps, sample_every)
#     freqs, ifc2 = setup_ifcs(ifc_path(n_uc[1], T), sys.masses[1], f_conv, ifc_conv, ifc2_unit)
#     OP = outpath(T)
#     mkpath(OP)
#     lambda_study(sys, sim, pot, ifc2, freqs, U0s[i], n_lambdas, OP)
# end


# NSTEPS STUDY
# for (i,T) in enumerate(temps)
#     @info "T = $(T)"
#     sys, pot = init_LJ_system(n_uc[1])
#     freqs, ifc2 = setup_ifcs(ifc_path(n_uc[1], T), sys.masses[1], f_conv, ifc_conv, ifc2_unit)
#     OP = outpath(T)
#     mkpath(OP)
#     init_nvt_n_steps = (N) -> NVT(T, damping, dt, n_steps_warmup, N, sample_every)
#     sims = init_nvt_n_steps.(n_steps_arr)
#     N_steps_study(sys, sims, pot, ifc2, freqs, U0s[i], n_lambda_optimal, OP)
# end

# SIZE EFFECT STUDY
for (i,T) in enumerate(temps)
    @info "T = $(T)"
    res = init_LJ_system.(n_uc)
    systems = [r[1] for r in res]
    pots = [r[2] for r in res]
    ifc_paths = ifc_path.(n_uc, Ref(T))
    res2 = setup_ifcs.(ifc_paths, Ref(m), Ref(f_conv), Ref(ifc_conv), Ref(ifc2_unit))
    all_freqs = [r[1] for r in res2]
    all_ifc2 = [r[2] for r in res2]
    OP = outpath(T)
    mkpath(OP)
    sim = NVT(T, damping, dt, n_steps_warmup, n_steps_optimal, sample_every)
    size_effects_study(systems, n_uc, sim, pots, all_ifc2, all_freqs, U0s_per_atom[i], n_lambda_optimal, OP)
end

########################################################
########################################################
########################################################


# SW Settings
# temps = [100, 1300]u"K"
# n_uc = [3,4,5]
# dt = 1.0u"fs"
# U0s = []u"eV"
# ifc_conv = 1.0
# f_conv = 9.82269474855602e13 # converts freqs from metal units --> rad / s
# m = 28.085u"u"

