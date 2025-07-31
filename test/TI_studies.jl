using Molly
using TDEP
using SimpleCrystals
using JLD2
using Plots
using DelimitedFiles
using ProgressMeter
using OhMyThreads

const N_THREADS_TOTAL = 40
const N_THREADS_PER_SIM = 4
const N_THREADS_PER_TASK = Int(N_THREADS_TOTAL / N_THREADS_PER_SIM)

function lambda_study(sys::System, sim::NVT, pair_pot, ifc2, freqs, U0, n_lambdas::Vector{Int}, outpath::String)
    F = []
    delU_data = Dict()
    # p = Progress(length(n_lambdas); dt = 5, desc = "Lambda Study:")
    @tasks for N in n_lambdas
        @info "Starting $N"
        @set ntasks = N_THREADS_PER_TASK
        delF, F0, mean_delUs = TI(sys, pair_pot, sim, ifc2, freqs, U0, N; nthreads = N_THREADS_PER_SIM)
        push!(F, delF + F0)
        delU_data[N] = mean_delUs
        # next!(p)
    end
    # finish!(p)

    for N in n_lambdas
        s1 = scatter(collect(range(1,N)), ustrip.(delU_data[N]) ./ length(sys); xlabel = "Lambda", ylabel = "<ΔU> [eV / atom]")
        savefig(s1, joinpath(outpath, "lambda$(N).png"))
    end

    open(joinpath(outpath, "F.txt"), "w") do f
        writedlm(f, F, ",")
    end
    s2 = scatter(n_lambdas, F ./ length(sys), xlabel = "Number of Quadrature Points", ylabel = "F [eV / atom]")
    savefig(s2, joinpath(outpath, "F_vs_Nlambda.png"))
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

# ifc_path = "/mnt/mntsdb/emeitz/ForceConstants/LJ_ALM/LJ_10K_residual.jld2"
# ifc_path = raw"C:/Users/ejmei/Desktop/LJ_TDEP_80K_5UC.jld2"
# ifc_path = raw"C:\Users\ejmei\Desktop\LJ_10K_residual.jld2"



# Global Settings
damping = 0.5u"ps^-1"
n_steps_warmup = 100_000
n_steps = 1_000_000
sample_every = 50
n_lambdas = collect(range(5, 15, step = 2))
energy_unit = u"eV"
length_unit = u"Å"
ifc2_unit = energy_unit / length_unit^2
outpath = (T) -> "/mnt/merged/emeitz/TDEP_TI_Test/$(ustrip(T))K/lambda_study"


#LJ Settings
ifc_path = (n,T) -> "/mnt/mntsdb/emeitz/ForceConstants/LJ_ALM/LJ_TDEP_$(n)UC_$(ustrip(T))K.jld2"
temps = [10, 80]u"K"
n_uc = [4,5,6]
dt = 2.0u"fs"
U0s = [-19.903232, -20.06716]u"eV"
ifc_conv = 23.060541945 # converts kcal/mol/Ang^2 to eV/Ang^2
f_conv = sqrt(418.4) * 1e12 # converts freqs from real units --> rad /s 

# LAMBDA STUDY
# CHOOSE LARGE SIM LENGTHS AND LARGEST
# SYSTEM SIZE TO MINIMIZE THEIR EFFECTS
for (i,T) in enumerate(temps)
    @info "T = $(T)"
    sys, pot = init_LJ_system(n_uc[1])
    sim = NVT(T, damping, dt, n_steps_warmup, n_steps, sample_every)
    freqs, ifc2 = setup_ifcs(ifc_path(n_uc[1], T), sys.masses[1], f_conv, ifc_conv, ifc2_unit)
    OP = outpath(T)
    mkpath(OP)
    lambda_study(sys, sim, pot, ifc2, freqs, U0s[i], n_lambdas, OP)
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

