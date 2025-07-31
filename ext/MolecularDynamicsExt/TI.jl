export TI, HarmonicPotential

# classical
function F_harmonic(freqs, ħ, kB, T)

    kBT = kB * T
    f = (freq) -> log(ħ*freq/kBT)

    E_preferred = unit(kBT)

    F₀ = mapreduce(+, freqs) do ω
        if ustrip(ω) > 1e-6
            return f(ω)
        else
            return 0.0
        end
    end

    return uconvert(E_preferred, kBT * F₀)

end

# Takes freqs with units of dynmat
# Assumes everything has length_units = u"Å", energy_units =  u"eV"
function TDEP.TI(
        sys::System{3},
        pair_pot,
        sim::NVT,
        ifc2::Matrix,
        freqs,
        U0,
        lambdas;
        weights,
        nthreads = Threads.nthreads(),
    )

    if !all(0.0 .<= lambdas .<= 1.0)
        throw(ArgumentError("Lambdas must be in the range [0, 1]."))
    end

    ET = typeof(1.0u"eV")
    mean_ΔUs = zeros(ET, length(lambdas)) # ⟨ΔU⟩

    f = (x) -> TI_core(sys, pair_pot, sim, x, ifc2, U0; nthreads = nthreads)

    p = Progress(length(lambdas); dt=1.0)
    for (i, λ) in enumerate(lambdas)
        mean_ΔUs[i] = f(λ)
        next!(p)
    end
    finish!(p)
    
    ΔF = dot(weights, mean_ΔUs)

    freqs_rad_s = nothing
    try
        freqs_rad_s = uconvert.(u"rad / s", freqs)
    catch
        throw(ArgumentError("Could not convert freuqencies to rad / s. This usually happens when the dynamical matrix has molar units. Try multiplying freqs_sq by Unitful.Na then taking sqrt."))
    end

    # Calculate harmonic F
    ħ = uconvert(u"eV * s", Unitful.ħ) # Planck's constant
    kB = sys.k # Boltzmann constant
    F0 = F_harmonic(freqs_rad_s, ħ, kB, temperature(sim))

    return ΔF, F0, mean_ΔUs

end


function TDEP.TI(
        sys::System{3},
        pair_pot,
        sim::NVT,
        ifc2::Matrix,
        freqs,
        U0,
        n_lambda::Integer;
        nthreads = Threads.nthreads(),
        quadrature_rule::Function = FastGaussQuadrature.gausslegendre
    )

    λ, w = quadrature_rule(n_lambda)

    if !all(-1.0 .<= λ .<= 1.0)
        throw(ArgumentError("Quadrature rule must emit points in [-1, 1]."))
    end

    # Map to [0,1]
    λ = (λ .+ 1) ./ 2
    w ./= 2
    # println(λ)
    # return

    return TI(sys, pair_pot, sim, ifc2, freqs, U0, λ; weights = w, nthreads = nthreads)

end

# Implement Molly.jl "General Interaction"
#* TODO JUST CALL libolle.jl --> load IFC obj from path, call pot_eng
struct HarmonicPotential{R,D,T,U,F}
    ifc2::Matrix{T}
    force_buffer::F
    U₀::U
end

function HarmonicPotential(ifc2, U₀, coord_T::DataType, force_buffer; dims = 3)
    T = eltype(ifc2)
    return HarmonicPotential{coord_T, dims, T, typeof(U₀), typeof(force_buffer)}(ifc2, force_buffer, U₀)
end

# function HarmonicPotential(ifc2_path::String)
#  #TODO
# end

function AtomsCalculators.forces(sys, inter::HarmonicPotential{R,D}; kwargs...) where {R,D}
    FT = float_type(sys)
    inter.force_buffer .= ustrip_vec.(-(inter.ifc2 * reinterpret(R, sys.loggers.disp.last_displacements)))
    return  reinterpret(SVector{D, FT}, inter.force_buffer) # return without units for now
end

function AtomsCalculators.potential_energy(sys, inter::HarmonicPotential{R}; kwargs...) where R
    FT = float_type(sys)
    half =  FT(0.5)
    #* THIS ASSUMES PE LOGGER RUN AFTER DISP LOGGER
    u_1D = reinterpret(R, sys.loggers.disp.last_displacements)
    return inter.U₀ + (half * dot(u_1D, inter.ifc2 * u_1D))
end


mutable struct MixedHamiltonian{A, B, C <: HarmonicPotential, U}
    const λ::A
    const H1::B #* FOR NOW ASSUME ITS A SINGLE PAIR POTENTIAL
    const H0::C
    const ΔU_storage::U
    current_step::Int
end

function MixedHamiltonian(λ, H1, H0, n_steps, sampling_rate; energy_units = u"eV")
    ΔU_storage = zeros(div(n_steps, sampling_rate) + 2) * energy_units
    return MixedHamiltonian{typeof(λ), typeof(H1), typeof(H0), typeof(ΔU_storage)}(λ, H1, H0, ΔU_storage, 1)
end

function AtomsCalculators.forces(sys, inter::MixedHamiltonian; 
                                    step_n = 0, n_threads=Threads.nthreads(), kwargs...)

    forces_nounits = Molly.ustrip_vec.(zero(sys.coords))
    forces_buffer = Molly.init_forces_buffer!(sys, forces_nounits, n_threads)

    F0 = AtomsCalculators.forces(sys, inter.H0; kwargs...)

    #* ASSUMES PAIR POTENTIAL
    neighbors = Molly.find_neighbors(sys; n_threads=n_threads)
    pairwise_inters_nonl = filter(!use_neighbors, (inter.H1,))
    pairwise_inters_nl   = filter( use_neighbors, (inter.H1,))
    if n_threads > 1
        F1 = Molly.pairwise_forces_threads!(forces_nounits, forces_buffer, sys.atoms, sys.coords, sys.velocities,
                                    sys.boundary, neighbors, sys.force_units, length(sys),
                                    pairwise_inters_nonl, pairwise_inters_nl, n_threads, step_n)
    else
        F1 = Molly.pairwise_forces!(forces_nounits, sys.atoms, sys.coords, sys.velocities, sys.boundary,
                            neighbors, sys.force_units, length(sys), pairwise_inters_nonl,
                            pairwise_inters_nl, step_n)
    end


    return ((1-inter.λ) .* (F0 .* sys.force_units)) .+ (inter.λ .* (F1 .* sys.force_units))
end

function AtomsCalculators.potential_energy(sys, inter::MixedHamiltonian; 
                                            step_n = 0, n_threads=Threads.nthreads(), kwargs...)
    
    U0 = AtomsCalculators.potential_energy(sys, inter.H0; step_n=step_n, n_threads=n_threads, kwargs...)

    T = typeof(ustrip(zero(eltype(eltype(sys.coords)))))
    #* ASSUMES PAIR POTENTIAL
    neighbors = Molly.find_neighbors(sys; n_threads=n_threads)
    pairwise_inters_nonl = filter(!use_neighbors, (inter.H1,))
    pairwise_inters_nl   = filter( use_neighbors, (inter.H1,))
    if n_threads > 1
        U1 = Molly.pairwise_pe_threads(sys.atoms, sys.coords, sys.velocities, sys.boundary,
                                    neighbors, sys.energy_units, length(sys), pairwise_inters_nonl,
                                    pairwise_inters_nl, T, n_threads, step_n)
    else
        U1 = Molly.pairwise_pe(sys.atoms, sys.coords, sys.velocities, sys.boundary, neighbors,
                            sys.energy_units, length(sys), pairwise_inters_nonl,
                            pairwise_inters_nl, T, step_n)
    end

    inter.ΔU_storage[inter.current_step] = U1 - U0
    inter.current_step += 1

    return ((1-inter.λ) * U0) + (inter.λ*U1)
end

# Runs MD simulation to approximate ΔU
function TI_core(sys::System{D}, pair_pot, sim::NVT, λ, ifc2, U0; nthreads = Threads.nthreads()) where D

    if length(sys.loggers) > 0
        @warn "Found $(length(sys.loggers)) loggers, removing from system"
    end

    if length(sys.pairwise_inters) > 0 || length(sys.general_inters) > 0 || length(sys.specific_inter_lists) > 0
        @warn "System passed with interactions, these will be ignored."
    end

    length_units = u"Å"
    energy_units =  u"eV"

    ET = typeof(1.0*energy_units)
    LT = typeof(1.0*length_units)

    # The order of these matters! We use displacements
    # to calculate PE, so that must be first.
    new_loggers = (
        # we don't actually need these saved, but need them available during sim
        disp = DisplacementsLogger(100_000, sys.coords),
        pe = PotentialEnergyLogger(ET, sim.sample_every),
    )

    coord_T = eltype(first(sys.coords))
    force_buffer = zeros(D*length(sys))
    hp = HarmonicPotential(ifc2, U0, coord_T, force_buffer; dims = D)

    mixed_H = MixedHamiltonian(
        λ, 
        pair_pot,
        hp,
        sim.n_steps + sim.n_steps_warmup,
        sim.sample_every; 
        energy_units = energy_units
    )

    new_system = System(
        deepcopy(sys);
        loggers = deepcopy(new_loggers),
        general_inters = (mixed_H,),
        force_units = energy_units / length_units,
        energy_units = energy_units,
        k = Molly.default_k(energy_units)
    )

    run_sim!(new_system, sim; run_loggers_warmup = true, nthreads = nthreads)

    n_samples_warmup = div(sim.n_steps_warmup, sim.sample_every) + 1

    println()
    println("λ = $(mixed_H.λ)")
    println("Max Dispalcement, ", maximum(norm.(new_system.loggers.disp.last_displacements)))
    println()

    return @views mean(mixed_H.ΔU_storage[n_samples_warmup : end])

end 