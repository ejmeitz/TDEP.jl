export TI, HarmonicPotential

function F_harmonic(freqs, ħ, kB, T)

    kBT = kB * T
    f = (freq) -> (0.5*ħ*freq) + kBT * log(1 - exp(-ħ*freq/kBT))

    F_zero = zero(ħ*first(freqs))
    E_preferred = unit(kBT)

    F₀ = mapreduce(+, freqs) do ω
        if ustrip(ω) > 1e-6
            return f(ω)
        else
            return F_zero
        end
    end

    return uconvert(E_preferred, F₀)

end

# Takes freqs with units of dynmat
# Assumes everything has length_units = u"Å", energy_units =  u"eV"
function TDEP.TI(
        sys::System{3},
        pair_pot,
        sim::NVT,
        ifc2::Matrix,
        freqs,
        lambdas::AbstractVector;
        weights = ones(length(lambdas)) ./ length(lambdas)
    )

    if !all(0.0 .<= lambdas .<= 1.0)
        throw(ArgumentError("Lambdas must be in the range [0, 1]."))
    end

    f = (x) -> TI_core(sys, pair_pot, sim, x, ifc2)

    ET = typeof(1.0u"eV")
    mean_ΔUs = zeros(ET, length(lambdas)) # ⟨ΔU⟩

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

    return ΔF, F0 

end


function TDEP.TI(
        sys::System{3},
        pair_pot,
        sim::NVT,
        ifc2::Matrix,
        freqs,
        n_lambda::Integer,
        quadrature_rule::Function = FastGaussQuadrature.gausslegendre
    )

    λ, w = quadrature_rule(n_lambda)

    # Map to [0,1]
    λ = (λ .+ 1) ./ 2
    w = w ./ 2

    return TI(sys, pair_pot, sim, ifc2, freqs, λ; weights = w)

end

# Implement Molly.jl "General Interaction"
#* TODO JUST CALL libolle.jl --> load IFC obj from path, call pot_eng
struct HarmonicPotential{R,D,T}
    ifc2::Matrix{T}
    r₀::AbstractVector{SVector{D, R}}
    u_storage::AbstractVector{SVector{D, R}}
end

function HarmonicPotential(ifc2, r₀)
    T = eltype(ifc2)
    return HarmonicPotential{eltype(first(r₀)), length(first(r₀)), T}(ifc2, r₀, similar(r₀))
end

# function HarmonicPotential(ifc2_path::String)
#  #TODO
# end

function AtomsCalculators.forces(sys, inter::HarmonicPotential{R,D}; kwargs...) where {R,D}
    inter.u_storage .= sys.coords .- inter.r₀ #* only update in forces
    FT = float_type(sys)
    one_T = one(float_type(sys)) # get 1.0 with proper floating point precision
    F = -one_T .* (inter.ifc2 * reinterpret(R, inter.u_storage))
    return  reinterpret(SVector{D, FT}, F) # return without units for now
end

function AtomsCalculators.potential_energy(sys, inter::HarmonicPotential{R}; kwargs...) where R
    FT = float_type(sys)
    half =  FT(0.5)
    u_1D = reinterpret(R, inter.u_storage)
    return half * ((transpose(u_1D) * inter.ifc2) * u_1D)
end


struct MixedHamiltonian{A, B, C <: HarmonicPotential, U}
    λ::A
    H1::B #* FOR NOW ASSUME ITS A SINGLE PAIR POTENTIAL
    H0::C
    ΔU_storage::U
end

function MixedHamiltonian(λ, H1, H0, n_steps; energy_units = u"eV")
    ΔU_storage = zeros(n_steps+1) * energy_units
    return MixedHamiltonian{typeof(λ), typeof(H1), typeof(H0), typeof(ΔU_storage)}(λ, H1, H0, ΔU_storage)
end

function AtomsCalculators.forces(sys, inter::MixedHamiltonian; 
                                    step_n = 0, n_threads=Threads.nthreads(), kwargs...)

    forces_nounits = Molly.ustrip_vec.(zero(sys.coords))
    forces_buffer = Molly.init_forces_buffer!(sys, forces_nounits, n_threads)

    F0 = AtomsCalculators.forces(sys, inter.H0; kwargs...)

    #* ASSUMES PAIR POTENTIAL
    neighbors = Molly.find_neighbors(sys; n_threads=n_threads)
    pairwise_inters_nonl = filter(!use_neighbors, values(sys.pairwise_inters))
    pairwise_inters_nl   = filter( use_neighbors, values(sys.pairwise_inters))
    if n_threads > 1
        F1 = Molly.pairwise_forces_threads!(forces_nounits, forces_buffer, sys.atoms, sys.coords, sys.velocities,
                                    sys.boundary, neighbors, sys.force_units, length(sys),
                                    pairwise_inters_nonl, pairwise_inters_nl, n_threads, step_n)
    else
        F1 = Molly.pairwise_forces!(forces_nounits, sys.atoms, sys.coords, sys.velocities, sys.boundary,
                            neighbors, sys.force_units, length(sys), pairwise_inters_nonl,
                            pairwise_inters_nl, step_n)
    end


    return ((1-inter.λ) .* (F0 .* sys.force_units)) .+ inter.λ .* (F1 .* sys.force_units)
end

function AtomsCalculators.potential_energy(sys, inter::MixedHamiltonian; 
                                            step_n = 0, n_threads=Threads.nthreads(), kwargs...)
    
    U0 = AtomsCalculators.potential_energy(sys, inter.H0; step_n=step_n, n_threads=n_threads, kwargs...)
    T = typeof(ustrip(zero(eltype(eltype(sys.coords)))))
    #* ASSUMES PAIR POTENTIAL
    neighbors = Molly.find_neighbors(sys; n_threads=n_threads)
    pairwise_inters_nonl = filter(!use_neighbors, values(sys.pairwise_inters))
    pairwise_inters_nl   = filter( use_neighbors, values(sys.pairwise_inters))
    if n_threads > 1
        U1 = Molly.pairwise_pe_threads(sys.atoms, sys.coords, sys.velocities, sys.boundary,
                                    neighbors, sys.energy_units, length(sys), pairwise_inters_nonl,
                                    pairwise_inters_nl, T, n_threads, step_n)
    else
        U1 = Molly.pairwise_pe(sys.atoms, sys.coords, sys.velocities, sys.boundary, neighbors,
                            sys.energy_units, length(sys), pairwise_inters_nonl,
                            pairwise_inters_nl, T, step_n)
    end

    inter.ΔU_storage[step_n+1] = U1 - U0

    return ((1-inter.λ) * U0) + inter.λ*U1
end

# Runs MD simulation to approximate ΔU
function TI_core(sys::System{3}, pair_pot, sim::NVT, λ, ifc2)

    if length(sys.loggers) > 0
        @warn "Found $(length(sys.loggers)) loggers, removing from system"
    end

    if length(sys.pairwise_inters) > 0 || length(sys.general_inters) > 0 || length(sys.specific_inter_lists) > 0
        @warn "System passed with interactions, these will be ignored."
    end

    length_units = u"Å"
    energy_units =  u"eV"

    LT = typeof(1.0*length_units)
    ET = typeof(1.0*energy_units)

    new_loggers = (
        # coords = CoordinatesLogger(LT, sim.sample_every), 
        pe = PotentialEnergyLogger(ET, sim.sample_every), 
    )

    #* HOW TO LOAD IFCS
    hp = HarmonicPotential(ifc2, copy(sys.coords))

    mixed_H = MixedHamiltonian(
        λ, 
        pair_pot,
        hp,
        sim.n_steps;
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

    run_sim!(new_system, sim)
    
    return mean(mixed_H.ΔU_storage)

end