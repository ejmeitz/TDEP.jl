export TI, HarmonicPotential

function TDEP.TI(
        sys::System{3},
        pair_pot,
        sim::NVT,
        ifc2::Matrix,
        lambdas::AbstractVector;
        weights = ones(length(lambdas)) ./ length(lambdas)
    )

    f = (x) -> TI_core(sys, pair_pot, sim, x, ifc2)
    return dot(weights, f.(lambdas))

end


function TDEP.TI(
        sys::System{3},
        pair_pot,
        sim::NVT,
        ifc2::Matrix,
        n_lambda::Integer,
        quadrature_rule::Function = FastGaussQuadrature.gausslegendre
    )

    λ, w = quadrature_rule(n_lambda)
    return TI(sys, pair_pot, sim, ifc2, λ; weights = w)

end

# Implement Molly.jl "General Interaction"
#* TODO JUST CALL libolle.jl --> load IFC obj from path, call pot_eng
struct HarmonicPotential{T, R}
    ifc2::Matrix{T}
    r₀::R
    u_storage::R
end

function HarmonicPotential(ifc2, r₀)
    T = eltype(ifc2)
    return HarmonicPotential{T, typeof(r₀)}(ifc2, r₀, similar(r₀))
end

# function HarmonicPotential(ifc2_path::String)
#  #TODO
# end

function AtomsCalculators.forces(sys, inter::HarmonicPotential; kwargs...)
    inter.u_storage .= sys.coords .- inter.r₀
    one_T = one(float_type(sys)) # get 1.0 with proper floating point precision
    return -one_T .* (inter.ifc2 * u_storage)
end

function AtomsCalculators.potential_energy(sys, inter::HarmonicPotential; kwargs...)
    FT = float_type(sys)
    half =  FT(0.5)
    return half * ((transpose(inter.u_storage) * ifc2) * inter.u_storage)
end


struct MixedHamiltonian{A, B, C <: HarmonicPotential, U}
    λ::A
    H1::B #* FOR NOW ASSUME ITS A SINGLE PAIR POTENTIAL
    H0::C
    ΔU_storage::U
end

function MixedHamiltonian(λ, H1, H0, n_steps; energy_units = u"eV")
    ΔU_storage = zeros(n_steps) * energy_units
    return MixedHamiltonian{typeof(λ), typeof(H1), typeof(H0), typeof(ΔU_storage)}(λ, H1, H0, ΔU_storage)
end

function AtomsCalculators.forces(sys, inter::MixedHamiltonian; 
                                    n_threads=Threads.nthreads(), kwargs...)

    F0 = AtomsCalculators.forces(sys, inter.H0; kwargs...)

    #* ASSUMES PAIR POTENTIAL
    neighbors = find_neighbors(sys; n_threads=n_threads)
    pairwise_inters_nonl = filter(!use_neighbors, values(sys.pairwise_inters))
    pairwise_inters_nl   = filter( use_neighbors, values(sys.pairwise_inters))
    if n_threads > 1
        pairwise_forces_threads!(fs_nounits, fs_chunks, sys.atoms, sys.coords, sys.velocities,
                                    sys.boundary, neighbors, sys.force_units, length(sys),
                                    pairwise_inters_nonl, pairwise_inters_nl, n_threads, step_n)
    else
        pairwise_forces!(fs_nounits, sys.atoms, sys.coords, sys.velocities, sys.boundary,
                            neighbors, sys.force_units, length(sys), pairwise_inters_nonl,
                            pairwise_inters_nl, step_n)
    end

    return ((1-inter.λ) * F0) + inter.λ*F1
end

function AtomsCalculators.potential_energy(sys, inter::MixedHamiltonian; 
                                            step_n = 0, n_threads=Threads.nthreads(), kwargs...)
    
    U0 = AtomsCalculators.potential_energy(sys, inter.H0; kwargs...)

    #* ASSUMES PAIR POTENTIAL
    neighbors = find_neighbors(sys; n_threads=n_threads)
    pairwise_inters_nonl = filter(!use_neighbors, values(sys.pairwise_inters))
    pairwise_inters_nl   = filter( use_neighbors, values(sys.pairwise_inters))
    if n_threads > 1
        U1 = pairwise_pe_threads(sys.atoms, sys.coords, sys.velocities, sys.boundary,
                                    neighbors, sys.energy_units, length(sys), pairwise_inters_nonl,
                                    pairwise_inters_nl, T, n_threads, step_n)
    else
        U1 = pairwise_pe(sys.atoms, sys.coords, sys.velocities, sys.boundary, neighbors,
                            sys.energy_units, length(sys), pairwise_inters_nonl,
                            pairwise_inters_nl, T, step_n)
    end

    inter.ΔU_storage[step_n] = U1 - U0

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

    return mixed_H.ΔU_storage

end