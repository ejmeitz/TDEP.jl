export TI, HarmonicPotential

function TI(
        sys::System{3},
        sim::NVT,
        lambdas::AbstractVector;
        weights = ones(length(lambdas)) ./ length(lambdas)
    )

    f = (x) -> TI_core(sys, sim, x)
    return dot(weights, f.(lambdas))

end


function TI(
        sys::System{3},
        sim::NVT,
        n_lambda::Integer,
        quadrature_rule::Function = FastGaussQuadrature.gausslegendre
    )

    λ, w = quadrature_rule(n_lambda)
    return TI(sys, sim, λ; weights = w)

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
    return HarmonicPotential{T}(ifc2, r₀, similar(r₀))
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
    H1::B #* DOES THIS HAVE TO BE A GENERAL INTERACTION?? SUPPORT MULTIPLE INTERACTIONS
    H0::C
    ΔU_storage::U
end

function AtomsCalculators.forces(sys, inter::MixedHamiltonian; kwargs...)

    F1 = ...
    F0 = AtomsCalculators.forces(sys, inter.H0; kwargs...)

    return ((1-inter.λ) * F0) + inter.λ*F1
end

function AtomsCalculators.potential_energy(sys, inter::MixedHamiltonian; step_n = 0, kwargs...)
    
    U1 = ...
    U0 = AtomsCalculators.potential_energy(sys, inter.H0; kwargs...)

    inter.ΔU_storage[step_n] = U1 - U0

    return ((1-inter.λ) * U0) + inter.λ*U1

end


#* HAVE USER DEFINE "REAL" POTENTIAL IN sys OBJECT, I ADD THE HARMONIC POTENTIAL IN 
#* BUT YOU CANT APPLY WEIGHTS THSI WAY

# Runs MD simulation to approximate ΔU
function TI_core(sys::System{3}, sim::NVT, λ)

    if length(sys.loggers) > 0
        @warn "Found $(length(sys.loggers)) loggers, removing from system"
    end

    length_units = u"Å"
    energy_units =  u"eV"

    LT = typeof(1.0*length_units)
    ET = typeof(1.0*energy_units)

    new_loggers = (
        coords = CoordinatesLogger(LT, sim.sample_every), 
        pe = PotentialEnergyLogger(ET, sim.sample_every), 
    )

    #* HOW TO LOAD IFCS
    hp = HarmonicPotential(ifc2, copy(sys.coords))

    mixed_H = MixedHamiltonian(
        λ, 
        #* HOW,
        hp
    )


    new_system = System(
        deepcopy(sys);
        loggers = deepcopy(new_loggers),
        force_units = energy_units / length_units,
        energy_units = energy_units,
        k = Molly.default_k(energy_units)
    )

    run_sim!(new_system, sim)



end