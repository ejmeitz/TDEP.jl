module MolecularDynamicsExt

using TDEP
using Molly
using Unitful
using Printf

"""
Runs a molecular dynamics simulation with [Molly.jl](https://github.com/JuliaMolSim/Molly.jl) given a Molly `System` object. This system could
be generated from any other library compatible with the [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl) interface. For TDEP a common pattern would be to use the native [SimpleCrystals.jl](https://github.com/ejmeitz/SimpleCrystals.jl) library or [ASEconert.jl](https://github.com/mfherbst/ASEconvert.jl) to
generate a AbstractSystem object which can be converted to a Molly System object. 

The simulation is equilibrated for `n_steps_warmup` and then run for `n_steps` with the `sim` simulator (must be Langevin). If the `n_seeds` key-word argument
is passed `n_seeds` simulations are run each initialized with different velocities. Data is saved every `sample_every` steps
during the main run. 

The required data loggers are automatically attached to the System. You only need to define the structure and interaction. 

This command will save the infile.positions, infile.forces, infile.stat, infile.meta and infile.ssposcar to `outdir`. Only the
infile.ucposcar must be made to extract force constants.
"""
function TDEP.generate_MDTDEP_dataset(sys::System{3}, sim, n_steps_warmup::Integer, n_steps::Integer, 
                                    sample_every::Integer, outdir::String; n_seeds::Integer = 1)

    @assert Molly.n_infinite_dims(sys) == 0 "Box must be fully periodic"

    # Warn if stripping existing loggers
    if length(sys.loggers) > 0
        @warn "Found $(length(sys.loggers)) loggers, removing from system"
    end

    if typeof(sim) <: Langevin
        #* MAKE ERROR()
        @error "Simulator must be Langevin. If Molly has a new NVT simulator, open a PR."
    end

    # By default this is 1, just to be safe check its on
    if sim.remove_CM_motion == 0
        @warn "You are not removing COM motion. You should probably set remove_CM_motion in the Langevin simulator."
    end

    length_units = u"Å"
    energy_units =  u"eV"
    force_units = energy_units / length_units

    FT = typeof(1.0*force_units)
    LT = typeof(1.0*length_units)
    ET = typeof(1.0*energy_units)
    TT = typeof(sim.temperature)

    new_loggers = (
        forces = ForcesLogger(FT, sample_every), 
        coords = CoordinatesLogger(LT, sample_every), 
        pe = PotentialEnergyLogger(ET, sample_every), 
        ke = KineticEnergyLogger(ET, sample_every),
        temps = TemperatureLogger(TT, sample_every)
    )

    #* TODO CHECK THIS IS RIGHT MATH, Off by 1?
    n_samples = n_seeds * div(n_steps, sample_every)
    dt_fs = unit(sim.dt) == Unitful.NoUnits ? sim.dt : uconvert(u"fs", sim.dt)
    n_atoms = length(sys)
    
    posns_file = joinpath(outdir, "infile.positions")
    forces_file = joinpath(outdir, "infile.forces")

    # Write ssposcar before we do dynamics
    cv = hcat(Molly.cell_vectors(sys)...) # matrix with vec as cols
    write_ssposcar(outdir, cv, sys.coords, Molly.atomic_symbol(sys))
    write_meta(temperature, n_samples, dt_fs, n_atoms)

    for s in 1:n_seeds
        
        @info "Starting seed $s"

        new_system = System(
            sys;
            loggers = copy(new_loggers),
            force_units = u"eV * Å^-1",
            energy_units = u"eV"
        )
    
        random_velocities!(new_system, sim.temperature)

        simulate!(sys, sim, n_steps_warmup; run_loggers=false)
        simulate!(sys, sim, n_steps)

        @info "Seed $s complete. Writing to disk."

        open(posns_file, "a") do pf
            c = values(sys.loggers.coords)
            for i in eachindex(c)
                @printf pf "%.15f %.15f %.15f" ustrip.(c[i])...
            end
        end

        open(forces_file, "a") do ff
            f = values(sys.loggers.forces)
            for i in eachindex(c)
                @printf ff "%.15f %.15f %.15f" ustrip.(f[i])...
            end
        end

        write_partial_stat(
            outdir, 
            values(sys.loggers.pe),
            values(sys.loggers.ke),
            values(sys.loggers.temps);
            sampled_every = sample_every,
            dt_fs = dt_fs, 
            file_mode = "a"
        )

    end

end

function generate_STDEP_dataset(sys::System{3}, n_samples::Integer)
    # calls single_point_forces a bunch
end

function single_point_forces(sys::System{3})
    # for s-TDEP
end

# Generates configurations to seed MD simulations
function generate_initial_configurations(sys::System{3}, temperature, n_seeds::Integer)
    # Run canonical_configurations to sample phase space
    cc = CanonicalConfiguration(temperature = temperature, nconf = n_seeds)
end


end