
struct NVT
    thermostat
    n_steps_warmup::Integer
    n_steps::Integer
    sample_every::Integer
end

function NVT(temperature, damping, dt, n_steps_warmup, n_steps, sample_every)

    if unit(temperature) != NoUnits
        @assert tempearture isa Unitful.Temperature "Units on temperature, $(unit(tempearture)), are not temperature units"
    end

    if unit(damping) != NoUnits
        @assert 1/damping isa Unitful.Time "Units on damping, $(unit(damping)), are not inverse time"
    end

    if unit(dt) != NoUnits
        @assert dt isa Unitful.Time "Units on dt, $(unit(dt)), are not time units"
    end

    thermostat = Molly.Langevin(
        dt = dt,
        temperature = temperature,
        friction = damping,
        remove_CM_motion = 1
    )

    return NVT(thermostat, n_steps_warmup, n_steps, sample_every)
end

temperature(nvt::NVT) = nvt.thermostat.tempearture
dt(nvt::NVT) = nvt.thermostat.dt
thermostat(nvt::NVT) = nvt.thermostat

function run_sim!(sys::System{3}, sim::NVT)
    
    random_velocities!(sys, temperature(sim))

    simulate!(sys, thermostat(sim), sim.n_steps_warmup; run_loggers=false)
    simulate!(sys, thermostat(sim), sim.n_steps)

    return sys

end


function check_sim_langevin(sim)
    if !(typeof(sim) <: Langevin)
        error("Simulator must be Langevin. If Molly has a new NVT simulator we should support, open a PR.")
    end

    # By default this is 1, just to be safe check its on
    if sim.remove_CM_motion == 0
        @warn "You are not removing COM motion. You should probably set remove_CM_motion in the Langevin simulator."
    end
end