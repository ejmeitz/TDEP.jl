export NVT


function TDEP.NVT(temp, damping, dt, n_steps_warmup, n_steps, sample_every)

    if unit(temp) != NoUnits
        @assert temp isa Unitful.Temperature "Units on temperature, $(unit(temp)), are not temperature units"
    end

    if unit(damping) != NoUnits
        @assert 1/damping isa Unitful.Time "Units on damping, $(unit(damping)), are not inverse time"
    end

    if unit(dt) != NoUnits
        @assert dt isa Unitful.Time "Units on dt, $(unit(dt)), are not time units"
    end

    thermostat = Molly.Langevin(
        dt = dt,
        temperature = temp,
        friction = damping,
        remove_CM_motion = 1
    )

    return NVT(thermostat, n_steps_warmup, n_steps, sample_every)
end

temperature(nvt::TDEP.NVT) = nvt.thermostat.temperature
dt(nvt::TDEP.NVT) = nvt.thermostat.dt
thermostat(nvt::TDEP.NVT) = nvt.thermostat

function run_sim!(sys::Molly.System{3}, sim::TDEP.NVT; run_loggers_warmup = false)
    
    random_velocities!(sys, temperature(sim))

    simulate!(sys, thermostat(sim), sim.n_steps_warmup; run_loggers=run_loggers_warmup)
    simulate!(sys, thermostat(sim), sim.n_steps)

    return sys

end


# function check_sim_langevin(sim)
#     if !(typeof(sim) <: Langevin)
#         error("Simulator must be Langevin. If Molly has a new NVT simulator we should support, open a PR.")
#     end

#     # By default this is 1, just to be safe check its on
#     if sim.remove_CM_motion == 0
#         @warn "You are not removing COM motion. You should probably set remove_CM_motion in the Langevin simulator."
#     end
# end