export generate_MDTDEP_dataset, NVT, TI


# This file contains function handles which are implemented
# as part of package extensions. These are separated as they
# require a heavy dependency (Molly.jl) which not all
# users may want to download.

function generate_MDTDEP_dataset end

struct NVT
    thermostat
    n_steps_warmup::Integer
    n_steps::Integer
    sample_every::Integer
end

function TI end