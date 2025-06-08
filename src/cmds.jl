export execute, ExtractForceConstants, CanonicalConfiguration

abstract type TDEP_Command{T<:Real} end

Base.@kwdef struct ExtractForceConstants{T} <: TDEP_Command{T}
    secondorder_cutoff::T
    thirdorder_cutoff::T = -1.0
    fourthorder_cutoff::T = -1.0
    polar::Bool = false
    stride::Int = 1
    firstorder::Bool = false
    temperature::T = -1.0
    norotational::Bool = false
    nohuang::Bool = false
    nohermitian::Bool = false
    potential_energy_differences::Bool = true
end

"""
    CanonicalConfiguration(temperature, nconf, quantum = false, outputformat = 1,
                            mindist = -1, debye_temperature = -1.0, maximum_frequency = -1.0)

# Arguments 
- `temperature::T`: Temperature to emulate.
- `nconf::Int`: Number of configurations to generate.
- `quantum::Bool = false` : Use Bose-Einstein statistics instead of Maxwell-Boltzmann.
- `output_format::Int = 1`: (1) VASP (2) Abinit (4) FHI-Aims (5) Siesta
- `mindist::T = -1`: Smallest distance between two atoms allowed, in units of the nearest neighbour distance.
- `debye_temperature::T = -1.0`: Generate forceconstants that match a Debye temperature,
                                 and build displacements according to these. 
- `maximum_frequency::T = -1.0`: Generate forceconstants that match a maximum frequency (in THz),
                                 and build displacements according to these.

# Description
This code takes a second order forceconstant and generates a supercell
with the atoms displaced corresponding to a harmonic canonical ensemble.
These configurations can be used directly, to generate force-displacement
statistics to determine force constants and other thermal properties, 
or as the starting point for molecular dynamics.

TDEP Documentatoin: https://tdep-developers.github.io/tdep/program/canonical_configuration/
"""
Base.@kwdef struct CanonicalConfiguration{T} <: TDEP_Command{T}
    temperature::T
    nconf::Int
    quantum::Bool = false
    output_format::Int = 1
    mindist::T = -1.0
    debye_temperature::T = -1.0
    maximum_frequency::T = -1.0
end

cmd_name(::ExtractForceConstants) = "extract_forceconstants"
cmd_name(::CanonicalConfiguration) = "canonical_configuration"

# BinaryBuilder creates Julia functions for each 
# ExecutableProduct. This function returns a
# handle to that function so that I can invoke it.
# Relevant Docs: https://docs.binarybuilder.org/stable/jll/#ExecutableProduct
function str_to_fn(str::TDEP_Command)
    return getfield(TDEP_jll, Symbol(cmd_name(str)))
end

function build_args(cmd::TDEP_Command)
    T = typeof(cmd)
    args = []
    for (arg, arg_type) in zip(fieldnames(T), T.types)
        if arg_type != Bool # if not a flag
            push!(args, "--$(arg)", "$(getproperty(cmd, arg))")
	elseif getproperty(cmd, arg) == true # only push flag if set to true
	    push!(args, "--$(arg)")
        end
    end

    return args
end

#* TODO ADD MPI SUPPORT
function execute(cmd::TDEP_Command, rundir::String = pwd())

    cmd_str = cmd_name(cmd)
    args = build_args(cmd)
    log_file = cmd_str * ".log"
    f = str_to_fn(cmd)

    cd(rundir) do
	shell_cmd = `$(f()) $(args)`
        run(pipeline(shell_cmd; stdout = log_file))
    end
    
    return nothing

end
