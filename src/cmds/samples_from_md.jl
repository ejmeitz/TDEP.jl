export SamplesFromMD

"""
    SamplesFromMD(; nsamples = 50, output_format = 1)

# Arguments
- `nsamples::Int = 50`: Number of representative samples to select.
- `output_format::Int = 1`: Output format.  
  - `1` = VASP  
  - `2` = Abinit  
  - `4` = FHI-AIMS  
  - `5` = Siesta

# Description
Selects approximately uncorrelated, evenly spaced samples from a molecular dynamics simulation that reproduce the mean and standard deviation of both potential and kinetic energies.

TDEP Documentation: https://tdep-developers.github.io/tdep/program/samples_from_md/
"""
Base.@kwdef struct SamplesFromMD <: TDEP_Command{Nothing}
    nsamples::Int = 50
    output_format::Int = 1
end

cmd_name(::SamplesFromMD) = "samples_from_md"