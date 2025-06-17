export sTDEP

"""
    sTDEP(
        sys::AbstractSystem{3},
        calc,
        basedir::String,
        niter::Integer,
        nconf::Union{Integer, AbstractVector{Integer}},
        rc2,
        temperature,
        maximum_frequency;
        quantum::Bool = false,
        ncores::Integer = Threads.nthreads(),
        verbose::Bool = false
    )

Calculates force constants self-consistently by sampling configurations from the canonical ensemble,
fitting force constants and iterating until convergence. 

This function is written to take in any `AbstractSystem` from the AtomsBase interface and any force calculator
from the AtomsCalculators interface. This means you can create systems using ASE, SimpleCrystals etc. and
use Molly.jl, DFTK.jl, ASEconvert.jl etc. to calculate forces.

# Arguments:
- `sys` : AbstractSystem containing at a minimum the positions and cell definition.
- `calc` : Force calculator (e.g. DFTKCalculator, MollyCalculator, ASECalculator)
- `basedir::String` : Directory where all output files are written
- `niter::Integer` : Number of self-consistent iterations to perform. Includes iter to initialize from freuency.
- `nconf` : If Integer, the number of configurations to generate. If list of integers, the number
     of configuritons to generate at each iteration.
- `rc2` : Force constant cutoff in Angstroms when calculating second-order IFCs from configurations
- `temperature` : Temperature to generate configurations at
- `maximum_frequency` : Estimate of maximum freuqency in THz to generate initial configurations from
- `quantum::Bool = false` : Generate configurations from the quantum canonical ensemble
- `ncores::Integer = Threads.nthreads()` : Number of cores used by IFC calculation
- `verbose::Bool = false` : Enable extra printing
"""
function sTDEP(
        sys::AbstractSystem{3},
        calc,
        basedir::String,
        niter::Integer,
        nconf::Union{Integer, AbstractVector{Integer}},
        rc2,
        temperature,
        maximum_frequency;
        quantum::Bool = false,
        ncores::Integer = Threads.nthreads(),
        verbose::Bool = false
    )
    #* Add pre-conditioning ?

    @assert isfile(joinpath(basedir, "infile.ucposcar")) "sTDEP requires an infile.ucposcar in $(basedir)"
    @assert AtomsCalculators.length_unit(calc) == u"Å" "Expected angstroms as length unit got $(AtomsCalculators.length_unit(calc))"
    @assert AtomsCalculators.energy_unit(calc) == u"eV" "Expected eV as energy unit got $(AtomsCalculators.energy_unit(calc))"
    @assert all(periodicity(cell(sys))) "Must be periodic system"

    if typeof(nconf) <: AbstractVector
        @assert length(nconf) == niter "Length of config list, $(length(nconf)) does not match niter: $(niter)"
    else
        nconf = [nconf for _ in 1:niter]
    end

    # Make ssposcar
    cv = hcat(AtomsBase.cell_vectors(sys)...)
    write_ssposcar(basedir, cv, AtomsBase.position(sys, :), Symbol.(AtomsBase.atomic_symbol(sys, :)))

    get_path = (i) -> joinpath(basedir, "iter$(lpad(i,3,'0'))")

    efc = ExtractForceConstants(secondorder_cutoff = rc2)
    pd = PhononDispersionRelations(dos = true)

    init_dir = get_path(0)
    mkdir(init_dir)

    # Generate initial configurations with maximum freuqency
    cp(joinpath(basedir, "infile.ssposcar"), joinpath(init_dir, "infile.ssposcar"))
    cp(joinpath(basedir, "infile.ucposcar"), joinpath(init_dir, "infile.ucposcar"))

    cc_init = CanonicalConfiguration(
        temperature = temperature,
        nconf = nconf[1],
        quantum = quantum, 
        maximum_frequency = maximum_frequency
    )
    generate_configs(sys, cc_init, calc, init_dir, verbose)

    # Generate remaining configurations with IFCs from prior iteration
    p = Progress(niter - 1)
    for i in 1:(niter - 1)

        cc = CanonicalConfiguration(
            temperature = temperature,
            nconf = nconf[i+1],
            quantum = quantum, 
        )

        outdir = get_path(i)
        # Move IFCs from last iter to current dir
        prepare_next_dir(get_path(i-1), outdir, i == 1)
        # Generate Configs Given Current IFCs
        generate_configs(sys, cc, calc, outdir, verbose)
        # Calculate IFCs to Generate Next Set of Configs
        execute(efc, outdir, ncores, verbose)
        # Generate DOS and Dispersion Data
        execute(pd, outdir, ncores, verbose)
        next!(p)
    end
    finish!(p)

end

function prepare_next_dir(current_dir, dest_dir, init_pass::Bool = false)
    mkdir(dest_dir)
    cp(joinpath(current_dir, "infile.ssposcar"), joinpath(dest_dir, "infile.ssposcar"))
    cp(joinpath(current_dir, "infile.ucposcar"), joinpath(dest_dir, "infile.ucposcar"))

    if init_pass
        cp(joinpath(current_dir, "outfile.fakeforceconstant"), joinpath(dest_dir, "infile.forceconstant"))
    else
        cp(joinpath(current_dir, "outfile.forceconstant"), joinpath(dest_dir, "infile.forceconstant"))
    end
end

function generate_configs(sys::AbstractSystem{3}, cc::CanonicalConfiguration, 
                            calc, outdir::String, verbose::Bool)

    execute(cc, outdir, 1, verbose)

    get_filepath = (i) -> joinpath(outdir, "contcar_conf$(lpad(i, 4, '0'))")

    E_units = (AtomsCalculators.energy_unit(calc) == NoUnits) ? NoUnits : u"eV"
    F_units = (E_units == NoUnits) ? NoUnits : u"eV / Å"

    # Parse coordinates into sys object and calculate forces
    for i in 1:cc.nconf
        filepath = get_filepath(i)
        x_frac, cell_vec = read_poscar_positions(filepath, length(sys))
        cell_vec = cell_vec * u"Å"

        x_cart = [SVector((cell_vec*xf)...) for xf in x_frac]
       
        # Construct AtomsBase.FastSystem object with new positions
        #* CURRENTLY WILL BREAK FOR LJ w/ Molly + MISSING VELOCITIES
        fs = FastSystem(
                AtomsBase.cell(sys),
                x_cart,
                AtomsBase.ChemicalSpecies.(AtomsBase.atomic_symbol(sys, :)), 
                AtomsBase.mass(sys, :)
            )

        F = AtomsCalculators.forces(fs, calc)
        println(F[1])
        PE = uconvert(E_units, AtomsCalculators.potential_energy(fs, calc))

        # Add data to infile.forces
        open(joinpath(outdir, "infile.forces"), "a") do ff
            for j in eachindex(F)
                @printf ff "%.15f %.15f %.15f\n" ustrip.(F_units, F[j])...
            end
        end

        # Add data to infile.positions
        open(joinpath(outdir, "infile.positions"), "a") do pf
            for j in eachindex(x_frac)
                @printf pf "%.15f %.15f %.15f\n" ustrip.(x_frac[j])...
            end
        end

        # Add PE to infile.stat
        write_partial_stat(
            outdir, 
            [PE],
            [0.0],
            [cc.temperature];
            file_mode = "a"
        )
    end

    # Write infile.meta
    write_meta(outdir, ustrip.(cc.temperature), cc.nconf, 1.0, length(sys))

end
