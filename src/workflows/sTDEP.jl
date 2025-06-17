export sTDEP

#* Add pre-conditioning ?
#* Generate ucposcar and ssposcar?
function sTDEP(
        sys::AbstractSystem{3},
        calc,
        basedir::String,
        niter::Integer,
        nconf::Integer,
        rc2,
        temperature,
        quantum::Bool,
        maximum_frequency;
        ncores::Integer = Threads.nthreads(),
        verbose = true
    )

    get_path = (i) -> joinpath(basedir, "iter$(lpad(i,3,'0'))")

    init_dir = get_path(0)
    mkdir(init_dir)

    efc = ExtractForceConstants(secondorder_cutoff = rc2)
    pd = PhononDispersionRelations(dos = true)

    # Generate initial configurations with maximum freuqency
    cc_init = CanonicalConfiguration(
                temperature = temperature,
                nconf = nconf,
                quantum = quantum, 
                maximum_frequency = maximum_frequency
            )

    generate_configs(sys, cc_init, calc, init_dir, ncores, verbose)

    # Generate remaining configurations with IFCs from prior iteration
    cc = CanonicalConfiguration(
            temperature = temperature,
            nconf = nconf,
            quantum = quantum, 
        )
    p = ProgressMeter(niter)
    for i in 1:niter
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
        cp(joinpath(current_dir, "outfile.forceconstant"), joinpath(dest_dir, "infile.forceconstant"))
    else
        cp(joinpath(current_dir, "outfile.fakeforceconstant"), joinpath(dest_dir, "infile.forceconstant"))
    end
end


function generate_configs(sys::AbstractSystem{3}, cc::CanoncaiConfiguration, 
                    calc, outdir::String, verbose::Bool)

    execute(cc, outdir, 1, verbose)

    # Parse coordinates into sys object and calculate forces

    # Generate necessary input files for extract_forceconstants
end
