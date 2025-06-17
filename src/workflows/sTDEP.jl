
function TDEP.sTDEP()


    # Generate initial configurations with maximum freuqency

    # Generate remaining configurations with IFCs from prior iteration

end

function sTDEP_iter(sys::System{3}, cc::CanoncaiConfiguration)

end


function generate_sTDEP_dataset(sys::System{3}, nconf, temperature, quantum; max_freq = nothing)

    if isnothing(max_freq)
        cc = CanonicalConfiguration(temperature = temperature, nconf = nconf, quantum = quantum)
    else
        cc = CanonicalConfiguration(temperature = temperature, nconf = nconf,
                                         quantum = quantum, maximum_frequency = max_freq)
    end

    return sTDEP_iter(sys, cc)
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

function TDEP.single_point_forces(sys::System{3})
    return Molly.forces(sys)
end
