export execute


"""
Base type for TDEP commands. All classes which inhert from `TDEP_Command`
implement `cmd_name` and other helper functions necessary to check 
required input files and parse output files.
"""
abstract type TDEP_Command{T} end

# BinaryBuilder creates Julia functions for each 
# ExecutableProduct. This function returns a
# handle to that function so that I can invoke it.
# Relevant Docs: https://docs.binarybuilder.org/stable/jll/#ExecutableProduct
function str_to_fn(str::TDEP_Command)
    return getfield(TDEP_jll, Symbol(cmd_name(str)))
end

#* UPDATE THIS TO HANDLE TUPLES 
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

function build_cmd(tdep_cmd::Cmd, tdep_args)
    return `$(tdep_cmd) $(tdep_args)`
end

function build_mpi_cmd(tdep_cmd::Cmd, tdep_args, ncores::Integer)

    # Combine environments from the two 
    # JLL command wrappers. Julia can't handle
    # two wrappers in the same Cmd automatically
    if !isnothing(tdep_cmd.env)
        mpi_exec_cmd = addenv(MPI.mpiexec(), tdep_cmd.env)
    else
        mpi_exec_cmd = MPI.mpiexec()
    end

    return `$(mpi_exec_cmd) -n $(ncores) $(tdep_cmd.exec) $(tdep_args)`
end

function execute(cmd::TDEP_Command, rundir::String = pwd(), ncores = Threads.nthreads(),
                    verbose::Bool = false)

    cmd_str = cmd_name(cmd)
    args = build_args(cmd)
    log_file = cmd_str * ".log"
    f = str_to_fn(cmd)

    if verbose
        args = [args; "--verbose"]
    end
    
    tdep_cmd = f()
    if ncores > 1
        tdep_cmd = build_mpi_cmd(tdep_cmd, args, ncores)
    elseif ncores == 1
        tdep_cmd = build_cmd(tdep_cmd, args)
    else
        @error "Got negative number of cores. Huh???"
    end

    cd(rundir) do
        verbose && @info "Running: $(join(tdep_cmd.exec, ' '))"
        run(pipeline(tdep_cmd; stdout = log_file))
    end
    
    return nothing

end
