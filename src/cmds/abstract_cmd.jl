export execute


"""
Base type for TDEP commands. All classes which inhert from `TDEP_Command`
implement `cmd_name` and other helper functions necessary to check 
required input files and parse output files.
"""
abstract type TDEP_Command{T} end

# BinaryBuilder creates Julia functions for each 
# ExecutableProduct. This function gets that handle and
# invokes it to return the resulting `Cmd` object.
# Relevant Docs: https://docs.binarybuilder.org/stable/jll/#ExecutableProduct
function to_julia_cmd(str::TDEP_Command)
    handle = getfield(TDEP_jll, Symbol(cmd_name(str)))
    cmd = handle()
    return cmd
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

@doc raw"""
When using an `ExecutableProduct` from a JLL the standard way to get the right environment is
to invoke the automatically generated wrapper. See the docs [here](https://docs.binarybuilder.org/stable/jll/#ExecutableProduct).

For example,
```julia
cmd = `$(mpiexec()) -n 4 some_program`
```

If you want to use two of these wrappers in the same `Cmd` object, Julia will object. For example, this does not work.
```julia
cmd = `$(mpiexec()) -n 4 $(another_program())`
```

To get around this, we can manually merge the environments of each `Cmd`. When the `ExecutableProduct`
wrapper is invoked it returns a `Cmd` object with the proper environment variables set to locate 
the artifacts that Julia downloaded (e.g. where mpiexec lives). The `addenv` function combines those
environments together and we can re-build the command with the correct paths. This should be relatively
safe as long as MPI is a dependency of the library which generated the second `ExecutableProduct`. 
"""
function build_mpi_cmd(tdep_cmd::Cmd, tdep_args, ncores::Integer)

    if isnothing(tdep_cmd.env)
        mpi_exec_cmd = MPI.mpiexec()
    else
        mpi_exec_cmd = addenv(MPI.mpiexec(), tdep_cmd.env)
    end

    return `$(mpi_exec_cmd) -n $(ncores) $(tdep_cmd.exec) $(tdep_args)`
end

function execute(cmd::TDEP_Command, rundir::String = pwd(), ncores = Threads.nthreads(),
                    verbose::Bool = false)

    cmd_str = cmd_name(cmd)
    args = build_args(cmd)
    log_file = cmd_str * ".log"
    tdep_cmd = to_julia_cmd(cmd)

    if verbose
        args = [args; "--verbose"]
    end
    
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
