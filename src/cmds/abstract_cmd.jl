export execute


"""
Base type for TDEP commands. All classes which inhert from `TDEP_Command`
implement `cmd_name` and other helper functions necessary to check 
required input files and parse output files.
"""
abstract type TDEP_Command{T<:Real} end

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
