export write_meta

function write_meta(outdir::String, temperature, n_samples, dt_fs, n_atoms)

    outpath = joinpath(outdir, "infile.meta")

    if isfile(outpath)
        @warn "Overwriting existing infile.meta at $(outpath)"
    end

    open(outpath, "w") do f
        println(f, "$(n_atoms) # Number of atoms")
        println(f, "$(n_samples) # Number of samples")
        println(f, "$(dt_fs) # Timestep in fs")
        println(f, "$(temperature) # Temperature in K")
    end
    
end

# Makes an infile.stat with the data that
# is easy to get access to
function write_partial_stat(outdir::String, PE::AbstractVector,
                             KE::AbstractVector, T::AbstractVector;
                             sampled_every::Integer = 1, dt_fs = 1.0,
                             file_mode = "w")

    outpath = joinpath(outdir, "infile.stat")

    if isfile(outpath)
        @warn "Overwriting existing infile.stat at $(outpath)"
    end

    N = length(PE)

    @assert N == length(KE) "Cannot make infile.stat from vectors of different length"
    @assert N == length(T) "Cannot make infile.stat from vectors of different length"

    TE = PE .+ KE
    steps = sampled_every .* collect(1:N)
    timesteps = steps .* dt_fs
    tmp = rand(N)

    open(outpath, file_mode) do f
        writedlm(f, [steps timesteps TE PE KE T tmp tmp tmp tmp tmp tmp tmp], " ")
    end

end

function poscar_symbol_block(symbols::AbstractVector{T}) where T
    counts = Dict{T,Int}()
    order  = T[]
    last_symbol = Symbol()
    for s in symbols
        counts[s] = get(counts, s, 0) + 1
        if last_symbol != s
            push!(order, s)
            last_symbol = s
        end
    end
    symbol_line = join(order, " ")
    count_line  = join([string(counts[s]) for s in order], " ")
    return symbol_line, count_line
end


function write_ssposcar(outdir::String, cell_vectors::AbstractMatrix{L},
                         positions::AbstractVector{<:AbstractVector{L}}, 
                         atomic_symbols) where L

    @assert length(positions) == length(atomic_symbols) "Cannot write ssposcar with mismatched lengths"

    outpath = joinpath(outdir, "infile.ssposcar")
    
    if isfile(outpath)
        @warn "Overwriting existing infile.ssposcar at $(outpath)"
    end
    
    # If there are units, convert to Angstroms
    # and then strip the units
    if units(L) != NoUnits
        cell_vectors = ustrip.(u"Å", cell_vectors)
        positions = [ustrip.(u"Å", p) for p in positions]
    else
        @warn "Got no length units on input. Assuming Angstroms."
    end

    frac = Vector{Vector{Float64}}(undef, length(positions))
    for i in eachindex(positions)
        frac[i] = mod.(cell \ positions[i], 1.0)
    end

    symbol_line, count_line = poscar_symbol_block(atomic_symbols)

    open(outpath, "w") do f
        # Write header
        println(f, "SSPOSCAR")
        @printf f "%.10f" 1.0
        @printf f "%.10f %.10f %.10f" cell_vectors[1,:]...
        @printf f "%.10f %.10f %.10f" cell_vectors[2,:]...
        @printf f "%.10f %.10f %.10f" cell_vectors[3,:]...
        println(f, symbol_line)
        println(f, count_line)
        println(f, "Direct coordinates")
        for i in eachindex(frac)
            @printf f "%.15f %.15f %.15f" frac[i]...
        end
    end
end