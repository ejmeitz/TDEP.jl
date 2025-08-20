export write_meta, write_partial_stat, write_ssposcar

function write_meta(outdir::String, temperature, n_samples, dt_fs, n_atoms)

    outpath = joinpath(outdir, "infile.meta")

    if isfile(outpath)
        @warn "Overwriting existing infile.meta at $(outpath)"
    end

    open(outpath, "w") do f
        println(f, "$(n_atoms) # Number of atoms")
        println(f, "$(n_samples) # Number of samples")
        println(f, "$(ustrip(dt_fs)) # Timestep in fs")
        println(f, "$(ustrip(temperature)) # Temperature in K")
    end
    
end

# Makes an infile.stat with the data that
# is easy to get access to
function write_partial_stat(outdir::String, PE::AbstractVector,
                             KE::AbstractVector, T::AbstractVector;
                             sampled_every::Integer = 1, dt_fs = 1.0,
                             file_mode = "w")

    outpath = joinpath(outdir, "infile.stat")

    if isfile(outpath) && file_mode == "w"
        @warn "Overwriting existing infile.stat at $(outpath)"
    end

    N = length(PE)

    @assert N == length(KE) "Cannot make infile.stat from vectors of different length"
    @assert N == length(T) "Cannot make infile.stat from vectors of different length"

    TE = ustrip.(PE) .+ ustrip.(KE)
    steps = Int.(sampled_every .* collect(1:N))
    timesteps = Int.(steps .* ustrip(dt_fs))

    tmp = rand(7)

    open(outpath, file_mode) do f
        for i in 1:N
            @printf f "%d %d %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n" steps[i] timesteps[i] ustrip(TE[i]) ustrip(PE[i]) ustrip(KE[i]) ustrip(T[i]) tmp...
        end
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

function read_poscar_symbol_block(filepath::String)
    species_line = ""
    count_line = ""
    open(filepath, "r") do f
        for _ in 1:5 
            readline(f)
        end
        species_line = readline(f)
        count_line = readline(f)
    end
    symbols = Symbol.(split(strip(species_line)))
    counts = parse.(Int, split(strip(count_line)))

    return symbols, counts
end

function to_frac_coords(cell::AbstractMatrix{L}, position::AbstractVector{L}) where L
    return mod.(cell \ position, L(1.0))
end

function write_ssposcar(outdir::String, cell_vectors::AbstractMatrix{L},
                         positions::AbstractVector{<:AbstractVector{L}}, 
                         atomic_symbols) where L

    @assert length(positions) == length(atomic_symbols) "Cannot write ssposcar with mismatched lengths"

    # Molly hard-codes unknown if can't parse symbol
    @assert isnothing(findfirst(x -> x == :unknown, atomic_symbols)) "Cannot make ssposcar with unknown symbols."

    outpath = joinpath(outdir, "infile.ssposcar")
    
    if isfile(outpath)
        @warn "Overwriting existing infile.ssposcar at $(outpath)"
    end
    
    # If there are units, convert to Angstroms
    # and then strip the units
    if unit(L) != NoUnits
        cell_vectors = ustrip.(u"Å", cell_vectors)
        positions = [ustrip.(u"Å", p) for p in positions]
    else
        @warn "Got no length units on input. Assuming Angstroms."
    end

    frac = Vector{Vector{Float64}}(undef, length(positions))
    for i in eachindex(positions)
        frac[i] = to_frac_coords(cell_vectors, positions[i])
    end

    symbol_line, count_line = poscar_symbol_block(atomic_symbols)

    open(outpath, "w") do f
        # Write header
        println(f, "SSPOSCAR")
        @printf f "%.10f\n" 1.0
        @printf f "%.10f %.10f %.10f\n" cell_vectors[1,:]...
        @printf f "%.10f %.10f %.10f\n" cell_vectors[2,:]...
        @printf f "%.10f %.10f %.10f\n" cell_vectors[3,:]...
        println(f, symbol_line)
        println(f, count_line)
        println(f, "Direct coordinates")
        for i in eachindex(frac)
            @printf f "%.15f %.15f %.15f\n" frac[i]...
        end
    end
end

function read_poscar_cell(path; n_atoms = nothing)

    cell = zeros(Float64, 3, 3)

    open(path, "r") do f
        readline(f)
        scale = parse(Float64, readline(f))
        lv1 = scale .* parse.(Float64, split(strip(readline(f))))
        lv2 = scale .* parse.(Float64, split(strip(readline(f))))
        lv3 = scale .* parse.(Float64, split(strip(readline(f))))

        cell .= hcat(lv1, lv2, lv3) # cell vecs as columns

        readline(f) # skip species line

        natoms_file = sum(parse.(Int, split(strip(readline(f)))))
        if !isnothing(n_atoms) && natoms_file != n_atoms
            error(ArgumentError("Poscar has $(natoms_file) but you told me it would have $(natoms)"))
        end
        n_atoms = natoms_file

    end

    return cell, n_atoms

end

# Just reads the positions and cell from POSCAR
function read_poscar_positions(path; n_atoms = nothing,
                                ssposcar_is_frac::Bool = true,
                                store_frac_coords::Bool = false)

                                
    cell, n_atoms = read_poscar_cell(path; n_atoms = n_atoms)

    positions = zeros(SVector{3, Float64}, n_atoms)

    return read_poscar_positions!(positions, path;
                                    n_atoms = n_atoms, 
                                    ssposcar_is_frac = ssposcar_is_frac,
                                    store_frac_coords = store_frac_coords)

end

function read_poscar_positions!(positions::Vector{SVector{3,T}}, path;
                                n_atoms = nothing, 
                                ssposcar_is_frac::Bool = true,
                                store_frac_coords::Bool = false) where T

    cell, n_atoms = read_poscar_cell(path; n_atoms = n_atoms)

    convert_to_cart = (!store_frac_coords && ssposcar_is_frac)

    if convert_to_cart
        parse_line = (line) -> SVector(cell * parse.(T, split(strip(line))[1:3])...)
    else
        parse_line = (line) -> SVector(parse.(T, split(strip(line))[1:3])...)
    end

    open(path, "r") do f
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f) # skip species line
        readline(f) # natoms line
        readline(f) # skip "direct coordinates" line

        for i in 1:n_atoms
            positions[i] = parse_line(readline(f))
        end
    end

    return positions, cell

end

function read_poscar_positions!(positions::AbstractMatrix{T}, path;
                                n_atoms = nothing, 
                                ssposcar_is_frac::Bool = true,
                                store_frac_coords::Bool = false) where T

    if size(positions, 1) != 3
        raise(ArgumentError("read_poscar_positions only supports 3D systems. Expected 3xN matrix."))
    end

    cell, n_atoms = read_poscar_cell(path; n_atoms = n_atoms)

    convert_to_cart = (!store_frac_coords && ssposcar_is_frac)

    if convert_to_cart
        parse_line = (line) -> cell * parse.(T, split(strip(line))[1:3])
    else
        parse_line = (line) -> parse.(T, split(strip(line))[1:3])
    end

    open(path, "r") do f
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f) # skip species liner
        readline(f) # atoms line
        readline(f) # skip "direct coordinates" line

        for i in 1:n_atoms
            positions[:,i] = parse_line(readline(f))
        end
    end

    return positions, cell

end