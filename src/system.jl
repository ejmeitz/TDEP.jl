export TDEPSystem

"""
    TDEPSystem(ssposcar_path::String; 
               ssposcar_is_frac::Bool = true,
               store_frac_coords::Bool = false)

Constructs an AtomsBase compatible system object given the path to
an ssposcar file. `ssposcar_is_frac` sets if the poscar file is in 
fractional coordinates or not. `store_frac_coords` dicates where
this object stores the coordinates as fractional or Cartesian. 
"""
struct TDEPSystem{TCELL, L <: Unitful.Length, M <: Unitful.Mass, S} <: AbstractSystem{3}
    cell::TCELL
    position::Vector{SVector{3, L}}
    species::Vector{S}
    mass::Vector{M}
    print_str::String
end


function TDEPSystem(ssposcar_path::String; 
                    ssposcar_is_frac::Bool = true,
                    store_frac_coords::Bool = false)

    posns, cell_vec = read_poscar_positions(ssposcar_path; 
                                            store_frac_coords = store_frac_coords,
                                            ssposcar_is_frac = ssposcar_is_frac)
    cell_vec = cell_vec * u"Å"
    posns = posns * u"Å"

    symbols, counts = read_poscar_symbol_block(ssposcar_path)
    str = join(["$(c) $(string(s)) atoms" for (s, c) in zip(symbols, counts)])
    print_str = "TDEP system with $(str)"

    syms = reduce(vcat, [fill(s, c) for (s, c) in zip(symbols, counts)])
    m = [periodic_table[s].atomic_mass for s in syms]

    syms = AtomsBase.ChemicalSpecies.(syms)

    c = (SVector(cell_vec[1:3, 1]...),
            SVector(cell_vec[1:3, 2]...),
            SVector(cell_vec[1:3, 3]...)
        )
    cell = PeriodicCell(; cell_vectors = c, periodicity = (true, true, true))
    
    L = eltype(first(posns))
    return TDEPSystem{typeof(cell), L, eltype(m), eltype(syms)}(
            cell,
            posns,
            syms, 
            m,
            print_str
        )

end

# Just trust the user knows what they're doing
function TDEPSystem(sys::TDEPSystem{T,L,M,S}, new_positions::Vector{SVector{3, L}}) where {T,L,M,S}
    return TDEPSystem(
        sys.cell,
        new_positions,
        sys.species,
        sys.mass,
        sys.print_str
    )
end

Base.length(sys::TDEPSystem)         = length(sys.position)
Base.size(sys::TDEPSystem)           = size(sys.position)

Base.getindex(sys::TDEPSystem, i::Integer)  = AtomView(sys, i)

# System property access
function Base.getindex(system::TDEPSystem, x::Symbol)
    if x === :cell_vectors
        cell_vectors(system)
    elseif x === :periodicity
        periodicity(system)
    else
        throw(KeyError(x))
    end
end
Base.haskey(::TDEPSystem, x::Symbol) = x in (:cell_vectors, :periodicity)
Base.keys(::TDEPSystem) = (:cell_vectors, :periodicity)

# Atom and atom property access
AtomsBase.atomkeys(::TDEPSystem) = (:position, :species, :mass)
AtomsBase.cell(ts::TDEPSystem) = ts.cell

AtomsBase.hasatomkey(system::TDEPSystem, x::Symbol) = x in atomkeys(system)

function Base.getindex(system::TDEPSystem, i::Union{Integer,AbstractVector}, x::Symbol)
    getfield(system, x)[i]
end

Base.getindex(system::TDEPSystem, ::Colon, x::Symbol) = getfield(system, x)

AtomsBase.position(s::TDEPSystem, ::Colon) = s.position
AtomsBase.position(sys::TDEPSystem, i::Union{Integer, AbstractVector}) = sys.position[i]

AtomsBase.mass(s::TDEPSystem, ::Colon) = s.mass
AtomsBase.mass(sys::TDEPSystem, i::Union{Integer, AbstractVector}) = sys.mass[i]

AtomsBase.species(s::TDEPSystem, ::Colon) = s.species
AtomsBase.species(sys::TDEPSystem, i::Union{Integer, AbstractVector}) = sys.species[i]

AtomsBase.atomic_symbol(s::TDEPSystem, ::Colon) = Symbol.(s.species)
AtomsBase.atomic_symbol(s::TDEPSystem, i::Union{Integer, AbstractVector}) = Symbol(s.species[i])

AtomsBase.atomic_number(s::TDEPSystem, ::Colon) = atomic_number.(s.species)
AtomsBase.atomic_number(s::TDEPSystem, i::Union{Integer, AbstractVector}) = atomic_number(s.species[i])


function Base.show(io::IO, sys::TDEPSystem)
    print(io, sys.print_str)
end

function Base.show(io::IO, ::MIME"text/plain", sys::TDEPSystem)
    print(io, sys.print_str)
end