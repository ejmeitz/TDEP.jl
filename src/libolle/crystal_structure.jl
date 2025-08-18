# Simplified wrapper for TDEP crystal structure using opaque pointers


struct CrystalStructure
    handle::Ptr{Cvoid}
end

function CrystalStructure(filename::String)
    empty_handle = Ptr{Cvoid}()
    return @ccall "TDEP".__type_crystalstructure_MOD_readfromfile(
        filename::Cstring
    )::CrystalStructure
end


# function CrystalStructure(lattice_vectors::AbstractMatrix,
#                           x_frac::AbstractMatrix,
#                           atomic_numbers::Vector{<:Integer},
#                           length_unit::Symbol = :Angstrom) where T

#     L = 0
#     if length_unit == :Angstrom
#         L = 0
#     elseif length_unit == :Bohr
#         L = 1
#     else
#         error("Invalid length unit: $length_unit. Expected :Angstrom or :Bohr")
#     end

#     empty_handle = Ptr{CVoid}()
#     return CrystalStructure(
#         @ccall "TDEP".create_crystalstructure(
#             lattice_vectors::Ptr{Cvoid},
#             x_frac::Ptr{Cvoid},
#             atomic_numbers::Ptr{Cint},
#             L::Ref{Int32}
#         )::CrystalStructure
#     )

# end