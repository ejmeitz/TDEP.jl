# lo_ptr_api.jl
const liblo = "./liblo_ptr.so"

struct LO_Mem
    p::Ptr{Cvoid}
end

struct LO_Crystal
    p::Ptr{Cvoid}
end

struct LO_FC2
    p::Ptr{Cvoid}
end

# --- Mem ---
function LO_Mem()
    out = Ref{Ptr{Cvoid}}(); st = Ref{Cint}(0)
    ccall((:lo_mem_new_ptr, liblo), Cvoid, (Ref{Ptr{Cvoid}}, Ref{Cint}), out, st)
    st[] == 0 || error("lo_mem_new_ptr failed: ", st[])
    obj = LO_Mem(out[])
    finalizer(obj) do x
        ccall((:lo_mem_destroy_ptr, liblo), Cvoid, (Ptr{Cvoid},), x.p)
    end
    return obj
end

# --- Crystal ---
function LO_Crystal(path::AbstractString, mem::LO_Mem; verbosity::Integer=1)
    out = Ref{Ptr{Cvoid}}(); st = Ref{Cint}(0)
    ccall((:lo_crystal_from_file_ptr, liblo), Cvoid,
          (Cstring, Ptr{Cvoid}, Cint, Ref{Ptr{Cvoid}}, Ref{Cint}),
          path, mem.p, Cint(verbosity), out, st)
    st[] == 0 || error("lo_crystal_from_file_ptr failed: ", st[])
    obj = LO_Crystal(out[])
    finalizer(obj) do x
        ccall((:lo_crystal_destroy_ptr, liblo), Cvoid, (Ptr{Cvoid},), x.p)
    end
    return obj
end

# --- FC2 ---
function LO_FC2(path::AbstractString, crys::LO_Crystal, mem::LO_Mem; verbosity::Integer=1)
    out = Ref{Ptr{Cvoid}}(); st = Ref{Cint}(0)
    ccall((:lo_fc2_from_file_ptr, liblo), Cvoid,
          (Cstring, Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ref{Ptr{Cvoid}}, Ref{Cint}),
          path, crys.p, mem.p, Cint(verbosity), out, st)
    st[] == 0 || error("lo_fc2_from_file_ptr failed: ", st[])
    obj = LO_FC2(out[])
    finalizer(obj) do x
        ccall((:lo_fc2_destroy_ptr, liblo), Cvoid, (Ptr{Cvoid},), x.p)
    end
    return obj
end

# (Optional) example method
function potential_energy(fc::LO_FC2, u::Vector{Float64})
    e  = Ref{Cdouble}(0)
    st = Ref{Cint}(0)
    ccall((:lo_fc2_potential_energy_ptr, liblo), Cvoid,
          (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Ref{Cdouble}, Ref{Cint}),
          fc.p, u, Cint(length(u)), e, st)
    st[] == 0 || error("lo_fc2_potential_energy_ptr failed: ", st[])
    return e[]
end
