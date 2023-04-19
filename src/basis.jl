import Base.@property
import Base.@cfunction
import Base.getproperty
import Base.@static
using LinearAlgebra
using StaticArrays
using Combinatorics: doublefactorial
# Define the struct for Basis
mutable struct Basis
    origin::Vector{Float64}
    shell::Vector{Int64}
    num_exps::Int64
    exps::Ptr{Float64}
    coefs::Ptr{Float64}
    norm::Ptr{Float64}
end

# Define functions to allocate and deallocate Basis objects
@inline function basis_alloc(origin::Vector{Float64}, shell::Vector{Int64}, num_exps::Int64, exps::Vector{Float64}, coefs::Vector{Float64})::Ptr{Basis}
    b = ccall(:malloc, Ptr{Basis}, (UInt64,), sizeof(Basis))
    if b == C_NULL
        throw(MemoryError())
    end
    b = unsafe_pointer_to_objref(b)
    b.origin = origin
    b.shell = shell
    b.num_exps = num_exps
    b.exps = ccall(:malloc, Ptr{Float64}, (UInt64,), num_exps*sizeof(Float64))
    b.coefs = ccall(:malloc, Ptr{Float64}, (UInt64,), num_exps*sizeof(Float64))
    b.norm = ccall(:malloc, Ptr{Float64}, (UInt64,), num_exps*sizeof(Float64))
    for i in 1:3
        b.origin[i] = origin[i]
        b.shell[i] = shell[i]
    end
    for i in 1:num_exps
        b.exps[i-1] = exps[i]
        b.coefs[i-1] = coefs[i]
        b.norm[i-1] = 0.0
    end
    basis_normalize(b)
    return b
end



#instead of using ccall to call malloc and sizeof, we use the malloc and sizeof functions directly from the Base module.
# To set the fields of the Basis struct, we use unsafe_store! and unsafe_load functions to manipulate pointers.
function basis_alloc(origin::Vector{Float64}, shell::Vector{Int64}, num_exps::Int64, exps::Vector{Float64}, coefs::Vector{Float64})::Ptr{Basis}
    b = Base.unsafe_convert(Ptr{Basis}, malloc(sizeof(Basis)))
    if b == C_NULL
        throw(MemoryError())
    end
    unsafe_store!(b, Basis(origin, shell, num_exps, RefValue{Ptr{Float64}}(), RefValue{Ptr{Float64}}(), RefValue{Ptr{Float64}}()))
    b_val = unsafe_load(b)
    b_val.exps = Base.unsafe_convert(Ptr{Float64}, malloc(num_exps * sizeof(Float64)))
    b_val.coefs = Base.unsafe_convert(Ptr{Float64}, malloc(num_exps * sizeof(Float64)))
    b_val.norm = Base.unsafe_convert(Ptr{Float64}, malloc(num_exps * sizeof(Float64)))
    for i in 1:3
        b_val.origin[i] = origin[i]
        b_val.shell[i] = shell[i]
    end
    for i in 1:num_exps
        b_val.exps[i-1] = exps[i]
        b_val.coefs[i-1] = coefs[i]
        b_val.norm[i-1] = 0.0
    end
    basis_normalize(b)
    return b
end


#This implementation uses Julia's finalizer function to register a finalizer that gets called when the object is garbage collected.
# The finalizer function _basis_finalizer does the actual deallocation of the memory blocks allocated by the Basis struct. 
#The unsafe_load function is used to load the struct pointed to by p, and unsafe_free is used to deallocate the memory blocks pointed to by the pointers in the Basis struct.
function basis_dealloc(b::Ptr{Basis})
    if !isnull(b)
        finalizer(b, _basis_finalizer)
    end
end

function _basis_finalizer(p::Ptr{Basis})
    b = unsafe_load(p)
    if !isnull(b.origin)
        unsafe_free(b.origin)
    end
    if !isnull(b.shell)
        unsafe_free(b.shell)
    end
    if !isnull(b.exps)
        unsafe_free(b.exps)
    end
    if !isnull(b.coefs)
        unsafe_free(b.coefs)
    end
    if !isnull(b.norm)
        unsafe_free(b.norm)
    end
    unsafe_free(p)
end



function basis_normalise(basis::Basis)
    l, m, n = basis.shell
    L = l + m + n
    # normalize primitives first (PGBFs)
    for ia in 1:basis.num_exps
        basis.norm[ia] = sqrt(2^(2*(l+m+n)+1.5) *
                              (basis.exps[ia])^(l+m+n+1.5) /
                              (fact2(2*l-1) *
                               fact2(2*m-1) *
                               fact2(2*n-1) *
                               pi^1.5))
    end

    # now normalize the contracted basis functions (CGBFs)
    # Eq. 1.44 of Valeev integral whitepaper
    prefactor = (pi^1.5 *
                 fact2(2*l - 1) *
                 fact2(2*m - 1) *
                 fact2(2*n - 1) /
                 2.0^L)

    N = 0.0
    for ia in 1:basis.num_exps
        for ib in 1:basis.num_exps
            N += basis.norm[ia] * basis.norm[ib] *
                 basis.coefs[ia] * basis.coefs[ib] /
                 (basis.exps[ia] + basis.exps[ib])^(L+1.5)
        end
    end

    N *= prefactor
    N = N^(-0.5)
    for ia in 1:basis.num_exps
        basis.coefs[ia] *= N
    end

    return basis.norm
end

