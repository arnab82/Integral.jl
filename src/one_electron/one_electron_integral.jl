include("./overlap.jl")
#include("./../basis.jl")
using SpecialFunctions
using Combinatorics: doublefactorial
using LinearAlgebra: norm, eigen
using StaticArrays
const ang2bohr = 1.8897261246257702
using PyCall
sp=pyimport("scipy.special")
hyp1f1=sp.hyp1f1
fact2=sp.factorial2


function kinetic(a::Float64, lmn1::Vector{Float64}, A::Vector{Float64}, b::Float64, lmn2::Vector{Float64}, B::Vector{Float64})
    """
    Evaluate kinetic energy integral between two Gaussians
    Returns a float.
    a:    orbital exponent on Gaussian 'a' 
    b:    orbital exponent on Gaussian 'b' 
    lmn1: vector containing orbital angular momentum (e.g. (1,0,0))
          for Gaussian 'a'
    lmn2: vector containing orbital angular momentum for Gaussian 'b'
    A:    vector containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
    B:    vector containing origin of Gaussian 'b'
    C:    vector containing origin of nuclear center 'C'
    """
    l1, m1, n1 = lmn1
    l2, m2, n2 = lmn2
    term0 = b * (2 * (l2 + m2 + n2) + 3) *
            overlap(a, lmn1, A, b, lmn2, B)
    term1 = -2 * b^2 *
            (overlap(a, [l1, m1, n1], A, b, [l2 + 2, m2, n2], B) +
             overlap(a, [l1, m1, n1], A, b, [l2, m2 + 2, n2], B) +
             overlap(a, [l1, m1, n1], A, b, [l2, m2, n2 + 2], B))
    term2 = -0.5 * (l2 * (l2 - 1) *
             overlap(a, [l1, m1, n1], A, b, [l2 - 2, m2, n2], B) +
             m2 * (m2 - 1) *
             overlap(a, [l1, m1, n1], A, b, [l2, m2 - 2, n2], B) +
             n2 * (n2 - 1) *
             overlap(a, [l1, m1, n1], A, b, [l2, m2, n2 - 2], B))
    return term0 + term1 + term2
end
function T(aexps::Vector{Float64}, acoefs::Vector{Float64}, ashell::Vector{Float64}, anorm::Vector{Float64}, aorigin::Vector{Float64},
    bexps::Vector{Float64}, bcoefs::Vector{Float64}, bshell::Vector{Float64}, bnorm::Vector{Float64}, borigin::Vector{Float64})
    """
    Evaluates kinetic energy between two contracted Gaussians
    Returns a float.
    
    Arguments:
    aexps,bexps: exponenets of contracted Gaussian 'a'and 'b', BasisFunction object
    acoefs,bcoefs: coefficients of contracted Gaussian 'a'and 'b', BasisFunction object
    ashell,bshell: exponenets of contracted Gaussian 'a'and 'b', BasisFunction object
    anorm,bnorm: norms of contracted Gaussian 'a'and 'b', BasisFunction object
    aorigin,borigin: origins of contracted Gaussian 'a'and 'b', BasisFunction object
    """
    t = 0.0
    for (ia, ca) in enumerate(acoefs)
        for (ib, cb) in enumerate(bcoefs)
            coef1 = acoefs[ia]
            coef2 = bcoefs[ib]
            t += anorm[ia] * bnorm[ib] * ca * cb *
                 kinetic(aexps[ia], ashell, aorigin,
                         bexps[ib], bshell, borigin)
        end
    end
    return t
end

function T_mat(exps::Array{Any}, coefs::Vector{Any}, origins::Vector{Any}, shells::Vector{Any}, norms::Vector{Any})
    nbasis = length(exps)
    tmat = zeros(nbasis, nbasis)
    tmat_view = view(tmat, :, :)
    s = 0.0
    for i in 1:nbasis
        for j in 1:i
            s = T(exps[i], coefs[i], shells[i], norms[i], origins[i], exps[j], coefs[j], shells[j], norms[j], origins[j])
            tmat_view[i, j] = tmat_view[j, i] = s
            
        end
    end
    
    return tmat
    display(tmat)
end



function R_aux_Hermite_coloumb(t, u, v, n, p, PCx, PCy, PCz, RPC)
    """
    Returns the Coulomb auxiliary Hermite integrals
    Returns a float.
    Arguments:
    t, u, v:   order of Coulomb Hermite derivative in x, y, z
    n:         order of Boys function
    PCx, PCy, PCz: Cartesian vector distance between Gaussian
                   composite center P and nuclear center C
    RPC:       Distance between P and C
    """
    T = p * RPC * RPC
    val = 0.0
    if t == u == v == 0
        val += ((-2 * p)^n) * boys(n, T)
    elseif t == u == 0
        if v > 1
            val += (v - 1) * R_aux_Hermite_coloumb(t, u, v - 2, n + 1, p, PCx, PCy, PCz, RPC)
        end
        val += PCz * R_aux_Hermite_coloumb(t, u, v - 1, n + 1, p, PCx, PCy, PCz, RPC)
    elseif t == 0
        if u > 1
            val += (u - 1) * R_aux_Hermite_coloumb(t, u - 2, v, n + 1, p, PCx, PCy, PCz, RPC)
        end
        val += PCy * R_aux_Hermite_coloumb(t, u - 1, v, n + 1, p, PCx, PCy, PCz, RPC)
    else
        if t > 1
            val += (t - 1) * R_aux_Hermite_coloumb(t - 2, u, v, n + 1, p, PCx, PCy, PCz, RPC)
        end
        val += PCx * R_aux_Hermite_coloumb(t - 1, u, v, n + 1, p, PCx, PCy, PCz, RPC)
    end
    return val
end
function boys(n, T)
    "Boys Function"
    return hyp1f1(n + 0.5, n + 1.5, -T) / (2.0 * n + 1.0)
end

function gaussian_product_center(a, A, b, B)
    return (a * A + b * B) / (a + b)
end
function norm_RPC(a::Float64, b::Float64, c::Float64)::Float64
    s = a^2 + b^2 + c^2
    return sqrt(s)
end

function nuclear_attraction(a::Float64, lmn1::Vector{Float64}, A::Vector{Float64}, b::Float64, lmn2::Vector{Float64}, B::Vector{Float64}, C)
    """
    Evaluates nuclear attraction integral between two Gaussians
    Returns a float.
    a:    orbital exponent on Gaussian 'a' 
    b:    orbital exponent on Gaussian 'b' 
    lmn1: vector containing orbital angular momentum for Gaussian 'a'
    lmn2: vector containing orbital angular momentum for Gaussian 'b'
    A:    vector containing origin of Gaussian 'a'
    B:    vector containing origin of Gaussian 'b'
    C:    vector containing origin of nuclear center 'C'
    """
    l1, m1, n1 = lmn1 
    l2, m2, n2 = lmn2
    p = a + b
    RPCdash1 = C[1] - gaussian_product_center(a, A[1], b, B[1])
    RPCdash2 = C[2] - gaussian_product_center(a, A[2], b, B[2])
    RPCdash3 = C[3] - gaussian_product_center(a, A[3], b, B[3])
    #println(P) # Gaussian composite center
    RPC = norm_RPC(RPCdash1, RPCdash2, RPCdash3)
    #println(RPC)
    val = 0.0
    for t in 0:convert(Int64,(l1+l2))
        for u in 0:convert(Int64,(m1+m2))
            for v in 0:convert(Int64,(n1+n2))
                val += Expansion_coeff(l1, l2, t, A[1] - B[1], a, b) * 
                       Expansion_coeff(m1, m2, u, A[2] - B[2], a, b) * 
                       Expansion_coeff(n1, n2, v, A[3] - B[3], a, b) * 
                       R_aux_Hermite_coloumb(t, u, v, 0, p,  -RPCdash1, -RPCdash2, -RPCdash3, RPC)
            end
        end
    end
    val *= 2*π / p 
    return val
end


function V(aexps::Vector{Float64}, acoefs::Vector{Float64}, ashell::Vector{Float64}, anorm::Vector{Float64}, aorigin::Vector{Float64},
    bexps::Vector{Float64}, bcoefs::Vector{Float64}, bshell::Vector{Float64}, bnorm::Vector{Float64}, borigin::Vector{Float64}, C)
    noa_coeffs = length(acoefs)
    nob_coeffs = length(bcoefs)
    v = 0.0
    for ia in 1:noa_coeffs
        for ib in 1:nob_coeffs
            norm1 = anorm[ia]
            norm2 = bnorm[ib]
            coef1 = acoefs[ia]
            coef2 = bcoefs[ib]
            exp1 = aexps[ia]
            exp2 = bexps[ib]

            v += norm1 * norm2 * coef1 * coef2 * nuclear_attraction(exp1, ashell, aorigin
                ,exp2,bshell , 
                borigin, C)
        end
    end
    return v
end

function V_mat(exps::Array{Any}, coefs::Vector{Any}, origins::Vector{Any}, shells::Vector{Any}, norms::Vector{Any}, atomic_nos::Vector{Any}, geom)

    nbasis = length(exps)
    no_of_atoms = length(atomic_nos)
    println(no_of_atoms)
    P1 = zeros((nbasis,nbasis))
    P2 = zeros((nbasis, nbasis))
    v = 0.0
    for i in 1:nbasis
        for j in 1:nbasis
            
            for k in 1:no_of_atoms
                v += atomic_nos[k]*V(exps[i], coefs[i], shells[i], norms[i], origins[i], exps[j], coefs[j], shells[j], 
                                    norms[j], origins[j], geom[k])
            P1[i,j] = -v
            end
        end
    end
    for m in 1:nbasis
        for n in 1:nbasis
            if m == 1 && n == 1
                P2[m,n] = P1[1,1]
            elseif n>1
                P2[m,n] = -(P1[m,n-1] - P1[m,n])
            end
        end
        P2[m,1] = P2[1,m]
    end

    return P2
    display(P2)
end