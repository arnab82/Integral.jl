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
function E(i::Float64, j::Float64, t::Int64, Qx::Float64, a::Float64, b::Float64)
    """
    Recursive definition of Hermite Gaussian coefficients.
    Returns a float.
    a: orbital exponent on Gaussian 'a'
    b: orbital exponent on Gaussian 'b' 
    i, j: orbital angular momentum number on Gaussian 'a' and 'b'
    t: number nodes in Hermite (depends on type of integral,
       e.g. always zero for overlap integrals)
    Qx: distance between origins of Gaussian 'a' and 'b'
    """
    p = a + b
    q = (a * b )/ p
    e = 2.718281828459045
    if (t < 0) || (t > (i + j))
        return 0.0
    elseif i == j == 0.0 && t==0
        # base case
        return e^(-q * Qx * Qx) # K_AB
    elseif j == 0
        # decrement index i
        return (1.0 / (2 * p)) * E(i - 1, j, t - 1, Qx, a, b) -((q * Qx / a) * E(i - 1, j, t, Qx, a, b)) +((t + 1) * E(i - 1, j, t + 1, Qx, a, b))
    else
        # decrement index j
        return (1 / (2 * p)) * E(i, j - 1, t - 1, Qx, a, b) +
               (q * Qx / b) * E(i, j - 1, t, Qx, a, b) +
               (t + 1) * E(i, j - 1, t + 1, Qx, a, b)
    end
end

function overlap(a::Float64, lmn1::Vector{Float64}, A::Vector{Float64}, b::Float64, lmn2::Vector{Float64}, B::Vector{Float64})
    """
    Evaluates overlap integral between two Gaussians
    Returns a float.
    a:    orbital exponent on Gaussian 'a' 
    b:    orbital exponent on Gaussian 'b' 
    lmn1: tuple containing orbital angular momentum (e.g. (1,0,0))
          for Gaussian 'a'
    lmn2: tuple containing orbital angular momentum for Gaussian 'b'
    A:    vector containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
    B:    vector containing origin of Gaussian 'b'
    """
    l1, m1, n1 = lmn1 # shell angular momentum on Gaussian 'a'
    l2, m2, n2 = lmn2 # shell angular momentum on Gaussian 'b'
    S1 = E(l1, l2, 0, A[1] - B[1], a, b) # X
    S2 = E(m1, m2, 0, A[2] - B[2], a, b) # Y
    S3 = E(n1, n2, 0, A[3] - B[3], a, b) # Z
    return S1 * S2 * S3 * ((pi / (a + b))^(3/2))
end

function S(aexps::Vector{Float64}, acoefs::Vector{Float64}, ashell::Vector{Float64}, anorm::Vector{Float64}, aorigin::Vector{Float64},
    bexps::Vector{Float64}, bcoefs::Vector{Float64}, bshell::Vector{Float64}, bnorm::Vector{Float64}, borigin::Vector{Float64})

    noa_coeffs = length(acoefs)
    nob_coeffs = length(bcoefs)
    s = 0.0

    for (ia,ca) in enumerate(acoefs)
        for (ib,cb) in enumerate(bcoefs)
            #println("the typeof of ashell is",typeof(ashell))
            #println(ashell)
            s += anorm[ia] * bnorm[ib] *ca *cb* overlap(aexps[ia],ashell, aorigin,bexps[ib],bshell, borigin)
        end
    end
    return s
end
function S_mat(exps::Array{Any}, coefs::Vector{Any}, origins::Vector{Any}, shells::Vector{Any}, norms::Vector{Any})
    nbasis = length(exps)
    smat = zeros(nbasis, nbasis)
    smat_view = view(smat, :, :)
    s = 0.0
    for i in 1:nbasis
        for j in 1:i
            if j == i
                s = 1.0
                smat_view[i,j]=smat_view[j,i]=s
            else
                s = S(exps[i], coefs[i], shells[i], norms[i], origins[i], exps[j], coefs[j], shells[j], norms[j], origins[j])
                #println("the shells are",shells[i, :],shells[j, :])
                smat_view[i, j] = smat_view[j, i] = s
            end
        end
    end
    
    return smat
    display(smat)
end