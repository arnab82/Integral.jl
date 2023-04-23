
using SpecialFunctions
using Combinatorics: doublefactorial
using LinearAlgebra: norm, eigen
using StaticArrays
const ang2bohr = 1.8897261246257702

function fact2(n::Float64)
    if n == 0
        return 1
    end
    return prod(collect(range(1,n,2)))
end

function E(i, j, t, Qx, a, b)
    """
    Recursive definition of Hermite Gaussian coefficients.
    Returns a float.
    a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
    b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
    i, j: orbital angular momentum number on Gaussian 'a' and 'b'
    t: number nodes in Hermite (depends on type of integral,
       e.g. always zero for overlap integrals)
    Qx: distance between origins of Gaussian 'a' and 'b'
    """
    p = a + b
    q = a * b / p
    if (t < 0) || (t > (i + j))
        # out of bounds for t
        return 0.0
    elseif i == j == t == 0
        # base case
        return exp(-q * Qx * Qx) # K_AB
    elseif j == 0
        # decrement index i
        return (1 / (2 * p)) * E(i - 1, j, t - 1, Qx, a, b) -
               (q * Qx / a) * E(i - 1, j, t, Qx, a, b) +
               (t + 1) * E(i - 1, j, t + 1, Qx, a, b)
    else
        # decrement index j
        return (1 / (2 * p)) * E(i, j - 1, t - 1, Qx, a, b) +
               (q * Qx / b) * E(i, j - 1, t, Qx, a, b) +
               (t + 1) * E(i, j - 1, t + 1, Qx, a, b)
    end
end
function overlap(a, lmn1, A, b, lmn2, B)
    """
    Evaluates overlap integral between two Gaussians
    Returns a float.
    a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)
    b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)
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
    return S1 * S2 * S3 * (pi / (a + b))^(3/2)
end


function S(a, b)
    """
    Evaluates overlap between two contracted Gaussians
    Returns float.
    Arguments:
    a: contracted Gaussian 'a', BasisFunction object
    b: contracted Gaussian 'b', BasisFunction object
    """
    s = 0.0
    for (ia, ca) in enumerate(a.coefs)
        for (ib, cb) in enumerate(b.coefs)
            s += a.norm[ia] * b.norm[ib] * ca * cb *
                 overlap(a.exps[ia], a.shell, a.origin,
                         b.exps[ib], b.shell, b.origin)
        end
    end
    return s
end
struct BasisFunction
    """
    A struct that contains all our basis function data
    Attributes:
    origin: array/list containing the coordinates of the Gaussian origin
    shell: tuple of angular momentum
    exps: list of primitive Gaussian exponents
    coefs: list of primitive Gaussian coefficients
    num_exps: number of primitive Gaussian exponents
    norm: list of normalization factors for Gaussian primitives
    """
    origin::Array{Float64}
    shell::Array{Float64}
    exps::Array{Float64}
    coefs::Array{Float64}
    num_exps::Int64
    norms::Array{Float64}
end

function BasisFunction(; origin=[0.0, 0.0, 0.0], shell=(0, 0, 0), exps=[], coefs=[], num_exps=0)
    if num_exps == 0
        num_exps = length(exps)
    end
    bf = BasisFunction(origin, shell, exps, coefs, num_exps, [0.0])
    normalize_basis!(bf)
    return bf
end

function normalize_basis!(bf::BasisFunction)
    """
    Routine to normalize the basis functions, in case they do not integrate to unity.
    """
    l, m, n = bf.shell
    println(l,m,n)
    L = l + m + n
    println(L)
    println(bf.exps)
    println(fact2(2*l-1))

    # norm is an array of length equal to number primitives
    # normalize primitives first (PGBFs)
    bf.norm = sqrt.(2.0^(2*(l + m + n) + 1.5) .* bf.exps.^(l + m + n + 1.5) /(fact2(2*l - 1) * fact2(2*m - 1) * fact2(2*n - 1) * pi^1.5))

    # now normalize the contracted basis functions (CGBFs)
    # Eq. 1.44 of Valeev integral whitepaper
    prefactor = (pi^1.5 * fact2(2*l - 1) * fact2(2*m - 1) * fact2(2*n - 1)) / 2.0^L

    N = 0.0
    for ia in 1:bf.num_exps
        for ib in 1:bf.num_exps
            N += bf.norm[ia] * bf.norm[ib] * bf.coefs[ia] * bf.coefs[ib] /
                 (bf.exps[ia] + bf.exps[ib])^(L + 1.5)
        end
    end

    N *= prefactor
    N = 1 / sqrt(N)
    bf.coefs .*= N
end

function kinetic(a, lmn1, A, b, lmn2, B)
    """
    Evaluates kinetic energy integral between two Gaussians. Returns a float.
    
    a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
    b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
    lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
          for Gaussian 'a'
    lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
    A: list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
    B: list containing origin of Gaussian 'b'
    """
    l1, m1, n1 = lmn1
    l2, m2, n2 = lmn2
    term0 = b * (2 * (l2 + m2 + n2) + 3) *
            overlap(a, (l1, m1, n1), A, b, (l2, m2, n2), B)
    term1 = -2 * b^2 *
            (overlap(a, (l1, m1, n1), A, b, (l2 + 2, m2, n2), B) +
             overlap(a, (l1, m1, n1), A, b, (l2, m2 + 2, n2), B) +
             overlap(a, (l1, m1, n1), A, b, (l2, m2, n2 + 2), B))
    term2 = -0.5 * (l2 * (l2 - 1) *
             overlap(a, (l1, m1, n1), A, b, (l2 - 2, m2, n2), B) +
             m2 * (m2 - 1) *
             overlap(a, (l1, m1, n1), A, b, (l2, m2 - 2, n2), B) +
             n2 * (n2 - 1) *
             overlap(a, (l1, m1, n1), A, b, (l2, m2, n2 - 2), B))
    return term0 + term1 + term2
end
function T(a::BasisFunction, b::BasisFunction)
    """
    Evaluates kinetic energy between two contracted Gaussians
    Returns a float.
    
    Arguments:
    a: contracted Gaussian 'a', BasisFunction object
    b: contracted Gaussian 'b', BasisFunction object
    """
    t = 0.0
    for (ia, ca) in enumerate(a.coefs)
        for (ib, cb) in enumerate(b.coefs)
            t += a.norm[ia] * b.norm[ib] * ca * cb *
                 kinetic(a.exps[ia], a.shell, a.origin,
                         b.exps[ib], b.shell, b.origin)
        end
    end
    return t
end
function R(t, u, v, n, p, PCx, PCy, PCz, RPC)
    """
    Returns the Coulomb auxiliary Hermite integrals
    Returns a float.
    Arguments:
    t, u, v:   order of Coulomb Hermite derivative in x, y, z
               (see defs in Helgaker and Taylor)
    n:         order of Boys function
    PCx, PCy, PCz: Cartesian vector distance between Gaussian
                   composite center P and nuclear center C
    RPC:       Distance between P and C
    """
    T = p * RPC * RPC
    val = 0.0
    if t == u == v == 0
        val += (-2 * p)^n * boys(n, T)
    elseif t == u == 0
        if v > 1
            val += (v - 1) * R(t, u, v - 2, n + 1, p, PCx, PCy, PCz, RPC)
        end
        val += PCz * R(t, u, v - 1, n + 1, p, PCx, PCy, PCz, RPC)
    elseif t == 0
        if u > 1
            val += (u - 1) * R(t, u - 2, v, n + 1, p, PCx, PCy, PCz, RPC)
        end
        val += PCy * R(t, u - 1, v, n + 1, p, PCx, PCy, PCz, RPC)
    else
        if t > 1
            val += (t - 1) * R(t - 2, u, v, n + 1, p, PCx, PCy, PCz, RPC)
        end
        val += PCx * R(t - 1, u, v, n + 1, p, PCx, PCy, PCz, RPC)
    end
    return val
end
function boys(n, T)
    return hyp1f1(n + 0.5, n + 1.5, -T) / (2.0 * n + 1.0)
end

function gaussian_product_center(a, A, b, B)
    return (a * A + b * B) / (a + b)
end

function nuclear_attraction(a, lmn1, A, b, lmn2, B, C)
    l1, m1, n1 = lmn1 
    l2, m2, n2 = lmn2
    p = a + b
    P = gaussian_product_center(a, A, b, B) # Gaussian composite center
    RPC = norm(P - C)

    val = 0.0
    for t in 0:(l1 + l2)
        for u in 0:(m1 + m2)
            for v in 0:(n1 + n2)
                val += E(l1, l2, t, A[1] - B[1], a, b) * \
                       E(m1, m2, u, A[2] - B[2], a, b) * \
                       E(n1, n2, v, A[3] - B[3], a, b) * \
                       R(t, u, v, 0, p, P[1] - C[1], P[2] - C[2], P[3] - C[3], RPC)
            end
        end
    end
    val *= 2 * π / p
    return val
end

function V(a::BasisFunction, b::BasisFunction, C::Vector{Float64})
    """
    Evaluates overlap between two contracted Gaussians.
    Returns a float.
    Arguments:
    a: contracted Gaussian 'a', BasisFunction object
    b: contracted Gaussian 'b', BasisFunction object
    C: center of nucleus
    """
    v = 0.0
    for (ia, ca) in enumerate(a.coefs)
        for (ib, cb) in enumerate(b.coefs)
            v += a.norm[ia] * b.norm[ib] * ca * cb *
                nuclear_attraction(a.exps[ia], a.shell, a.origin,
                                    b.exps[ib], b.shell, b.origin, C)
        end
    end
    return v
end
function electron_repulsion(a, lmn1, A, b, lmn2, B, c, lmn3, C, d, lmn4, D)
    """
    Evaluates kinetic energy integral between two Gaussians.
    Returns a float.
    a,b,c,d:   orbital exponent on Gaussian 'a','b','c','d'
    lmn1,lmn2,lmn3,lmn4: int tuple containing orbital angular momentum
                   for Gaussian 'a','b','c','d', respectively
    A,B,C,D:   list containing origin of Gaussian 'a','b','c','d'
    """
    l1, m1, n1 = lmn1
    l2, m2, n2 = lmn2
    l3, m3, n3 = lmn3
    l4, m4, n4 = lmn4
    p = a + b # composite exponent for P (from Gaussians 'a' and 'b')
    q = c + d # composite exponent for Q (from Gaussians 'c' and 'd')
    alpha = p * q / (p + q)
    P = gaussian_product_center(a, A, b, B) # A and B composite center
    Q = gaussian_product_center(c, C, d, D) # C and D composite center
    RPQ = norm(P - Q)

    val = 0.0
    for t in 0:(l1 + l2)
        for u in 0:(m1 + m2)
            for v in 0:(n1 + n2)
                for tau in 0:(l3 + l4)
                    for nu in 0:(m3 + m4)
                        for phi in 0:(n3 + n4)
                            val += E(l1, l2, t, A[1] - B[1], a, b) *
                                E(m1, m2, u, A[2] - B[2], a, b) *
                                E(n1, n2, v, A[3] - B[3], a, b) *
                                E(l3, l4, tau, C[1] - D[1], c, d) *
                                E(m3, m4, nu, C[2] - D[2], c, d) *
                                E(n3, n4, phi, C[3] - D[3], c, d) *
                                (-1)^(tau + nu + phi) *
                                R(t + tau, u + nu, v + phi, 0, alpha,
                                    P[1] - Q[1], P[2] - Q[2], P[3] - Q[3], RPQ)
                        end
                    end
                end
            end
        end
    end

    val *= 2 * (pi^2.5) / (p * q * sqrt(p + q))
    return val
end
function ERI(a::BasisFunction, b::BasisFunction, c::BasisFunction, d::BasisFunction)
    eri = 0.0
    for ja in eachindex(a.coefs)
        for jb in eachindex(b.coefs)
            for jc in eachindex(c.coefs)
                for jd in eachindex(d.coefs)
                    eri += a.norm[ja]*b.norm[jb]*c.norm[jc]*d.norm[jd]*\
                             a.coefs[ja]*b.coefs[jb]*c.coefs[jc]*d.coefs[jd]*
                             electron_repulsion(a.exps[ja], a.shell, a.origin,b.exps[jb], b.shell, b.origin,c.exps[jc], c.shell, c.origin,d.exps[jd], d.shell, d.origin)
                end
            end
        end
    end
    return eri
end
