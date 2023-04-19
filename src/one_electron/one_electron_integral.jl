include("./../Intgl.jl")
include("./../basis.jl")


using LinearAlgebra

function T(a ,b)::Float64
    t = 0.0
    for (ia, ca) in enumerate(a.coefs)
        for (ib, cb) in enumerate(b.coefs)
            t += a.norm[ia]*b.norm[ib]*ca*cb*
                     kinetic(a.exps[ia],a.shell,a.origin,
                     b.exps[ib],b.shell,b.origin)
        end
    end
    return t
end
function V(a, b, C::Vector{Float64})
    """ Nuclear attraction integrals """
    v = 0.0
    for (ia, ca) in enumerate(a.coefs)
        for (ib, cb) in enumerate(b.coefs)
            v += a.norm[ia]*b.norm[ib]*ca*cb*
                     nuclear_attraction(a.exps[ia],a.shell,a.origin,
                     b.exps[ib],b.shell,b.origin,C)
        end
    end
    return v
end
function kinetic(a, lmn1, A, b, lmn2, B)
    # explicit kinetic in terms of "E" operator
    # generalized to include GIAO derivatives
    l1, m1, n1 = lmn1
    l2, m2, n2 = lmn2
    Ax, Ay, Az = (2 .* [lmn2...] .+ 1) .* b
    Bx = By = Bz = -2 .* b.^2  # redundant, I know
    Cx, Cy, Cz = -0.5 .* [lmn2...] .* ([lmn2...] .- 1)

    Tx = Ax .* E(l1, l2  , 0, A[1]-B[1], a, b) .+ \
         Bx .* E(l1, l2+2, 0, A[1]-B[1], a, b) .+ \
         Cx .* E(l1, l2-2, 0, A[1]-B[1], a, b)
    Tx .= Tx .* E(m1, m2, 0, A[2]-B[2], a, b)
    Tx .= Tx .* E(n1, n2, 0, A[3]-B[3], a, b)

    Ty = Ay .* E(m1, m2  , 0, A[2]-B[2], a, b) .+ \
         By .* E(m1, m2+2, 0, A[2]-B[2], a, b) .+ \
         Cy .* E(m1, m2-2, 0, A[2]-B[2], a, b)
    Ty .= Ty .* E(l1, l2, 0, A[1]-B[1], a, b)
    Ty .= Ty .* E(n1, n2, 0, A[3]-B[3], a, b)

    Tz = Az .* E(n1, n2  , 0, A[3]-B[3], a, b) .+ \
         Bz .* E(n1, n2+2, 0, A[3]-B[3], a, b) .+ \
         Cz .* E(n1, n2-2, 0, A[3]-B[3], a, b)
    Tz .= Tz .* E(l1, l2, 0, A[1]-B[1], a, b)
    Tz .= Tz .* E(m1, m2, 0, A[2]-B[2], a, b)

    return (Tx + Ty + Tz) .* (pi ./ (a+b)).^1.5
end

function angular(a, lmn1, A, b, lmn2, B, C, direction)
    l1, m1, n1 = lmn1
    l2, m2, n2 = lmn2
    P = gaussian_product_center(a, A, b, B)

    XPC = P[1] - C[1]
    YPC = P[2] - C[2]
    ZPC = P[3] - C[3]

    S0x = E(l1, l2, 0, A[1]-B[1], a, b)
    S0y = E(m1, m2, 0, A[2]-B[2], a, b)
    S0z = E(n1, n2, 0, A[3]-B[3], a, b)

    S1x = E(l1, l2, 0, A[1]-B[1], a, b, 1, A[1]-C[1]) + XPC * E(l1, l2, 0, A[1]-B[1], a, b)
    S1y = E(m1, m2, 0, A[2]-B[2], a, b, 1, A[2]-C[2]) + YPC * E(m1, m2, 0, A[2]-B[2], a, b)
    S1z = E(n1, n2, 0, A[3]-B[3], a, b, 1, A[3]-C[3]) + ZPC * E(n1, n2, 0, A[3]-B[3], a, b)

    D1x = l2 * E(l1, l2-1, 0, A[1]-B[1], a, b) - 2b * E(l1, l2+1, 0, A[1]-B[1], a, b)
    D1y = m2 * E(m1, m2-1, 0, A[2]-B[2], a, b) - 2b * E(m1, m2+1, 0, A[2]-B[2], a, b)
    D1z = n2 * E(n1, n2-1, 0, A[3]-B[3], a, b) - 2b * E(n1, n2+1, 0, A[3]-B[3], a, b)

    if lowercase(direction) == "x"
        return -S0x * (S1y * D1z - S1z * D1y) * power(pi/(a+b), 1.5)
    elseif lowercase(direction) == "y"
        return -S0y * (S1z * D1x - S1x * D1z) * power(pi/(a+b), 1.5)
    elseif lowercase(direction) == "z"
        return -S0z * (S1x * D1y - S1y * D1x) * power(pi/(a+b), 1.5)
    else
        error("Invalid direction specified.")
    end
end
function nuclear_attraction(a, lmn1, A, b, lmn2, B, C)
    """ Returns nuclear attraction integral between two primitive Gaussians"""
    l1, m1, n1 = lmn1
    l2, m2, n2 = lmn2
    p = a + b
    P = gaussian_product_center(a, A, b, B)
    RPC = norm(P - C)

    val = 0.0
    for t in 0:(l1 + l2)
        for u in 0:(m1 + m2)
            for v in 0:(n1 + n2)
                val += E(l1, l2, t, A[1] - B[1], a, b) *
                       E(m1, m2, u, A[2] - B[2], a, b) *
                       E(n1, n2, v, A[3] - B[3], a, b) *
                       R(t, u, v, 0, p, P[1] - C[1], P[2] - C[2], P[3] - C[3], RPC)
            end
        end
    end
    val *= 2 * pi / p  # Pink book, Eq(9.9.40)
    return val
end
