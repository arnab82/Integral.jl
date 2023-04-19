include("./../basis.jl")
include("./../Intgl.jl")



function dipole(a, lmn1, A, b, lmn2, B, C, direction)
    l1, m1, n1 = lmn1
    l2, m2, n2 = lmn2
    P = gaussian_product_center(a, A, b, B)
    if direction == "x"
        XPC = P[1] - C[1]
        D = E(l1, l2, 1, A[1]-B[1], a, b) + XPC*E(l1, l2, 0, A[1]-B[1], a, b)
        S2 = E(m1, m2, 0, A[2]-B[2], a, b)
        S3 = E(n1, n2, 0, A[3]-B[3], a, b)
        return D*S2*S3*sqrt(pi/(a+b))^3
    elseif direction == "y"
        YPC = P[2] - C[2]
        S1 = E(l1, l2, 0, A[1]-B[1], a, b)
        D = E(m1, m2, 1, A[2]-B[2], a, b) + YPC*E(m1, m2, 0, A[2]-B[2], a, b)
        S3 = E(n1, n2, 0, A[3]-B[3], a, b)
        return S1*D*S3*sqrt(pi/(a+b))^3
    elseif direction == "z"
        ZPC = P[3] - C[3]
        S1 = E(l1, l2, 0, A[1]-B[1], a, b)
        S2 = E(m1, m2, 0, A[2]-B[2], a, b)
        D = E(n1, n2, 1, A[3]-B[3], a, b) + ZPC*E(n1, n2, 0, A[3]-B[3], a, b)
        return S1*S2*D*sqrt(pi/(a+b))^3
    end
end
function Mu(a, b, C, direction)
    mu = 0.0
    for (ia, ca) in enumerate(a.coefs)
        for (ib, cb) in enumerate(b.coefs)
            mu += a.norm[ia]*b.norm[ib]*ca*cb*\
                     dipole(a.exps[ia],a.shell,a.origin,
                     b.exps[ib],b.shell,b.origin,C,direction)
        end
    end
    return mu
end
function RxDel(a, b, C, direction)
    l = 0.0
    for ia in 1:length(a.coefs)
        for ib in 1:length(b.coefs)
            l += a.norm[ia] * b.norm[ib] * a.coefs[ia] * b.coefs[ib] *
                 angular(a.exps[ia], a.shell, a.origin,
                         b.exps[ib], b.shell, b.origin, C, direction)
        end
    end
    return l
end