include("./../basis.jl")
include("./../Intgl.jl")
function ERI(a::Basis, b::Basis, c::Basis, d::Basis)
    eri = 0.0
    for ja = 1:a.num_exps
        for jb = 1:b.num_exps
            for jc = 1:c.num_exps
                for jd = 1:d.num_exps
                    eri += a.norm[ja]*b.norm[jb]*c.norm[jc]*d.norm[jd]*
                            a.coefs[ja]*b.coefs[jb]*c.coefs[jc]*d.coefs[jd]*
                            electron_repulsion(a.exps[ja],a.shell,a.origin,
                                                b.exps[jb],b.shell,b.origin,
                                                c.exps[jc],c.shell,c.origin,
                                                d.exps[jd],d.shell,d.origin)
                end
            end
        end
    end
    return eri
end

function doERIs(N::Int64, TwoE::Array{Float64, 4}, bfs::Array{Float64, 2})
    for i in 1:N
        for j in 1:i
            ij = (i*(i-1)÷2 + j)
            for k in 1:N
                for l in 1:k
                    kl = (k*(k-1)÷2 + l)
                    if ij ≥ kl
                        val = ERI(bfs[i,:], bfs[j,:], bfs[k,:], bfs[l,:])
                        TwoE[i,j,k,l] = val
                        TwoE[k,l,i,j] = val
                        TwoE[j,i,l,k] = val
                        TwoE[l,k,j,i] = val
                        TwoE[j,i,k,l] = val
                        TwoE[l,k,i,j] = val
                        TwoE[i,j,l,k] = val
                        TwoE[k,l,j,i] = val
                    end
                end
            end
        end
    end
    return TwoE
end
function electron_repulsion(a::Float64, lmn1::Array{Int64,1}, A::Array{Float64,1}, 
    b::Float64, lmn2::Array{Int64,1}, B::Array{Float64,1}, 
    c::Float64, lmn3::Array{Int64,1}, C::Array{Float64,1}, 
    d::Float64, lmn4::Array{Int64,1}, D::Array{Float64,1})::Float64
    l1, m1, n1 = lmn1
    l2, m2, n2 = lmn2
    l3, m3, n3 = lmn3
    l4, m4, n4 = lmn4
    if (a < 0) || (b < 0) || (c < 0) || (d < 0)
        throw(ArgumentError("All exponents must be non-negative"))
    end
    if (l1 < 0) || (l2 < 0) || (l3 < 0) || (l4 < 0)
        throw(ArgumentError("All angular momentum values must be non-negative"))
    end
    if length(A) != 3 || length(B) != 3 || length(C) != 3 || length(D) != 3
        throw(DimensionMismatch("All coordinate vectors must have length 3"))
    end
    
    


    p = a + b
    q = c + d
    alpha = p * q / (p + q)
    Px = (a * A[1] + b * B[1]) / p
    Py = (a * A[2] + b * B[2]) / p
    Pz = (a * A[3] + b * B[3]) / p
    Qx = (c * C[1] + d * D[1]) / q
    Qy = (c * C[2] + d * D[2]) / q
    Qz = (c * C[3] + d * D[3]) / q
    RPQ = sqrt((Px - Qx)^2 + (Py - Qy)^2 + (Pz - Qz)^2)

    val = 0.0
    for t = 0:(l1 + l2)
        for u = 0:(m1 + m2)
            for v = 0:(n1 + n2)
                for tau = 0:(l3 + l4)
                    for nu = 0:(m3 + m4)
                        for phi = 0:(n3 + n4)
                            val += E(l1, l2, t, A[1] - B[1], a, b) * \
                                E(m1, m2, u, A[2] - B[2], a, b) * \
                                E(n1, n2, v, A[3] - B[3], a, b) * \
                                E(l3, l4, tau, C[1] - D[1], c, d) * \
                                E(m3, m4, nu, C[2] - D[2], c, d) * \
                                E(n3, n4, phi, C[3] - D[3], c, d) * \
                                (-1)^(tau + nu + phi) * \
                                R(t + tau, u + nu, v + phi, 0, alpha, Px - Qx, Py - Qy, Pz - Qz, RPQ)
                        end
                    end
                end
            end
        end
    end

    val *= 2 * (pi^2.5) / (p * q * sqrt(p + q))
    return val
end


