include("./../Intgl.jl")
include("./../basis.jl")

function overlap(a::Int64, lmn1::Tuple{Int64,Int64,Int64}, A::Vector{Float64}, b::Int64, lmn2::Tuple{Int64,Int64,Int64}, B::Vector{Float64})::Float64
    l1, m1, n1 = lmn1
    l2, m2, n2 = lmn2
    S1 = E(l1, l2, 0, A[1] - B[1], a, b)
    S2 = E(m1, m2, 0, A[2] - B[2], a, b)
    S3 = E(n1, n2, 0, A[3] - B[3], a, b)
    return S1 * S2 * S3 * (pi / (a + b))^1.5
end


function S(a::Basis, b::Basis)::Float64
    s = 0.0
    for ia in eachindex(a.coefs)
        for ib in eachindex(b.coefs)
            s += a.norm[ia] * b.norm[ib] * a.coefs[ia] * b.coefs[ib] * overlap(a.exps[ia], a.shell, a.origin, b.exps[ib], b.shell, b.origin)
        end
    end
    return s
end


