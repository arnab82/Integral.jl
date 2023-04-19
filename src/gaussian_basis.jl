

mutable struct Basis_Function
    origin::Array{Float64,1}
    shell::NTuple{3,Int}
    exps::Vector{Float64}
    coefs::Vector{Float64}
    num_exps::Int
    norm::Vector{Float64}
    
    function Basis_Function(origin::Array{Float64,1}, shell::NTuple{3,Int}, exps::Vector{Float64}, coefs::Vector{Float64})
        num_exps = length(exps)
        norm = zeros(num_exps)
        normalize!(exps, coefs, shell, norm)
        new(origin, shell, exps, coefs, num_exps, norm)
    end
end

function normalize!(exps::Vector{Float64}, coefs::Vector{Float64}, shell::NTuple{3,Int}, norm::Vector{Float64})
    l, m, n = shell
    L = l + m + n
    # `norm` is a vector of length equal to number primitives
    # normalize primitives first (PGBFs)
    norm .= sqrt.(2^(2*(l+m+n)+1.5) .* exps.^(l+m+n+1.5) ./ (fact2(2*l-1) .* fact2(2*m-1) .* fact2(2*n-1) .* pi^1.5))

    # now normalize the contracted basis functions (CGBFs)
    # Eq. 1.44 of Valeev integral whitepaper
    prefactor = pi^1.5 * fact2(2*l - 1) * fact2(2*m - 1) * fact2(2*n - 1) / 2.0^L

    N = 0.0
    for ia in 1:num_exps
        for ib in 1:num_exps
            N += norm[ia]*norm[ib]*coefs[ia]*coefs[ib] / (exps[ia] + exps[ib])^(L+1.5)
        end
    end

    N *= prefactor
    N = N^(-0.5)
    coefs .*= N
    nothing
end
