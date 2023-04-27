#include("./../basis.jl")
include(".//..//one_electron//overlap.jl")
include(".//..//one_electron//one_electron_integral.jl")
function compound_index(i::Int, j::Int, k::Int, l::Int)::Int
    if i > j
        ij = i*(i+1)÷2 + j
    else
        ij = j*(j+1)÷2 + i
    end

    if k > l
        kl = k*(k+1)÷2 + l
    else
        kl = l*(l+1)÷2 + k
    end

    if ij > kl
        ijkl = ij*(ij+1)÷2 + kl
    else
        ijkl = kl*(kl+1)÷2 + ij
    end

    return Int(ijkl)
end
function uniqueindex(n::Int)
    maxx = compound_index(n, n, n, n)
    minn = 0.0
    ll = []
    for i in 1:n, j in 1:n, k in 1:n, l in 1:n
        if compound_index(i-1, j-1, k-1, l-1) == minn
            push!(ll, [i, j, k, l])
            minn += 1
        end
        if minn >= maxx
            break
        end
    end
    return ll
end

function electron_repulsion(a::Float64, lmn1::Vector{Any}, A::Vector{Float64}, b::Float64, lmn2::Vector{Any}, B::Vector{Float64},c::Float64, lmn3::Vector{Any}, C::Vector{Float64}, d::Float64, lmn4::Vector{Any}, D::Vector{Float64})
    """
    Evaluates kinetic energy integral between two Gaussians.
    Returns a float.
    a,b,c,d:   orbital exponent on Gaussian 'a','b','c','d'
    lmn1,lmn2,lmn3,lmn4: int tuple containing orbital angular momentum
                   for Gaussian 'a','b','c','d', respectively
    A,B,C,D:   list containing origin of Gaussian 'a','b','c','d'
    """
    l1, m1, n1 = lmn1[1]
    l2, m2, n2 = lmn2[1]
    l3, m3, n3 = lmn3[1]
    l4, m4, n4 = lmn4[1]
    p = a + b # composite exponent for P (from Gaussians 'a' and 'b')
    q = c + d # composite exponent for Q (from Gaussians 'c' and 'd')
    alpha = p * q / (p + q)
    P = gaussian_product_center(a, A, b, B) # A and B composite center
    Q = gaussian_product_center(c, C, d, D) # C and D composite center
    RPQ = norm(P - Q)

    val = 0.0
    for t in 0:convert(Int64,l1 + l2)
        for u in 0:convert(Int64,m1 + m2)
            for v in 0:convert(Int64,n1 + n2)
                for tau in 0:convert(Int64,l3 + l4)
                    for nu in 0:convert(Int64,m3 + m4)
                        for phi in 0:convert(Int64,n3 + n4)
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

function ERI(  aexps::Vector{Float64}, acoefs::Vector{Float64}, ashell::Vector{Any}, anorm::Vector{Float64}, aorigin::Vector{Float64},  
    bexps::Vector{Float64}, bcoefs::Vector{Float64}, bshell::Vector{Any}, bnorm::Vector{Float64}, borigin::Vector{Float64},
    cexps::Vector{Float64}, ccoefs::Vector{Float64}, cshell::Vector{Any}, cnorm::Vector{Float64}, corigin::Vector{Float64},
    dexps::Vector{Float64}, dcoefs::Vector{Float64}, dshell::Vector{Any}, dnorm::Vector{Float64}, dorigin::Vector{Float64})
    
    noa_coeffs = length(acoefs)
    nob_coeffs = length(bcoefs)
    noc_coeffs = length(ccoefs)
    nod_coeffs = length(dcoefs)
    norm1, norm2, norm3, norm4 = 0.0, 0.0, 0.0, 0.0
    coef1, coef2, coef3, coef4 = 0.0, 0.0, 0.0, 0.0
    exp1, exp2, exp3, exp4 = 0.0, 0.0, 0.0, 0.0
    eri = 0.0
    for ja in 1:noa_coeffs
        norm1,coef1,exp1 = anorm[ja],acoefs[ja],aexps[ja]
        for jb in 1:nob_coeffs
            norm2,coef2,exp2 = bnorm[jb],bcoefs[jb],bexps[jb]
            for jc in 1:noc_coeffs
                norm3,coef3,exp3 = cnorm[jc],ccoefs[jc],cexps[jc]
                for jd in 1:nod_coeffs
                    norm4,coef4,exp4 = dnorm[jd],dcoefs[jd],dexps[jd]
                    eri += norm1*norm2*norm3*norm4*coef1*coef2*coef3*coef4*electron_repulsion(exp1,ashell,aorigin,exp2,bshell,borigin,exp3,cshell,corigin,exp4,dshell,dorigin)
                end

            end
        end
    return eri
    #println(eri)
    end
end
function Eri_mat(exps::Vector{Any}, coefs::Vector{Any}, origins::Vector{Any}, shells::Vector{Any}, norms::Vector{Any})
    unique = uniqueindex(length(exps))
    length_unique = length(unique)
    Temp_mat = zeros(length_unique)
    a, b, c, d, i = 0, 0, 0, 0, 0
    exp1, exp2, exp3, exp4 = 0.0, 0.0, 0.0, 0.0
    coefs1, coefs2, coefs3, coefs4 = 0.0, 0.0, 0.0, 0.0
    norm1, norm2, norm3, norm4 = 0.0, 0.0, 0.0, 0.0
    #println(unique)
    for i in 1:length_unique
        a, b, c, d = unique[i][1], unique[i][2], unique[i][3], unique[i][4]
        exp1, exp2, exp3, exp4 = exps[a], exps[b], exps[c], exps[d]
        coefs1, coefs2, coefs3, coefs4 = coefs[a], coefs[b], coefs[c], coefs[d]
        norm1, norm2, norm3, norm4 = norms[a], norms[b], norms[c], norms[d]
        Temp_mat[i] = ERI(exp1, coefs1, shells[a,:], norm1, origins[a],
                           exp2, coefs2, shells[b,:], norm2, origins[b],
                           exp3, coefs3, shells[c,:], norm3, origins[c],
                           exp4, coefs4, shells[d,:], norm4, origins[d])
    end
        nbasis = length(exps)
        Twoe_mat = zeros(nbasis, nbasis, nbasis, nbasis)
        ij, kl = 0.0, 0.0
        m, j, k, l = 1, 1, 1, 1

        for i in 1:nbasis
            for j in 1:i
                for k in 1:nbasis
                    for l in 1:k
                        ij = i * (i + 1) / 2 + j
                        kl = k * (k + 1) / 2 + l
                        if ij >= kl
                            Twoe_mat[i,j,k,l] = Temp_mat[m]
                            Twoe_mat[i,j,l,k] = Temp_mat[m]
                            Twoe_mat[j,i,k,l] = Temp_mat[m]
                            Twoe_mat[j,i,l,k] = Temp_mat[m]
                            Twoe_mat[k,l,i,j] = Temp_mat[m]
                            Twoe_mat[l,k,i,j] = Temp_mat[m]
                            Twoe_mat[k,l,j,i] = Temp_mat[m]
                            Twoe_mat[l,k,j,i] = Temp_mat[m]

                            m += 1
                        end
                    end
                end
            end
        end

    return Twoe_mat, Temp_mat
end
