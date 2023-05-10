include(".//..//one_electron//overlap.jl")
include(".//..//one_electron//one_electron_integral.jl")
function compound_index(i::Int, j::Int, k::Int, l::Int)::Int
    if i > j
        ij = i*(i+1)/2 + j
    else
        ij = j*(j+1)/2 + i
    end

    if k > l
        kl = k*(k+1)/2 + l
    else
        kl = l*(l+1)/2 + k
    end

    if ij > kl
        ijkl = ij*(ij+1)/2 + kl
    else
        ijkl = kl*(kl+1)/2 + ij
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

function electron_repulsion(a::Float64, lmn1::Vector{Float64}, A::Vector{Float64}, b::Float64, lmn2::Vector{Float64}, B::Vector{Float64},c::Float64, lmn3::Vector{Float64}, C::Vector{Float64}, d::Float64, lmn4::Vector{Float64}, D::Vector{Float64})
    """
    Evaluates kinetic energy integral between two Gaussians.
    Returns a float.
    a,b,c,d:   orbital exponent on Gaussian 'a','b','c','d'
    lmn1,lmn2,lmn3,lmn4: vectors containing orbital angular momentum for Gaussian 'a','b','c','d', respectively
    A,B,C,D: vector  containing origin of Gaussian 'a','b','c','d'
    """
    l1, m1, n1 = lmn1
    l2, m2, n2 = lmn2
    l3, m3, n3 = lmn3
    l4, m4, n4 = lmn4
    p = a + b # composite exponent for P (from Gaussians 'a' and 'b')
    q = c + d # composite exponent for Q (from Gaussians 'c' and 'd')
    alpha = p * q / (p + q)
    P1 = gaussian_product_center(a, A[1], b, B[1]) 
    P2 = gaussian_product_center(a, A[2], b, B[2]) 
    P3=gaussian_product_center(a,A[3],b,B[3])
    Q1 = gaussian_product_center(c, C[1], d, D[1]) 
    Q2 = gaussian_product_center(c, C[2], d, D[2]) 
    Q3=gaussian_product_center(c,C[3],d,D[3])
    PQ1=P1-Q1
    PQ2=P2-Q2
    PQ3=P3-Q3
    RPQ = norm_RPC(PQ1,PQ2,PQ3)

    val = 0.0
    for t in 0:convert(Int64,(l1 + l2))
        for u in 0:convert(Int64,(m1 + m2))
            for v in 0:convert(Int64,(n1 + n2))
                for tau in 0:convert(Int64,(l3 + l4))
                    for nu in 0:convert(Int64,(m3 + m4))
                        for phi in 0:convert(Int64,(n3 + n4))
                            val += Expansion_coeff(l1, l2, t, A[1] - B[1], a, b) *
                                Expansion_coeff(m1, m2, u, A[2] - B[2], a, b) *
                                Expansion_coeff(n1, n2, v, A[3] - B[3], a, b) *
                                Expansion_coeff(l3, l4, tau, C[1] - D[1], c, d) *
                                Expansion_coeff(m3, m4, nu, C[2] - D[2], c, d) *
                                Expansion_coeff(n3, n4, phi, C[3] - D[3], c, d) *((-1)^(tau + nu + phi)) *
                                R_aux_Hermite_coloumb(t + tau, u + nu, v + phi, 0, alpha,PQ1,PQ2,PQ3, RPQ)
                        end
                    end
                end
            end
        end
    end

    val *= 2 * (pi^2.5) / (p * q * sqrt(p + q))
    return val
end

function ERI(  aexps::Vector{Float64}, acoefs::Vector{Float64}, ashell::Vector{Float64}, anorm::Vector{Float64}, aorigin::Vector{Float64},  
    bexps::Vector{Float64}, bcoefs::Vector{Float64}, bshell::Vector{Float64}, bnorm::Vector{Float64}, borigin::Vector{Float64},
    cexps::Vector{Float64}, ccoefs::Vector{Float64}, cshell::Vector{Float64}, cnorm::Vector{Float64}, corigin::Vector{Float64},
    dexps::Vector{Float64}, dcoefs::Vector{Float64}, dshell::Vector{Float64}, dnorm::Vector{Float64}, dorigin::Vector{Float64})
    eri=0.0
    for (ja,ca) in enumerate(acoefs)
        for (jb,cb) in enumerate(bcoefs)
            for (jc,cc) in enumerate(ccoefs)
                for (jd,cd) in enumerate(dcoefs)
                    eri += anorm[ja]*bnorm[jb]*cnorm[jc]*dnorm[jd]*ca*cb*cc*cd*electron_repulsion(aexps[ja],ashell,aorigin,bexps[jb],bshell,borigin,cexps[jc],cshell,corigin,dexps[jd],dshell,dorigin)
                end
            end
        end

    end
    return eri
    #println(eri)
end




#a=ERI([130.7093214, 23.80886605, 6.443608313], [0.1543289673, 0.5353281423, 0.4446345422], [0.0, 0.0, 0.0],[27.551167822078394, 7.681819989204459, 2.882417873168662],[0.0, -0.143225816552, 0.0], 
#[5.033151319, 1.169596125, 0.38038896], [-0.09996722919, 0.3995128261, 0.7001154689], [0.0, 0.0, 0.0],[2.394914882501622, 0.8015618386293725, 0.34520813393821864],[0.0, -0.143225816552, 0.0],
#[5.033151319, 1.169596125, 0.38038896], [0.155916275, 0.6076837186, 0.3919573931], [1.0, 0.0, 0.0],[10.745832634231425, 1.7337440707285057, 0.42581893344677013],[0.0, -0.143225816552, 0.0],
#[5.033151319, 1.169596125, 0.38038896], [0.155916275, 0.6076837186, 0.3919573931], [0.0, 1.0, 0.0],[10.745832634231425, 1.7337440707285057, 0.42581893344677013],[0.0, -0.143225816552, 0.0])
#println(a)
#error("gcg")
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
            #println(a,b,c,d)
            exp1, exp2, exp3, exp4 = exps[a], exps[b], exps[c], exps[d]
            #println(exp1,exp2,exp3,exp4)
            coefs1, coefs2, coefs3, coefs4 = coefs[a], coefs[b], coefs[c], coefs[d]
            #println(coefs1,coefs2,coefs3,coefs4)
            norm1, norm2, norm3, norm4 = norms[a], norms[b], norms[c], norms[d]
            #println(norm1,norm2,norm3,norm4)
            #println(origins[a],origins[b],origins[c],origins[d])
            Temp_mat[i] = ERI(exp1, coefs1, shells[a], norm1, origins[a],
                            exp2, coefs2, shells[b], norm2, origins[b],
                            exp3, coefs3, shells[c], norm3, origins[c],
                            exp4, coefs4, shells[d], norm4, origins[d])
        #println(Temp_mat)
    end
        nbasis = length(exps)
        Twoe_mat = zeros(nbasis, nbasis, nbasis, nbasis)
        ij, kl = 0.0, 0.0
        m=1
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
function doERIs(exps::Vector{Any}, coefs::Vector{Any}, origins::Vector{Any}, shells::Vector{Any}, norms::Vector{Any})
    N=length(exps)
    TwoE= zeros(N,N,N,N)
    for i in 1:N
        for j in 1:i
            ij = (i*(i-1)/2 + j)
            for k in 1:N 
                for l in 1:k
                    kl = (k*(k-1)/2 + l)
                    if ij ≥ kl
                        val = ERI(exps[i], coefs[i], shells[i], norms[i], origins[i],
                        exps[j], coefs[j], shells[j], norms[j], origins[j],
                        exps[k], coefs[k], shells[k], norms[k], origins[k],
                        exps[l], coefs[l], shells[l], norms[l], origins[l])
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
