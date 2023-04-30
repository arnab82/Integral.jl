include(".//rhf.jl")
using Einsum
function make_fock_uhf(D::Matrix{Float64}, D_V::Matrix{Float64}, H::Matrix{Float64}, twoe::Array{Float64, 4})
    nbasis = size(D)[1]
    fock = zeros(nbasis, nbasis)
    
    for i in 1:nbasis
        for j in 1:nbasis
            for k in 1:nbasis
                for l in 1:nbasis
                    fock[i, j] += (D[l, k] * twoe[i, j, k, l]) - (D_V[l, k] * twoe[i, l, k, j])
                end
            end
        end
    end
    
    return fock + H
end
function scf_energy_uhf(D::Matrix{Float64}, Fock_alpha::Matrix{Float64}, Fock_beta::Matrix{Float64}, D_alpha::Matrix{Float64}, D_beta::Matrix{Float64}, core_h::Matrix{Float64})
    nbasis = size(D)[1]
    new_energy = 0.0 
    for i in 1:nbasis
        for j in 1:nbasis
            new_energy += 0.5 * ((D[j,i] * core_h[i,j]) + (D_alpha[j,i] * Fock_alpha[i,j]) + (D_beta[j,i] * Fock_beta[i,j]))
        end
    end
    return new_energy
end
function scf_uhf(s::Matrix{Float64},T_mat::Matrix{Float64},V_mat::Matrix{Float64},Eri::Array{Float64, 4},nuclear_repul::Float64,n_elec::Int64,mult::Int64)
    """
    Attributes:
        s: overlap matrix
        T_mat: kinetic energy integral
        V_mat: potential energy integral
        eri: two electron integral
        nuclear_repul: Nuclear repulsion energy
        no of electrons
    return:
        scf energy: Float64 
        c:  Mo co-efficient matrix
        fock: Fock Matrix
        nbasis
    """#@einsum temp1[p,j,k,l]=(c[i,p]*eri[i,j,k,l])
    n_alpha=convert(Int64,(n_elec+mult-1)/2)
    println(n_alpha)
    n_beta=(n_elec-n_alpha)
    println(n_beta)
    h,nbasis=size(s)
    h1e=T_mat+V_mat
    fock= zeros(Float64,nbasis,nbasis)
    init_fock=zeros(Float64,nbasis,nbasis)
    s_half=make_s_half(s)
    init_fock= (h1e*s_half')
    init_fock=(s_half*init_fock)
    init_fock=(init_fock+init_fock')/2
    Hartree_fock_energy=0.0
    del_e=0.0
    no=n_elec/2
    e,c0 = LinearAlgebra.eigen(init_fock)
    c = (s_half*c0)
    list_e=[0.0,]
    eps_list_alpha=[]
    eps_list_beta=[]
    fock= zeros(Float64,nbasis,nbasis)
    c_alpha=zeros(Float64,nbasis,nbasis)
    c_beta=zeros(Float64,nbasis,nbasis)
    D_alpha=make_density(c,n_alpha)
    D_beta=make_density(c,n_beta)
    D=D_alpha+D_beta
    for n in 1:100
        fock_alpha=make_fock_uhf(D,D_alpha,h1e,Eri)
        fock_beta=make_fock_uhf(D,D_beta,h1e,Eri)
        f_dash_alpha= (fock_alpha*s_half')
        f_dash_alpha=(s_half*f_dash_alpha)
        f_dash_alpha=(f_dash_alpha+transpose(f_dash_alpha))/2
        eps_alpha,c_dash_alpha=LinearAlgebra.eigen(f_dash_alpha)
        f_dash_beta= (fock_beta*s_half')
        f_dash_beta=(s_half*f_dash_beta)
        f_dash_beta=(f_dash_alpha+transpose(f_dash_beta))/2
        eps_beta,c_dash_beta=LinearAlgebra.eigen(f_dash_beta)
        push!(eps_list_alpha,eps_alpha)
        push!(eps_list_beta,eps_beta)
        c_alpha=(s_half*c_dash_alpha)
        c_beta=(s_half*c_dash_beta)
        D_alpha=make_density(c_alpha,n_alpha)
        D_beta=make_density(c_beta,n_beta)
        D=D_alpha+D_beta
        hf_energy=scf_energy_uhf(D,f_dash_alpha,f_dash_beta,D_alpha,D_beta,h1e)
        #println(hf_energy)
        push!(list_e,sum(hf_energy)) 
        #println(list_e)
        for i in list_e
            del_e=list_e[lastindex(list_e)]-list_e[lastindex(list_e)-1]
        end
        #println(del_e)
        if (abs(del_e))<=10^(-12)
            break
        end
        Hartree_fock_energy= sum(hf_energy)+nuclear_repul
        println("iteration= ",n,"    energy= ", Hartree_fock_energy,"         delta_e= ",round(del_e,digits=12))
    end
    #println(c_alpha)
    return Hartree_fock_energy,c_alpha,c_beta
end