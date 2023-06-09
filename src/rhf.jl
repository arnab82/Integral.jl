using PyCall
time=PyCall.pyimport("time")
start=time.time()
using LinearAlgebra,Statistics

n_elec=10
no=Int64((n_elec)/2)
function scf_energy(D::Matrix{Float64},fock::Matrix{Float64},h1e::Matrix{Float64})
    """
    Attributes:
        D: Density matrix 
        fock: Fock matrix
        h1e: core hamiltonian matrix
    return:
        scf energy: Float64 
    """
    h=[]
    h,nbasis=size(D)
    n_hf_energy=0.0
    for i in 1:nbasis
        for j in 1:nbasis
            n_hf_energy+=(D[i,j]*(fock[i,j]+h1e[i,j]))
        end
    end
    return n_hf_energy
end

function make_fock(D::Matrix{Float64},h1e::Matrix{Float64},eri::Array{Float64, 4})
    """
    Attributes:
        D: Density matrix 
        h1e: core hamiltonian matrix
        eri: two electron matrix of size(nbasis,nbasis,nbasis,nbasis)
    return:
        Fock matrix
    """
    h=[]
    h,nbasis=size(D)
    fock= zeros(Float64,nbasis,nbasis)
    for i in 1:nbasis
        for j in 1:nbasis
            for k in 1:nbasis
                for l in 1:nbasis
                    fock[i,j]+=(D[k,l]*((2*eri[i,j,k,l])-eri[i,k,j,l]))
                #push!(new_fock[i,j],(*(D[k,l],(*((2*twoe[i,j,k,l])-twoe[i,l,k,j])))))
                end
            end
        end
    end 
    return fock+h1e
end
function make_density(c::Matrix{Any},no::Int64)
    """
    Attributes:
        C: Mo co-efficient matrix 
        no: no of occupied orbitals
    return:
        Density matrix
    """
    h=[]
    h,nbasis=size(c)
    D=zeros(Float64,nbasis,nbasis)
    for i in 1:nbasis
        for j in 1:nbasis
            for m in 1:no
                D[i,j]+=(c[i,m]*c[j,m])
                #push!(D,(*(c[i,m],c[j,m])))
            end
        end
    end   
    return D
end

#println("the one electron integral is ",h1e ,"\n")
#println("the two electron integral is ",eri,"\n")
#println(" shape of two electron integral is",size(eri),"\n")
#println("the nuclear nuclear repulsion term is ",constant,"\n")
#println("the overlap integral is",S)
function make_s_half(S::Matrix{Float64}) 
    """
    Attributes:
        S: overlap matrix
    return:
        S^(-1/2)
    """
    s=(S+S')/2
    #println(nbasis)
    q,L=LinearAlgebra.eigen(s)
    q_half=[]
    for i in 1:lastindex(q)
        push!(q_half,q[i]^(-0.5))
    end
    q_half=Diagonal(q_half)
    s_half=(L*q_half)
    s_half=(s_half*L')
    #println("the value of s_half is    ",s_half)
    return s_half
end 
#scf initialisation 

function scf(s::Matrix{Float64},T_mat::Matrix{Float64},V_mat::Matrix{Float64},Eri::Array{Float64, 4},nuclear_repul::Float64,no_of_e::Int64)
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
    """
    h,nbasis=size(s)
    h1e=T_mat+V_mat
    fock= zeros(Float64,nbasis,nbasis)
    s_half=make_s_half(s)
    init_fock= (h1e*s_half')
    init_fock=(s_half*init_fock)
    init_fock=(init_fock+init_fock')/2
    Hartree_fock_energy=0.0
    del_e=0.0
    no=Int(no_of_e/2)
    e,c0 = LinearAlgebra.eigen(init_fock)
    c = (s_half*c0)
    list_e=[0.0,]
    eps_list=[]
    fock= zeros(Float64,nbasis,nbasis)
    for n in 1:100
        if n==1
            D=zeros(Float64,nbasis,nbasis)
            for i in 1:no
                D[i,i]+=1.0
            end
        end
        D=make_density(c,no)
        fock=make_fock(D,h1e,Eri)
        f_dash= (fock*make_s_half(s)')
        f_dash=(make_s_half(s)*f_dash)
        f_dash=(f_dash+transpose(f_dash))/2
        eps,c_dash=LinearAlgebra.eigen(f_dash)
        push!(eps_list,eps)
        c=(make_s_half(s)*c_dash)
        D=make_density(c,no)
        hf_energy=scf_energy(D,fock,h1e)
        push!(list_e,sum(hf_energy)) 
        #println(list_e)
        for i in list_e
            del_e=list_e[lastindex(list_e)]-list_e[lastindex(list_e)-1]
        end
        #println(del_e)
        if (abs(del_e))<=10^(-12)
            break
        end
        Hartree_fock_energy=sum(hf_energy)+nuclear_repul
        println("iteration= ",n,"    energy= ", sum(hf_energy)+nuclear_repul,"         delta_e= ",round(del_e,digits=12))
    end
    return Hartree_fock_energy,c,fock,nbasis,eps_list[lastindex(eps_list)]
end
