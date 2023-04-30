
using TensorOperations
using Einsum

function compute_mp2(c::Matrix{Any},eri::Array{Float64, 4},no_of_e::Int64,eps)
    new_eri=zeros(Float64,nbasis,nbasis,nbasis,nbasis)
    temp1=zeros(Float64,nbasis,nbasis,nbasis,nbasis)
    temp2=zeros(Float64,nbasis,nbasis,nbasis,nbasis)
    temp3=zeros(Float64,nbasis,nbasis,nbasis,nbasis)
    @einsum temp1[p,j,k,l]=(c[i,p]*eri[i,j,k,l])
    @einsum temp2[p,q,k,l]=(c[j,q]*temp1[p,j,k,l])
    @einsum temp3[p,q,r,l]=(c[k,r]*temp2[p,q,k,l])
    @einsum new_eri[p,q,r,s]=(c[l,s]*temp3[p,q,r,l])
    #println(typeof(new_eri))
    ndocc=Int(no_of_e/2)
    emp2=0.0

    for i in 1:ndocc
        for a in (ndocc+1):nbasis
            for j in 1:ndocc
                for b in (ndocc+1):nbasis
                    emp2+=((*(new_eri[i,a,j,b],(*(2,new_eri[i,a,j,b])-new_eri[i,b,j,a])))/(eps[i]+eps[j]-eps[a]-eps[b]))
                end
            end
        end
    end
    return emp2
end

