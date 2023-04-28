
using JSON3
function momentum2shell(momentum::Int64)
    """Routine to convert angular momentum to Cartesian shell pair in order
       to create the appropriate BasisFunction object (e.g. form px,py,pz)
    """
    shells = Dict(
        0 => [(0,0,0)],
        1 => [(1,0,0),(0,1,0),(0,0,1)],
        2 => [(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)],
        3 => [(3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,0,2),
               (0,3,0),(0,2,1),(0,1,2), (0,0,3)]
    )
    return shells[momentum]
end
# Read the contents of the file into a string
json_string = read("/home/arnabbachhar/Intgl/basis/3-21g.json", String)

# Parse the JSON string into a Julia data structure
bs = JSON3.read(json_string)
expon=[]
coeffic=[]
angu_mom=[]
attributes=[[],[],[],[],[]]
for i in 1:lastindex(bs.elements["8"].electron_shells)
    exps = bs.elements["8"].electron_shells[i]["exponents"]
    coeff = bs.elements["8"].electron_shells[i]["coefficients"]
    angm = bs.elements["8"].electron_shells[i]["angular_momentum"]
    println(length(angm))
    
    push!(expon,exps)
    push!(coeffic,coeff)
end
println(expon)
println(coeffic)
    for k in 1:lastindex(angu_mom)
        if angu_mom[k]==[(0,0,0)]
            push!(attributes[k],expon[k])
            push!(attributes[k],coeffic[k][1])
            push!(attributes[k],angu_mom[k])
        elseif angu_mom[k]==[(1,0,0)] || angu_mom[k]==[(0,1,0)] ||angu_mom[k]==[(0,0,1)] 
            push!(attributes[k],expon[2])
            push!(attributes[k],coeffic[2][2])
            push!(attributes[k],angu_mom[k])
        elseif 6<k<20
            println("thinking")
        end
    end


    for j in 1:lastindex(angm)
        if angm[j]==0
            angum=momentum2shell(angm[j])
            push!(angu_mom,angum)
        elseif angm[j]>0
            angum=momentum2shell(angm[j])
            for l in 1:lastindex(angum)
                push!(angu_mom,[angum[l]])
            end
        else
            println("angular_momentum count exceeded")

        end

    end

#println("the attributes are ",attributes)
#display(attributes)

println(size(angu_mom))
println(angu_mom)
nested_object_property = bs.elements["8"].electron_shells

#a=momentum2shell("1")[1]