using JSON
using JSON3
include("./molecule.jl")
include("./attributes.jl")
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
"""file_name= read("/home/arnabbachhar/Intgl/basis/sto-3g.json",String)
bs = JSON.parse(file_name)
expon=[]
coeffic=[]
angu_mom=[]
attributes=[]
for i in 1:lastindex(bs["elements"]["31"]["electron_shells"])
    exps = bs["elements"]["31"]["electron_shells"][i]["exponents"]
    coeff = bs["elements"]["31"]["electron_shells"][i]["coefficients"]
    angm = bs["elements"]["31"]["electron_shells"][i]["angular_momentum"]
    push!(angu_mom,angm)
    push!(expon,exps)
    push!(coeffic,coeff)
end
println(expon)
println(coeffic)
println(angu_mom)
for k in 1:lastindex(angu_mom)
    if angu_mom[k]==[0]
            push!(attributes,[expon[k],coeffic[k][1],[0,0,0]])
    elseif angu_mom[k]==[0,1]
            push!(attributes,[expon[k],coeffic[k][1],[0,0,0]])
            push!(attributes,[expon[k],coeffic[k][2],[1,0,0]])
            push!(attributes,[expon[k],coeffic[k][2],[0,1,0]])
            push!(attributes,[expon[k],coeffic[k][2],[0,0,1]])
    elseif angu_mom[k]==[2]
        push!(attributes,[expon[k],coeffic[k][1],[2,0,0]])
        push!(attributes,[expon[k],coeffic[k][1],[1,1,0]])
        push!(attributes,[expon[k],coeffic[k][1],[1,0,1]])
        push!(attributes,[expon[k],coeffic[k][1],[0,2,0]])
        push!(attributes,[expon[k],coeffic[k][1],[0,1,1]])
        push!(attributes,[expon[k],coeffic[k][1],[0,0,2]])

    else
        println("this is not supoorted")
    end
end
println(attributes)
display(attributes)"""
function sym2num(sym)
    """Routine that converts atomic symbol to atomic number"""
    symbol = [
        "H","He",
        "Li","Be","B","C","N","O","F","Ne",
        "Na","Mg","Al","Si","P","S","Cl","Ar",
        "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
        "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
        "Rh", "Pd", "Ag", "Cd",
        "In", "Sn", "Sb", "Te", "I", "Xe",
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",  "Eu",
        "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl","Pb","Bi","Po","At","Rn"]
    return findfirst(x -> x == string(sym), symbol)
end
#expon=[]
#coeffic=[]
#angu_mom=[]
#atoms=["O","H" ,"H"]
#geom=[[0.000000000000 , -0.143225816552  , 0.000000000000],[1.638036840407 ,  1.136548822547 , -0.000000000000],[-1.638036840407  , 1.136548822547  ,-0.000000000000]]
#basis_set="sto-3g"
#println(angu_mom)
#println(expon)
#println(coeffic)
function basis_fig_out(atoms, geom, basis_set)
    attributes = []
    expon=[]
    coeffic=[]
    angu_mom=[]
    coord=[]
    if basis_set=="sto-3g"
        file_name= read("/home/arnabbachhar/Intgl/basis/sto-3g.json",String)
        bs = JSON.parse(file_name)
    elseif basis_set=="sto-6g"
        file2= read("/home/arnabbachhar/Intgl/basis/sto-6g.json",String)
        bs = JSON.parse(file2)
    elseif basis_set=="3-21g"
        file_name= read("/home/arnabbachhar/Intgl/basis/3-21g.json",String)
        bs = JSON.parse(file_name)
    elseif basis_set=="6-31g"
        file_name= read("/home/arnabbachhar/Intgl/basis/6-31g.json",String)
        bs = JSON.parse(file_name)
    elseif basis_set=="6-31g"
        file_name= read("/home/arnabbachhar/Intgl/basis/6-31g.json",String)
        bs = JSON.parse(file_name)
    elseif basis_set=="6-31g(d,p)"
        file_name= read("/home/arnabbachhar/Intgl/basis/6-31g(d,p).json",String)
        bs = JSON.parse(file_name)
    elseif basis_set=="6-31g*"
        file_name= read("/home/arnabbachhar/Intgl/basis/6-31g_st.json",String)
        bs = JSON.parse(file_name)
    elseif basis_set=="6-31g**"
        file_name= read("/home/arnabbachhar/Intgl/basis/6-31g_st__st.json",String)
        bs = JSON.parse(file_name)
    elseif basis_set=="3-21g-uc"
        file_name= read("/home/arnabbachhar/Intgl/basis/3-21g-uc.json",String)
        bs = JSON.parse(file_name)
    else
        println("basis set not present")
    end
	#println(lastindex(atoms))
	for (j ,c) in enumerate(atoms)
        symbol_num=sym2num(c)
        symbol_num=string(symbol_num)
        #println(symbol_num)
        #println(typeof(symbol_num))
        for i in 1:lastindex(bs["elements"][symbol_num]["electron_shells"])
            exps = bs["elements"][symbol_num]["electron_shells"][i]["exponents"]
            coeff = bs["elements"][symbol_num]["electron_shells"][i]["coefficients"]
            angm = bs["elements"][symbol_num]["electron_shells"][i]["angular_momentum"]
            push!(angu_mom,[angm,geom[j]])
            push!(expon,exps)
            push!(coeffic,coeff)
        end
    end
	#println(expon)
    #println(coeffic)
    #println(angu_mom)
    for k in 1:lastindex(angu_mom)
        if angu_mom[k][1]==[0]
            push!(attributes,[expon[k],coeffic[k][1],[0,0,0],angu_mom[k][2]])
        elseif angu_mom[k][1]==[0,1]
            push!(attributes,[expon[k],coeffic[k][1],[0,0,0],angu_mom[k][2]])
            push!(attributes,[expon[k],coeffic[k][2],[1,0,0],angu_mom[k][2]])
            push!(attributes,[expon[k],coeffic[k][2],[0,1,0],angu_mom[k][2]])
            push!(attributes,[expon[k],coeffic[k][2],[0,0,1],angu_mom[k][2]])
        elseif angu_mom[k][1]==[2]
            push!(attributes,[expon[k],coeffic[k][1],[2,0,0],angu_mom[k][2]])
            push!(attributes,[expon[k],coeffic[k][1],[1,1,0],angu_mom[k][2]])
            push!(attributes,[expon[k],coeffic[k][1],[1,0,1],angu_mom[k][2]])
            push!(attributes,[expon[k],coeffic[k][1],[0,2,0],angu_mom[k][2]])
            push!(attributes,[expon[k],coeffic[k][1],[0,1,1],angu_mom[k][2]])
            push!(attributes,[expon[k],coeffic[k][1],[0,0,2],angu_mom[k][2]])

        else
            println("this is not supoorted")
        end
    end
    
    #display(attributes)
    #println(attributes)
    #println(typeof(attributes))
	#println(attributes)
	#println(size(attributes[2][1]))
	orbital_objects = [] 
	for i in 1:lastindex(attributes)
        exp = [parse(Float64, s) for s in attributes[i][1]]
        coeff = [parse(Float64, s) for s in attributes[i][2]]
		norms=BasisFunction(attributes[i][4], attributes[i][3],exp , coeff,size(attributes[i][1])[1],[0.0])
		N,norms=normalize_basis!(norms)
		push!(orbital_objects, BasisFunction(attributes[i][4], attributes[i][3],exp,coeff,size(attributes[i][1])[1],norms))
	end
    println("*****************************************************************************************")
	display(orbital_objects)
	return sort_attributes(orbital_objects)
end
#@time a=basis_fig_out(atoms, geom, basis_set)