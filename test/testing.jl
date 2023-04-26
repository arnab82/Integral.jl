
include(".//..//src//attributes.jl")
include(".//..//src/Intgl.jl")
basis_set="sto3g"
atoms=["O","H","H"]
geom=[[0.000000000000 , -0.143225816552  , 0.000000000000],[1.638036840407 ,  1.136548822547 , -0.000000000000],[-1.638036840407  , 1.136548822547  ,-0.000000000000]]
exp,coeff,origins,shells,norms=orbital_config(atoms,geom,basis_set)
println(typeof(origins))
origin=ang2bohr.*origins
println("the exponents are",exp,"\n")
display(exp)
println("the coefficients are",coeff,"\n")
display(coeff)
println("the origins are",origin,"\n")
display(origin)
println(typeof(origin))
println("the shells are",shells,"\n")
display(shells)
println("the norms are",norms,"\n")
display(norms)
S_mat(exp,coeff,origin,shells,norms)
T_mat(exp,coeff,origin,shells,norms)
#println(ang2bohr)
"""inFile = ARGS[1]
f = open(inFile, "r")
content = readlines(f)
input_file = []
for line in content
    v_line = split(strip(line))
    if length(v_line) > 0
        push!(input_file, v_line)
    end
end

Level_of_theory = input_file[1][1]
basis_set = input_file[1][2]
charge, multiplicity = input_file[2][1:end]

input_file = input_file[3:end]
geom_file = input_file
Atoms = []
for line in geom_file
    push!(Atoms, line[1])
end

geom_raw = []
for line in geom_file
    push!(geom_raw, [parse(Float64, x) for x in line[2:end]])
end

geom = convert(Array{Float64,2}, geom_raw)

# Convert geom_raw to a Julia array
geom = convert(Array{Float64, 2}, geom_raw)"""

no_of_e , atomic_nos = no_of_electrons(atoms)
println(atomic_nos)

Potential_mat = V_mat(exp,coeff,origin,shells,norms,atomic_nos,geom)
E_nuc=enuc(atomic_nos,geom)
println(E_nuc)
twoe,eri = Eri_mat(exp,coeff,origin,shells,norms)
display(twoe)