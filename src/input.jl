include(".//..//src//attributes.jl")
include(".//..//src/Intgl.jl")
include(".//../src//rhf.jl")
inFile = ARGS[1]
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
println(Atoms)
geom = []
for line in geom_file
    push!(geom, [parse(Float64, x) for x in line[2:end]])
end
println(geom)
println(basis_set)
exp,coeff,origin,shells,norms=orbital_config(Atoms,geom,basis_set)
if Level_of_theory=="hf"
    @time S_matrix=S_mat(exp,coeff,origin,shells,norms)
    @time kinetic_energy=T_mat(exp,coeff,origin,shells,norms)
    no_of_e , atomic_nos = no_of_electrons(Atoms)
    @time Potential_mat = V_mat(exp,coeff,origin,shells,norms,atomic_nos,geom)
    @time E_nuc=enuc(atomic_nos,geom)
    core_h=kinetic_energy+Potential_mat
    @time twoe,eri = Eri_mat(exp,coeff,origin,shells,norms)
    @time Hartree_fock_energy,c,fock,nbasis=scf(S_matrix,kinetic_energy,Potential_mat,twoe,E_nuc)
    println("final Hartree fock energy= ", Hartree_fock_energy)
end