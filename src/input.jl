include(".//..//src//attributes.jl")
include(".//..//src/Intgl.jl")
include(".//../src//rhf.jl")
include(".//../src//mp2.jl")
include(".//..//src//basis_build.jl")

#Take the first argument as the input file that contains information about the geometry 
# charge, multiplicity and level of theory
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
println(typeof(Atoms))
geom = []
for line in geom_file
    push!(geom, [parse(Float64, x) for x in line[2:end]])
end
println(geom)
println(typeof(geom))
println(basis_set)
function do_scf(Atoms::Vector{Any},coordinates::Vector{Any},basis_set)
    """
    A function  that does scf calculation
    Attributes:
    Atoms: array/Vector that contains the atoms 
    coordinates : vector/Array that  contains the x,y,z geomtry 
    basis_set: String
    Returns:
        Hartree Fock energy,mo coefficient matrix ,fock matrix,no of basisfunction
    """
    
    #exp,coeff,origin,shells,norms=orbital_config(Atoms,geom,basis_set)
    exp,coeff,origin,shells,norms=basis_fig_out(Atoms,geom,basis_set)
    if Level_of_theory=="hf"
        @time S_matrix=S_mat(exp,coeff,origin,shells,norms)
        @time kinetic_energy=T_mat(exp,coeff,origin,shells,norms)
        no_of_e , atomic_nos = no_of_electrons(Atoms)
        @time Potential_mat = V_mat(exp,coeff,origin,shells,norms,atomic_nos,geom)
        @time E_nuc=enuc(atomic_nos,geom)
        core_h=kinetic_energy+Potential_mat
        @time twoe,eri = Eri_mat(exp,coeff,origin,shells,norms)
        @time Hartree_fock_energy,c,fock,nbasis,eps=scf(S_matrix,kinetic_energy,Potential_mat,twoe,E_nuc,no_of_e)
        println("final Hartree fock energy= ", Hartree_fock_energy)
    end
    return Hartree_fock_energy,c,fock,nbasis,no_of_e,twoe,eps
end
Hartree_fock_energy,c,fock,nbasis,no_of_e,eri,eps =do_scf(Atoms,geom,basis_set)
@time mp2_energy=compute_mp2(c,eri,no_of_e,eps)
print("MP2 correlation Energy = ",mp2_energy ,"\n")
print("Total Energy = ",mp2_energy+Hartree_fock_energy,"\n")