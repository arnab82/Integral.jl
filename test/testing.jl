using PyCall
include(".//..//src//attributes.jl")
include(".//..//src/Intgl.jl")
include(".//../src//rhf.jl")
basis_set="sto3g"
atoms=["O","H","H"]
geom=[[0.000000000000 , -0.143225816552  , 0.000000000000],[1.638036840407 ,  1.136548822547 , -0.000000000000],[-1.638036840407  , 1.136548822547  ,-0.000000000000]]
exp,coeff,origin,shells,norms=orbital_config(atoms,geom,basis_set)
#println(typeof(origins))
#origin=ang2bohr.*origins
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
S_matrix=S_mat(exp,coeff,origin,shells,norms)
println(S_mat(exp,coeff,origin,shells,norms))

#display(S_mat(exp,coeff,origin,shells,norms))
kinetic_energy=T_mat(exp,coeff,origin,shells,norms)
#display(T_mat(exp,coeff,origin,shells,norms))
#println(ang2bohr)



pyscf=PyCall.pyimport("pyscf")
function make_molecule()
    mol=pyscf.gto.M()
    atoms="O 0.000000000000 -0.143225816552 0.000000000000;H 1.638036840407 1.136548822547 -0.000000000000;H -1.638036840407 1.136548822547 -0.000000000000"
    mol.charge=0
    mol.unit = "Bohr"
    mol.spin=0
    mol.build(
	    atom= atoms,
	    basis = "sto3g")
	    #basis = "cc-pVDZ")
    return mol
end
mol=make_molecule()

function pyscf_2e(mol)
    eri = mol.intor("int2e")
    return eri
end
#display(pyscf_2e(mol))

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
println(typeof(Potential_mat))
E_nuc=enuc(atomic_nos,geom)
display(Potential_mat)
println(E_nuc)
core_h=kinetic_energy+Potential_mat
display(core_h)
println(core_h)
twoe,eri = Eri_mat(exp,coeff,origin,shells,norms)
println(typeof(twoe))
"""[array([27.55116782,  7.68181999,  2.88241787]), array([2.39491488, 0.80156184, 0.34520813]), array([10.74583263,  1.73374407,  0.42581893]), array([10.74583263,  1.73374407,  0.42581893]), array([10.74583263,  1.73374407,  0.42581893]), array([1.79444183, 0.50032649, 0.18773546]), array([1.79444183, 0.50032649, 0.18773546])]"""

s_half=make_s_half(S_matrix)
display(s_half)
h,nbasis=size(S_matrix)
h1e=kinetic_energy+Potential_mat
fock= zeros(Float64,nbasis,nbasis)
init_fock= (h1e*make_s_half(S_matrix)')
init_fock=(make_s_half(S_matrix)*init_fock)
init_fock=(init_fock+init_fock')/2
display(init_fock)
e,c0 = LinearAlgebra.eigen(init_fock)
c = (make_s_half(S_matrix)*c0)
display(c)
D=make_density(c,no)
display(D)
initial_E=scf_energy(D,init_fock,core_h)
println(initial_E)
#Hartree_fock_energy,E,c,fock,nbasis=scf(S_matrix,kinetic_energy,Potential_mat,twoe,E_nuc)
#println("final Hartree fock energy= ", Hartree_fock_energy)
#println(E)
#println("the value of coefficient matrix is",c)
#println("the size of fock is ",size(fock))
#endtime=time.time()
#println("the runtime of the code is", endtime-start)
display(uniqueindex(7))
display(twoe)
println(size(twoe))