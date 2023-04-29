using PyCall
include(".//..//src//attributes.jl")
include(".//..//src/Intgl.jl")
include(".//../src//rhf.jl")
basis_set="sto3g"
atoms=["O","H","H"]
geom=[[0.000000000000 , -0.143225816552  , 0.000000000000],[1.638036840407 ,  1.136548822547 , -0.000000000000],[-1.638036840407  , 1.136548822547  ,-0.000000000000]]

println(typeof(geom))
println(typeof(atoms))
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
#-----------------------------------------
#       Overlap integrals
#-----------------------------------------

@time S_matrix=S_mat(exp,coeff,origin,shells,norms)
println(S_mat(exp,coeff,origin,shells,norms))
#display(S_mat(exp,coeff,origin,shells,norms))


	#----------------------------------------
	#     Kinetic integrals
	#----------------------------------------


@time kinetic_energy=T_mat(exp,coeff,origin,shells,norms)
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
    eri = mol.intor("int2e",aosym="s8")
    return eri
end


no_of_e , atomic_nos = no_of_electrons(atoms)
println(atomic_nos)

	#---------------------------------------
	#   Coulomb integrals
	#---------------------------------------
@time Potential_mat = V_mat(exp,coeff,origin,shells,norms,atomic_nos,geom)
println(typeof(Potential_mat))
display(Potential_mat)
#---------------------------------------
#     ENUC
#--------------------------------------
@time E_nuc=enuc(atomic_nos,geom)
println(E_nuc)

core_h=kinetic_energy+Potential_mat
display(core_h)
println(core_h)
	#----------------------------------------
	#  Two elec integrals 
	#----------------------------------------
@time twoe,eri = Eri_mat(exp,coeff,origin,shells,norms)
@time eri2=doERIs(exp,coeff,origin,shells,norms)
println(eri)

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
@time Hartree_fock_energy,c,fock,nbasis=scf(S_matrix,kinetic_energy,Potential_mat,twoe,E_nuc)
println("final Hartree fock energy= ", Hartree_fock_energy)
@time a=pyscf_2e(mol)
#display(doERIs(exp,coeff,origin,shells,norms))