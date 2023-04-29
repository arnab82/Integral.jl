
include(".//..//src//attributes.jl")
include(".//..//src/Intgl.jl")
include(".//../src//rhf.jl")

function do_scf(Atoms::Vector{String},coordinates::Vector{Vector{Float64}},basis_set::String)
    exp,coeff,origin,shells,norms=orbital_config(Atoms,coordinates,basis_set)

    @time S_matrix=S_mat(exp,coeff,origin,shells,norms)
    @time kinetic_energy=T_mat(exp,coeff,origin,shells,norms)
    no_of_e, atomic_nos = no_of_electrons(Atoms)
    @time Potential_mat = V_mat(exp,coeff,origin,shells,norms,atomic_nos,geom)
    @time E_nuc=enuc(atomic_nos,geom)
    core_h=kinetic_energy+Potential_mat
    @time twoe,eri = Eri_mat(exp,coeff,origin,shells,norms)
    @time Hartree_fock_energy,c,fock,nbasis=scf(S_matrix,kinetic_energy,Potential_mat,twoe,E_nuc,no_of_e)
    println("final Hartree fock energy= ", Hartree_fock_energy)
    return Hartree_fock_energy,c,fock,nbasis
end


function bond_distance(coords::Vector{Vector{Float64}}, atom1::Int64, atom2::Int64)
    r = sqrt(sum((coords[atom1] .- coords[atom2]).^2))
    return r
end
geom=[[0.0, 0.0, 0.000000000000],[0.745  , 0.0  ,0.000000000000]]
basis_set="sto3g"
atoms=["H","H"]
function get_geometry(r::Float64, Atoms::Vector{String})
    atom1=Atoms[1]
    atom2=Atoms[2]
    # Calculate coordinates of first atom
    x1 = 0.0
    y1 = 0.0
    z1 = 0.0
    
    # Calculate coordinates of second atom
    x2 = r
    y2 = 0.0
    z2 = 0.0
    
    # Return coordinates as a vector of vectors
    return [[x1, y1, z1], [x2, y2, z2]]
end
E_hartree,C,F,nbasis= do_scf(atoms,geom,basis_set)


function hartree_fock_optimization(Atoms::Vector{String},bond_distance::Float64,basis_set::String)
    coordinates=get_geometry(bond_distance,Atoms)
    println(coordinates)
    # Initialize variables
    δ = 0.001
    converged = false
    coordinate_D_Pos=get_geometry(bond_distance+δ,Atoms)
    coordinate_D_neg=get_geometry(bond_distance-δ,Atoms)
    # Perform Hartree-Fock calculation at starting distance R
    E_hartree,C,F,nbasis= do_scf(Atoms,coordinates,basis_set)
    println(E)
    E_δ_pos,C,F,nbasis= do_scf(Atoms,coordinate_D_Pos,basis_set)
    E_δ_neg,C,F,nbasis= do_scf(Atoms,coordinate_D_neg,basis_set)
    # Calculate the first and second derivatives of the energy
    dE_dr = (E_δ_pos-E_δ_neg) / (2 * δ)
    d2E_dr2 = (E_δ_pos+E_δ_neg - 2 * E_hartree) / δ^2
    if abs(dE_dr)>1e-6
        R =bond_distance-dE_dr / d2E_dr2
    end
    for i in 1:100
        println("********************************************************")
        println("iteration no =",i)
        println("********************************************************")

        coordinates=get_geometry(R,Atoms)

        # Perform Hartree-Fock calculation at starting distance R
        E_new,C,F,nbasis= do_scf(Atoms,coordinates,basis_set)

        # Calculate the first and second derivatives of the energy

        coordinate_D_Pos=get_geometry(R+δ,Atoms)
        coordinate_D_neg=get_geometry(R-δ,Atoms)
        E_δ_pos,C,F,nbasis= do_scf(Atoms,coordinate_D_Pos,basis_set)
        E_δ_neg,C,F,nbasis= do_scf(Atoms,coordinate_D_neg,basis_set)
        dE_dr = (E_δ_pos-E_δ_neg) / (2 * δ)
        d2E_dr2 = (E_δ_pos+E_δ_neg - (2 * E_new)) / δ^2
        # Take Newton-Raphson step
        if abs(dE_dr)>1e-6
            R =R-dE_dr / d2E_dr2
        end
        # Check if first derivative is small enough
        println("the value of dE_Dr is ",abs(dE_dr))
        println("Bond distance: ", R, " a.u.")
        if abs(dE_dr) < 1e-6
            converged = true
            break
        end
    end
    coordinates=get_geometry(R,Atoms)

    # Perform Hartree-Fock calculation at starting distance R
    E_new,C,F,nbasis= do_scf(Atoms,coordinates,basis_set)
    return d2E_dr2,R,E_new
end
d2E_dr2,bond_length,E_hartree=hartree_fock_optimization(atoms,0.9,basis_set)
println(d2E_dr2)
function harmonic_frequency(Atoms::Vector{String},bond_distance::Float64,basis_set::String, M1::Float64, M2::Float64)
    # Calculate the second derivative of the energy with respect to bond-distance
    d2E_dr2,bond_length,E_hartree=hartree_fock_optimization(Atoms,bond_distance,basis_set)

    # Compute the reduced mass
    μ = M1 * M2 / (M1 + M2)
    # Compute the harmonic force constant
    k = d2E_dr2
    # Convert the harmonic force constant to units of N/m
    fc = k * 4.35974417e-18 / 1.66053904e-27 * 1e10^2

    # Convert the Hessian energy from atomic units to J/m^2
    Hessian_Jm2 = d2E_dr2* (4.35974417e-18 / 1.66053904e-27 * 1e10^2)^2

    # Convert the Hessian energy from J/m^2 to cm^-1
    Hessian_cm1 = Hessian_Jm2 / (2 * pi * 3e10 * 100)

    # Calculate the vibrational frequency in units of cm^-1
    freq = sqrt(fc/μ )/ (2 * pi * 3e10 )
 
    
    # Report results
    println("Harmonic force constant: ", fc, " N/m")
    println("Harmonic vibrational frequency: ", freq, " cm^-1")
    return k,freq
end
force_constant, frequency=harmonic_frequency(atoms,0.9,basis_set,1.008,1.008)