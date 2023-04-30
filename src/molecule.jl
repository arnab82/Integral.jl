using PyCall
using LinearAlgebra
using DelimitedFiles
using LinearAlgebra
#using Base.Iterators: product
#using StatsBase
sp=pyimport("scipy.special")
hyp1f1=sp.hyp1f1
fact2=sp.factorial2
mutable struct my_Atom
    mass_number::Int64
    mass::Float64
    symbol::String
    xyz::Array{Float64,1}
end
mutable struct BasisFunction
    """
    A struct that contains all our basis function data
    Attributes:
    origin: array/list containing the coordinates of the Gaussian origin
    shell: tuple of angular momentum
    exps: list of primitive Gaussian exponents
    coefs: list of primitive Gaussian coefficients
    num_exps: number of primitive Gaussian exponents
    norm: list of normalization factors for Gaussian primitives
    """
    origin::Array{Float64}
    shell::Array{Float64}
    exps::Array{Float64}
    coefs::Array{Float64}
    num_exps::Int64
    norms::Array{Float64}
end

function BasisFunction(; origin=[0.0, 0.0, 0.0], shell=(0, 0, 0), exps=[], coefs=[], num_exps=0)
    if num_exps == 0
        num_exps = length(exps)
    end
    bf = BasisFunction(origin, shell, exps, coefs, num_exps, [0.0])
    normalize_basis!(bf)
    return bf
end

function normalize_basis!(bf::BasisFunction)
    """
    Routine to normalize the basis functions, in case they do not integrate to unity.
    """
    l, m, n = bf.shell
    l=convert(Int64,l)
    m=convert(Int64,m)
    n=convert(Int64,n)
    println(l,m,n)
    L = l + m + n
    #println(L)
    #println(bf.exps)
    #println(fact2(2*l-1))
    #println((2.0^(2*(l + m + n) + 1.5)).*bf.exps.^((l + m + n + 1.5) ))

    #println((fact2(2*l - 1) .* fact2(2*m - 1) .* fact2(2*n - 1) .* pi^1.5))
    # norm is an array of length equal to number primitives
    # normalize primitives first (PGBFs)
    bf.norms = sqrt.(2.0^(2*(l + m + n) + 1.5) .* bf.exps.^(l + m + n + 1.5) /(fact2(2*l - 1) .* fact2(2*m - 1) .* fact2(2*n - 1) .* pi^1.5))

    # now normalize the contracted basis functions (CGBFs)
    # Eq. 1.44 of Valeev integral whitepaper
    prefactor = (pi^1.5 .* fact2(2*l - 1) .* fact2(2*m - 1) .* fact2(2*n - 1)) / 2.0^L

    N = 0.0
    for ia in 1:bf.num_exps
        for ib in 1:bf.num_exps
            N += bf.norms[ia] * bf.norms[ib] * bf.coefs[ia] * bf.coefs[ib] /
                 (bf.exps[ia] + bf.exps[ib])^(L + 1.5)
        end
    end

    N *= prefactor
    N = 1 / sqrt(N)
    bf.coefs .*= N
    """for ia in 1:num_exps
        bf.coefs[ia]*= N
    end"""
    return N,bf.norms
end

struct Basis_int
    origin::Array{Float64,1}
    momentum::Array{Int64,1}
    nprims::Int64
    exps::Array{Float64,1}
    coefs::Array{Float64,1}
end

struct Molecule
    """
    A struct that contains type of molecule
    Attributes:
    atoms: array/list of struct type my_Atom
    charge: contains the total charge of the moecule
    multiplicity: contains the multiplicity of the molecule => 2S+1 ;S= spin
    bfs ::Array of basisfunction
    """
    atoms::Array{my_Atom,1}
    charge::Int64
    multiplicity::Int64
    nelec::Int64
    nocc::Int64
    basis_data::Dict{Int64,Array{Tuple{Tuple{Char,Int64,Int64},Array{Tuple{Float64,Float64},1}},1}}
    bfs::Array{Basis_int,1}
    nbasis::Int64
    center_of_charge::Array{Float64,1}
    is_built::Bool
    geometry_input::String

        
    function Molecule(geometry::String, basis::String="sto-3g")
        # geometry is now specified in input file
        charge, multiplicity, atomlist = read_molecule(geometry)
        println(charge)
        println(atomlist)
        new(charge, multiplicity, atomlist, sum([atom.charge for atom in atomlist])-charge, 
            (sum([atom.charge for atom in atomlist])-charge)÷2, false, geometry, getBasis(joinpath(dirname(@__FILE__),lowercase(basis)*".gbs")), nothing)
        formBasis!(this)
    end
end
    


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
function symbol_name(sym)
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
    return symbol[sym]
end
function no_of_electrons(atoms::Array{Any,1})
    symbol = ["H","He", "Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K", "Ca", "Sc",
     "Ti", "V", "Cr", "Mn", "Fe","Co", "Ni", "Cu", "Zn","Ga", "Ge", "As", "Se", "Br", "Kr","Rb", "Sr", "Y", "Zr",
     "Nb", "Mo", "Tc", "Ru","Rh", "Pd", "Ag", "Cd","In", "Sn", "Sb", "Te", "I", "Xe","Cs", "Ba", "La", "Ce", "Pr",
     "Nd", "Pm", "Sm",  "Eu","Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu","Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
     "Au", "Hg","Tl","Pb","Bi","Po","At","Rn"]
    
    mole_elec = 0
    atomic_nos = []
    for i in 1:length(atoms)
        mole_elec += findfirst(x -> x == atoms[i], symbol)[1]
        push!(atomic_nos, findfirst(x -> x == atoms[i], symbol)[1])
    end
    
    return mole_elec, atomic_nos
end
function no_of_electrons(atoms::Vector{String})
    symbol = ["H","He", "Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K", "Ca", "Sc",
     "Ti", "V", "Cr", "Mn", "Fe","Co", "Ni", "Cu", "Zn","Ga", "Ge", "As", "Se", "Br", "Kr","Rb", "Sr", "Y", "Zr",
     "Nb", "Mo", "Tc", "Ru","Rh", "Pd", "Ag", "Cd","In", "Sn", "Sb", "Te", "I", "Xe","Cs", "Ba", "La", "Ce", "Pr",
     "Nd", "Pm", "Sm",  "Eu","Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu","Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
     "Au", "Hg","Tl","Pb","Bi","Po","At","Rn"]
    
    mole_elec = 0
    atomic_nos = []
    for i in 1:length(atoms)
        mole_elec += findfirst(x -> x == atoms[i], symbol)[1]
        push!(atomic_nos, findfirst(x -> x == atoms[i], symbol)[1])
    end
    
    return mole_elec, atomic_nos
end


function momentum2shell(momentum::String)
    """Routine to convert angular momentum to Cartesian shell pair in order
       to create the appropriate BasisFunction object (e.g. form px,py,pz)
    """
    shells = Dict(
        "S" => [(0,0,0)],
        "P" => [(1,0,0),(0,1,0),(0,0,1)],
        "D" => [(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)],
        "F" => [(3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,0,2),
               (0,3,0),(0,2,1),(0,1,2), (0,0,3)]
    )
    return shells[momentum]
end




function save_integrals(self::Molecule, folder::String)
    """Routine to save integrals for SCF in Crawford group format"""
    if folder === nothing
        error("Please provide a folder to save the integrals.")
    else
        if !self.is_built
            self.build()
        end
        mkpath(folder) # careful! will overwrite.

        writedlm(joinpath(folder, "enuc.dat"), reshape([self.nuc_energy], 1))
        writedlm(joinpath(folder, "nbf.dat"), reshape([self.nbasis], 1))
        writedlm(joinpath(folder, "nelec.dat"), reshape(self.nelec, 1))
        writedlm(joinpath(folder, "s.dat"), self.S)
        writedlm(joinpath(folder, "t.dat"), self.T)
        writedlm(joinpath(folder, "v.dat"), self.V)

        open(joinpath(folder, "eri.dat"), "w") do f
            for (i, j, k, l) in product(0:self.nbasis-1, 0:self.nbasis-1, 0:self.nbasis-1, 0:self.nbasis-1)
                println(f, i+1, " ", j+1, " ", k+1, " ", l+1, " ", self.TwoE[i+1,j+1,k+1,l+1])
            end
        end

        open(joinpath(folder, "geometry.txt"), "w") do f
            println(f, self.geometry_input)
        end
    end
end
