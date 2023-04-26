
using LinearAlgebra
using DelimitedFiles
using LinearAlgebra
#using Base.Iterators: product
#using StatsBase

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
    

function formBasis( ::Molecule )
    """Routine to create the basis from the input molecular geometry and
       basis set. On exit, you should have a basis in self.bfs, which is a 
       list of BasisFunction objects. This routine also defines the center
       of nuclear charge. 
    """
    # initialize the list of Basis objects
    bfs = []
    for atom in molecule.atoms
        for (momentum, prims) in molecule.basis_data[atom.charge]
            exps = [e for (e,c) in prims]
            coefs = [c for (e,c) in prims]
            for shell in momentum2shell(momentum)
                push!(bfs, Basis_int(atom.origin, shell, length(exps), exps, coefs))
            end
        end
    end
    
    molecule.bfs = bfs
    molecule.nbasis = length(bfs)
    
    # create masking vector for geometric derivatives
    idx = 1
    for atom in molecule.atoms
        atom.mask = zeros(molecule.nbasis)
        for (momentum, prims) in molecule.basis_data[atom.charge]
            for shell in momentum2shell(momentum)
                atom.mask[idx] = 1.0
                idx += 1
            end
        end
    end

    # note this is center of positive charge (atoms only, no electrons)
    charges = [atom.charge for atom in molecule.atoms]
    origins = [atom.origin for atom in molecule.atoms]
    center_of_charge = sum(origins .* charges', dims=1) ./ sum(charges)
    molecule.center_of_charge = center_of_charge
end
function build(self::Molecule, direct::Bool=false)
    """Routine to build necessary integrals"""
    one_electron_integrals!(self)
    if direct
        # populate dict for screening
        self.screen = Dict{Int, ERIs}()
        for p in 1:self.nbasis
            for q in 1:p
                pq = p*(p+1)÷2 + q
                self.screen[pq] = ERI(self.bfs[p],self.bfs[q],self.bfs[p],self.bfs[q])
            end
        end
    else
        two_electron_integrals!(self)
    end
    self.is_built = true
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
function no_of_electrons(atoms::Array{String, 1})
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

function getBasis(self, filename)
    """
    Routine to read the basis set files (EMSL Gaussian 94 standard)
    The file is first split into atoms, then iterated through (once).
    At the end we get a basis, which is a dictionary of atoms and their
    basis functions: a tuple of angular momentum and the primitives
    
    Return: Dict{Int, Tuple{String, Tuple{Float64, Float64}[]}[]}
    """
    basis = Dict{Int, Tuple{String, Tuple{Float64, Float64}[]}[]}()

    basisset = open(filename) do file
        read(file, String)
    end
    data = split(basisset, "****")

    # Iterate through all atoms in basis set file
    for i in 2:lastindex(data)
        atomData = [split(x) for x in split(data[i], '\n')[2:end-1]]
        for (idx, line) in enumerate(atomData)
            # Ignore empty lines
            if isempty(line)
                continue
            # first line gives atom
            elseif idx == 1
                @assert length(line) == 2
                atom = sym2num(line[1])
                basis[atom] = []
                # now set up primitives for particular angular momentum
                newPrim = true
            # Perform the set up once per angular momentum
            elseif idx > 1 && newPrim
                momentum  = line[1]
                numPrims  = parse(Int, line[2])
                newPrim   = false
                count     = 0
                prims     = []
                prims2    = [] # need second list for 'SP' case
            else
                # Combine primitives with its angular momentum, add to basis
                if momentum == "SP"
                    # Many basis sets share exponents for S and P basis
                    # functions so unfortunately we have to account for this.
                    push!(prims, (parse(Float64, replace(line[1], "D" => "E")), parse(Float64, replace(line[2], "D" => "E"))))
                    push!(prims2, (parse(Float64, replace(line[1], "D" => "E")), parse(Float64, replace(line[3], "D" => "E"))))
                    count += 1
                    if count == numPrims
                        push!(basis[atom], ("S", prims))
                        push!(basis[atom], ("P", prims2))
                        newPrim = true
                    end
                else
                    push!(prims, (parse(Float64, replace(line[1], "D" => "E")), parse(Float64, replace(line[2], "D" => "E"))))
                    count += 1
                    if count == numPrims
                        push!(basis[atom], (momentum, prims))
                        newPrim = true
                    end
                end
            end
        end
    end

    return basis
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

function read_molecule(geometry::String)
    # atomic masses (isotop avg)
    masses = [1.008,4.003,6.941,9.012,10.812,12.011,14.007,5.999,
        18.998,20.180,22.990,24.305,26.982,28.086,30.974,32.066,
        35.453,39.948]
    f = split(geometry, '\n')
        # remove any empty lines
    f = filter(x -> !isempty(x), f)            
    atomlist = []
    for (line_number, line) in enumerate(f)
        if line_number == 1
            @assert length(split(line)) == 2
            charge = parse(Int, split(line)[1])
            multiplicity = parse(Int, split(line)[2])
            #println(charge)
            #println(multiplicity)
        else 
            if length(split(line)) == 0
                break
            end
            @assert length(split(line)) == 4
            sym = sym2num(split(line)[1])
            name=symbol_name(sym)
            println(name)
            mass = masses[sym]
            println(sym)
            println(mass)
            # Convert Angstrom to Bohr (au)
            x   = parse(Float64, split(line)[2])/0.52917721092
            y   = parse(Float64, split(line)[3])/0.52917721092
            z   = parse(Float64, split(line)[4])/0.52917721092
            #println([x,y,z])
            # Convert amu to atomic units
            mass *= 1822.8885
            atom = my_Atom(sym,mass,name,[x,y,z])
            #println(atom)
            push!(atomlist, atom)
            #println(atomlist)
        end
    end
    #println(atomlist)
    return  0,1,atomlist
end

function one_electron_integrals(obj)
    """
    Routine to set up and compute one-electron integrals
    """
    N = obj.nbasis

    # core integrals
    obj.S = zeros(N,N)
    obj.V = zeros(N,N)
    obj.T = zeros(N,N)

    # dipole integrals
    obj.M = zeros(3,N,N)
    obj.mu = zeros(3,Complex)

    # angular momentum
    obj.L = zeros(3,N,N)

    obj.nuc_energy = 0.0

    # Get one electron integrals
    for i in 1:N
        for j in 1:i
            obj.S[i,j] = obj.S[j,i] = S(obj.bfs[i], obj.bfs[j])
            obj.T[i,j] = obj.T[j,i] = T(obj.bfs[i], obj.bfs[j])

            obj.M[1,i,j] = obj.M[1,j,i] = Mu(obj.bfs[i], obj.bfs[j], obj.center_of_charge, "x")
            obj.M[2,i,j] = obj.M[2,j,i] = Mu(obj.bfs[i], obj.bfs[j], obj.center_of_charge, "y")
            obj.M[3,i,j] = obj.M[3,j,i] = Mu(obj.bfs[i], obj.bfs[j], obj.center_of_charge, "z")

            for atom in obj.atoms
                obj.V[i,j] += -atom.charge * V(obj.bfs[i], obj.bfs[j], atom.origin)
            end
            obj.V[j,i] = obj.V[i,j]

            # RxDel is antisymmetric
            obj.L[1,i,j] = RxDel(obj.bfs[i], obj.bfs[j], obj.center_of_charge, "x")
            obj.L[2,i,j] = RxDel(obj.bfs[i], obj.bfs[j], obj.center_of_charge, "y")
            obj.L[3,i,j] = RxDel(obj.bfs[i], obj.bfs[j], obj.center_of_charge, "z")

            obj.L[:,j,i] = -1 * obj.L[:,i,j]
        end
    end

    # Compute nuclear repulsion energy
    for pair in combinations(obj.atoms,2)
        obj.nuc_energy += pair[1].charge * pair[2].charge / norm(pair[1].origin - pair[2].origin)
    end

    # Preparing for SCF
    obj.Core = obj.T + obj.V
    obj.X = pow(obj.S,-0.5)
    obj.U = pow(obj.S,0.5)
end
function two_electron_integrals(self)
    """
    Routine to setup and compute two-electron integrals
    """
    N = self.nbasis
    self.TwoE = zeros(N,N,N,N)  
    doERIs!(N, self.TwoE, self.bfs)
    self.TwoE = Array(self.TwoE)
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
function getBasis(self,filename)
    basis = Dict()
    f = open(filename, "r")
    data = split(read(f, String),"****")
    close(f)

    for i in 2:lastindex(data)
        atomData = split.(split(data[i],'\n')[2:end-1])
        newPrim = true
        for idx = 1:lastindex(atomData)
            line = atomData[idx]
            if isempty(line)
                continue
            elseif idx == 1
                @assert length(line) == 2
                atom = self.sym2num(line[1])
                basis[atom] = []
            elseif idx > 1 && newPrim
                momentum  = line[1]
                numPrims  = parse(Int64,line[2])
                newPrim   = false
                count     = 0
                prims     = []
                prims2    = []
            else
                if momentum == "SP"
                    push!(prims,(parse(Float64,replace(line[1],'D'=>'E')),parse(Float64,replace(line[2],'D'=>'E'))))
                    push!(prims2,(parse(Float64,replace(line[1],'D'=>'E')),parse(Float64,replace(line[3],'D'=>'E'))))
                    count += 1
                    if count == numPrims
                        push!(basis[atom],("S",prims))
                        push!(basis[atom],("P",prims2))
                        newPrim = true
                    end
                else
                    push!(prims,(parse(Float64,replace(line[1],'D'=>'E')),parse(Float64,replace(line[2],'D'=>'E'))))
                    count += 1
                    if count == numPrims
                        push!(basis[atom],(momentum,prims))
                        newPrim = true
                    end
                end
            end
        end
    end

    return basis
end
