using JSON3
#include("./../src/attributes.jl")
include("./../src/molecule.jl")
include("./../src/Intgl.jl")
# Read the contents of the file into a string
json_string = read("./../basis/sto-3g.json", String)

# Parse the JSON string into a Julia data structure
bs = JSON3.read(json_string)
nested_object_property = bs.elements["1"].electron_shells[1]


#function getBasis(self, basis_name)
    """
    Routine to read the basis set files (EMSL Gaussian 94 standard)
    The file is first split into atoms, then iterated through (once).
    At the end we get a basis, which is a dictionary of atoms and their
    basis functions: a tuple of angular momentum and the primitives
    Return: {atom: [('angmom',[(exp,coef),...]), ('angmom',[(exp,...}
    Return: Dict{Int, Tuple{String, Tuple{Float64, Float64}[]}[]}
    """

basis_name="sto-3g"
filename=(lowercase(basis_name))*".gbs"
basisset =open("./basis/sto-3g.gbs") do file
    read(file, String)
end

using DelimitedFiles

# Read input file
inFile = ARGS[1]
content = readdlm(inFile, String)
input_file = [split(strip(line)) for line in content if length(strip(line)) > 0]

# Extract calculation parameters
Level_of_theory = input_file[1][1]
basis_set = input_file[1][2]
charge, multiplicity = parse.(Int, input_file[2])

# Extract molecular geometry
geom_file = input_file[3:end]
Atoms = [line[1] for line in geom_file]
geom_raw = [[parse(Float64, x) for x in line[2:4]] for line in geom_file]
geom = hcat(geom_raw...)

# Print extracted information
println(Level_of_theory, " ", basis_set, " ", charge, " ", multiplicity)
println(Atoms)
println(geom)

exps, coefs, origins, shells, norms = orbital_config(Atoms, geom_raw)

