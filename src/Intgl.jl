include("./one_electron/overlap.jl")
include("./one_electron/one_electron_integral.jl")
include("./two_electron/two_electron.jl")
using SpecialFunctions
using Combinatorics: doublefactorial
using LinearAlgebra: norm, eigen
using StaticArrays
const ang2bohr = 1.8897261246257702
using PyCall
sp=pyimport("scipy.special")
hyp1f1=sp.hyp1f1
fact2=sp.factorial2
#println(fact2(-1))


function find_distance(a_coordi::Vector{Float64}, b_coordi::Vector{Float64})
    """attributes:
            coordinates of two atoms
        returns:
            return the bond-distance between them
    """
    x_1 = a_coordi[1]
    x_2 = b_coordi[1]
    y_1 = a_coordi[2]
    y_2 = b_coordi[2]
    z_1 = a_coordi[3]
    z_2 = b_coordi[3]
    R_square = (x_1 - x_2)^2 + (y_1 - y_2)^2 + (z_1 - z_2)^2
    R = sqrt(R_square)
    return R
end


function enuc(atomic_nos::Vector{Any}, geom)
    """
    attributes:
        atomic_nos: Vector of the atomic no's of the atoms
        geom: Array that contains the co-ordinates of the atoms
    return:
        Nuclear repulsion energy
    """
    n = length(atomic_nos)
    E_nuc = 0.0

    for i in 1:n
        for j in 1:i-1
            Z_a = atomic_nos[i]
            Z_b = atomic_nos[j]
            R_ab = find_distance(geom[i], geom[j])
            E_nuc += (Z_a * Z_b) / R_ab
        end
    end

    return E_nuc
end