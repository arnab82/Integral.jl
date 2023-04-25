# Read the contents of the file into a string
json_string = read("/home/arnabbachhar/Intgl/basis/3-21guc_H.json", String)

# Parse the JSON string into a Julia data structure
bs = JSON3.read(json_string)
nested_object_property = bs.elements["1"].electron_shells


function T_mat(exps::Vector{Float64}, coefs::Vector{Float64}, origins::Matrix{Float64}, shells::Matrix{Int64}, norms::Vector{Float64})
    nbasis = length(exps)
    tmat = zeros(nbasis, nbasis)
    exp1 = zeros(1)
    exp2 = zeros(1)
    coefs1 = zeros(1)
    coefs2 = zeros(1)
    norm1 = zeros(1)
    norm2 = zeros(1)
    s = 0.0
    for i in 1:nbasis
        for j in 1:i
            exp1[1] = exps[i]
            exp2[1] = exps[j]
            coefs1[1] = coefs[i]
            coefs2[1] = coefs[j]
            norm1[1] = norms[i]
            norm2[1] = norms[j]
            s = T(exp1, coefs1, shells[i, :], norm1, origins[i, :], exp2, coefs2, shells[j, :], norm2, origins[j, :])
            tmat[i, j] = tmat[j, i] = s
        end
    end
    return tmat
end