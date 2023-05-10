"""

function compute_gradient(geom::Vector{Vector{Float64}}, C, S, H, F)
    
    Argument=>
    R:: Nuclear co-ordinates
    C::co-efficient Matrix
    S::overlap Matrix
    H::Core Hamiltonian
    F::Fock Matrix
    Return=>
    gradient of the energy
    
    n = length(R)
    P = compute_density_matrix(C, S)
    grad = zeros(n)
    for i in 1:n
        dP = zeros(size(P))
        for j in 1:(size(C, 1))
            for k in 1:(size(C, 1))
                for a in 1:(size(C, 2))
                    dP[j, k] += (C[j, a] * C[k, a] * S[a, a] * (R[i] - R[a])) / norm(R - C[:, a])
                end
            end
        end
        grad[i] = sum(dP .* (H + F))
    end
    return grad
end

function compute_hessian(R, C, S, H, F)

    Argument=>
    R:: Nuclear co-ordinates
    C::co-efficient Matrix
    S::overlap Matrix
    H::Core Hamiltonian
    F::Fock Matrix
    Return=>
    Hessian of the energy
    
    n = length(R)
    P = compute_density_matrix(C, S)
    hess = zeros(n, n)
    for i in 1:n
        for j in i:n
            d2P = zeros(size(P))
            for k in 1:(size(C, 1))
                for l in 1:(size(C, 1))
                    for a in 1:(size(C, 2))
                        d2P[k, l] += (C[k, a] * C[l, a] * S[a, a] * (R[i] - R[a]) * (R[j] - R[a])) / norm(R - C[:, a])^3
                    end
                end
            end
            hess[i, j] = sum(d2P .* (H + F)) + sum((P .* (H + F)) .* (R[i] - R[j])^2) / norm(R)^2
            hess[j, i] = hess[i, j]
        end
    end
    return hess
end


function hartree_fock_optimize(R, C0, S, H, F, maxiter=100, tol=1e-6)
    C = copy(C0)
    P = compute_density_matrix(C, S)
    E = compute_energy(R, C, S, H, F)
    grad = compute_gradient(R, C, S, H, F)
    hess = compute_hessian(R, C, S, H, F)

    for iter in 1:maxiter
        # Solve the Hessian equations to get the displacement
        displ = -inv(hess) * grad

        # Update the nuclear coordinates
        R_new = R + displ

        # Recalculate the molecular orbitals, density matrix, energy, gradient, and Hessian at the new geometry
        C_new, F_new = hartree_fock_solve(S, H, R_new)
        P_new = compute_density_matrix(C_new, S)
        E_new = compute_energy(R_new, C_new, S, H, F_new)
        grad_new = compute_gradient(R_new, C_new, S, H, F_new)
        hess_new = compute_hessian(R_new, C_new, S, H, F_new)

        # Check for convergence
        if norm(grad_new) < tol
            return R_new, E_new
        end

        # Update the variables for the next iteration
        R = R_new
        C = C_new
        F = F_new
        P = P_new
        E = E_new
        grad = grad_new
        hess = hess_new
    end

    # If we get here, the optimization did not converge
    println("WARNING: Optimization did not converge after $maxiter iterations")
    return R, E
end


function optimize_geometry(R, C, S, H, F, max_iter=100, conv_thresh=1e-5)
    for iter in 1:max_iter
        E = compute_energy(R, C, S, H, F)
        grad = compute_gradient(R, C, S, H, F)
        hess = compute_hessian(R, C, S, H, F)
        if norm(grad) < conv_thresh
            println("Converged in $iter iterations.")
            return R, E
        end
        ΔR = -inv(hess) * grad
        R += ΔR
    end
    println("Did not converge within $max_iter iterations.")
    return R, E
end


function vibrational_analysis(R, C0, S, H, F)
    # Calculate the initial energy and gradient
    C = copy(C0)
    P = compute_density_matrix(C, S)
    E = compute_energy(R, C, S, H, F)
    grad = compute_gradient(R, C, S, H, F)
    hess = compute_hessian(R, C, S, H, F)

    # Diagonalize the Hessian matrix to get the normal modes and frequencies
    w2, V = eigen(hess)
    w = sqrt(w2)

    # Convert the frequencies to wavenumbers
    cminv = 1 / (2 * pi * 2.99792458e10 * 100)
    freq = sort(cminv * w)

    # Print out the vibrational modes and frequencies
    println("Vibrational analysis:")
    println("====================")
    for i in 1:length(freq)
        println("Mode $i: ${freq[i]} cm^-1")
        println("Displacement       X           Y           Z")
        for j in 1:size(R, 1)
            disp = V[(j-1)*3+1:j*3, i]
            println("Atom $j        $(disp[1])    $(disp[2])    $(disp[3])")
        end
        println("====================")
    end
end
"""
"""include("./../Intgl.jl")

#the dipole_moment is not implemented

function dipole(a, lmn1, A, b, lmn2, B, C, direction)
    l1, m1, n1 = lmn1
    l2, m2, n2 = lmn2
    P = gaussian_product_center(a, A, b, B)
    if direction == "x"
        XPC = P[1] - C[1]
        D = E(l1, l2, 1, A[1]-B[1], a, b) + XPC*E(l1, l2, 0, A[1]-B[1], a, b)
        S2 = E(m1, m2, 0, A[2]-B[2], a, b)
        S3 = E(n1, n2, 0, A[3]-B[3], a, b)
        return D*S2*S3*sqrt(pi/(a+b))^3
    elseif direction == "y"
        YPC = P[2] - C[2]
        S1 = E(l1, l2, 0, A[1]-B[1], a, b)
        D = E(m1, m2, 1, A[2]-B[2], a, b) + YPC*E(m1, m2, 0, A[2]-B[2], a, b)
        S3 = E(n1, n2, 0, A[3]-B[3], a, b)
        return S1*D*S3*sqrt(pi/(a+b))^3
    elseif direction == "z"
        ZPC = P[3] - C[3]
        S1 = E(l1, l2, 0, A[1]-B[1], a, b)
        S2 = E(m1, m2, 0, A[2]-B[2], a, b)
        D = E(n1, n2, 1, A[3]-B[3], a, b) + ZPC*E(n1, n2, 0, A[3]-B[3], a, b)
        return S1*S2*D*sqrt(pi/(a+b))^3
    end
end
function Mu(a, b, C, direction)
    mu = 0.0
    for (ia, ca) in enumerate(a.coefs)
        for (ib, cb) in enumerate(b.coefs)
            mu += a.norm[ia]*b.norm[ib]*ca*cb*\
                     dipole(a.exps[ia],a.shell,a.origin,
                     b.exps[ib],b.shell,b.origin,C,direction)
        end
    end
    return mu
end
function RxDel(a, b, C, direction)
    l = 0.0
    for ia in 1:length(a.coefs)
        for ib in 1:length(b.coefs)
            l += a.norm[ia] * b.norm[ib] * a.coefs[ia] * b.coefs[ib] *
                 angular(a.exps[ia], a.shell, a.origin,
                         b.exps[ib], b.shell, b.origin, C, direction)
        end
    end
    return l
end"""

