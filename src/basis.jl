function scf_hf(Atoms, geom)
    no_of_e, atomic_nos = no_of_electrons(Atoms)
    for i in 1:length(atomic_nos)
        atomic_nos[i] = float(atomic_nos[i])
    end
    atomic_nos = convert(Array{Float64}, atomic_nos)

    #---------------------------------------
    #     ENUC
    #--------------------------------------

    t1 = time()
    enuc = enuc_calculator(atomic_nos, geom)
    t2 = time()
    println(t2-t1)
    #----------------------------------------
    #    Orbitals
    #----------------------------------------
    exps, coefs, origins, shells, norms = orbital_config(Atoms, geom_raw)
    #println(exps)
    #-----------------------------------------
    #       Overlap integrals
    #-----------------------------------------
    t1 = time()
    Overlap_mat = S_mat(exps, coefs, origins, shells, norms)
    nbasis = size(Overlap_mat)[1]
    t2 = time()
    println("$t2-$t1 seconds for overlap integrals")
    #println(Overlap_mat)

    #----------------------------------------
    #     Kinetic integrals
    #----------------------------------------
    t1 = time()
    Kinetic_mat = T_mat(exps, coefs, origins, shells, norms)
    t2 = time()
    println("$t2-$t1 seconds for kinetic integrals")
    #println(Kinetic_mat)

    #---------------------------------------
    #   Coulomb integrals
    #---------------------------------------
    #println(size(atomic_nos)[1])
    t1 = time()
    Potential_mat = V_mat(exps, coefs, origins, shells, norms, atomic_nos, geom)
    t2 = time()
    println("$t2-$t1 seconds for coulomb integrals")
    core_hamil = [[Kinetic_mat[i,j] + Potential_mat[i,j] for j in 1:nbasis] for i in 1:nbasis]
    core_hamil = real(core_hamil)
    #println(core_hamil)
    #----------------------------------------
    #  Two elec integrals
    #----------------------------------------
    t1 = time()
    twoe, eri = Eri_mat(exps, coefs, origins, shells, norms)
    t2 = time()
    println("$t2-$t1 seconds for two elec integrals")

    EN, E, C, P, F = scf_iteration(1e-12, enuc, no_of_e, nbasis, Overlap_mat, core_hamil, twoe, false, true)

    #println(twoe)
end

t1 = time()
SCF = scf_hf(Atoms, geom)
t2 = time()
println("Total time: $(t2-t1)")
