# Intgl.jl
This program computes one and two electron integrals over gaussian functions. It can also perform Hartree fock and mp2 calculation.

Currently only sto-3g,sto-6g, def2-svp, def2-tzvp, 6-31g, 6-31g(d,p), 6-31g*, 6-31g**, 3-21g basis sets are implemented. All geometries should be included in Bohr unit.
# Example of Input File
h2o.inp

rhf sto-3g

0 1

O 0.000000000000 -0.143225816552 0.000000000000

H 1.638036840407 1.136548822547 -0.000000000000

H -1.638036840407 1.136548822547 -0.00000000000

For running the hartree-fock calculation clone the repo and change the root path of the basis set json files in the input.jl file,run (give proper path of the files) 

julia input.jl h2o.inp 

# Running H2 vibrational frequency 
Run  H2_frequency.jl file , the basis set available for this sto-3g,cc-pVDZ,3-21g-uc ,DZP
