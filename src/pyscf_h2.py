from pyscf import gto, scf, lib, cc,mp
import numpy as np

h2 = "H         0.000000    0.000000    0.0;\
H       0.7097959   0.000000   0.0"


mol = gto.M( atom = h2,
 basis =  "unc-321g")


mf = mol.RHF().run()

mf.MP2().run()