
import sys
import numpy as np
import pyMolecular.atomicUtils as au
#import pyMolecular.plotUtils as pu

# ---- setup
#  use:    python replace111_file.py ../common_data/Si538H.bas ../common_data/CH2OH.bas 1 14 2.0 1.83 1.0

fin       =     sys.argv[1]
fgroup    =     sys.argv[2]
typ       = int(sys.argv[3])
ofTyp     = int(sys.argv[4])
rcut      = float(sys.argv[5])
bond_length  = float(sys.argv[6])
prob         = float(sys.argv[7])

# ---- main

atoms = np.genfromtxt( fin, skip_header=1 )
group = np.genfromtxt( fgroup, skip_header=1 )

mask_Si = atoms[:,0]==ofTyp
mask_H  = atoms[:,0]==typ

bond_counts = au.countTypeBonds( atoms, atoms[mask_H], rcut );
mask_SiH1   = np.logical_and( ( bond_counts==1 ), mask_Si )
found, foundDict = au.findBondsTo( atoms, typ, atoms[mask_SiH1], rcut=rcut )
#atoms2 = au.replace( atoms.copy(), found, to=rep[0], bond_length=rep[1], prob=prob  )
atoms2 = au.replaceGroup( atoms, found, foundDict, group=group, Ups=(0.0,0.0,0.0), bond_length=bond_length, prob=prob )
au.saveAtoms( atoms2, "replaced_111.xyz" , xyz=True )
#au.saveAtoms( atoms2, "replaced_111.bas" , xyz=False )


