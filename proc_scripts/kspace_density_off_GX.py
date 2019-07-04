
import sys
sys.path.append('/home/dohnalov/git/nonAdiabaticCoupling')
basis_name="DZVP-MOLOPT-SR-GTH"
from functools import (partial, reduce)
from math import (pi, sqrt)
from multiprocessing import Pool
from nac.integrals.fourierTransform import calculate_fourier_trasform_cartesian, get_fourier_basis
from nac.schedule.components        import create_dict_CGFs
from qmworks.parsers.xyzParser import readXYZ
from qmworks.parsers.cp2KParser import readCp2KCoeff
from qmworks.utils import concat
from qmworks.common import InfoMO
from itertools import islice
from qmworks.utils import chunksOf
from nac import retrieve_hdf5_data
import h5py
import argparse
import itertools
import numpy as np
import os
import subprocess
from typing import (Callable, List, Tuple)


Vector = np.ndarray
Matrix = np.ndarray

# ============ Setup

angstroms_to_au = 1.889725989
lattice_const   = 5.40
nPoints         = 40
nOrb            = 200

MOfname   = 'point_0-mo_coeff_0.out-1_0.MOLog'
xyzfname  = "relaxed.xyz"
hdf5fname = "/home/dohnalov/CP2K/cp2k_basis/quantum.hdf5"

klines_glob = [
    (( 2.0,  0.0,  0.0),  ( 1.0,  1.0,  0.0)),
    (( 2.0,  0.0,  0.0),  ( 1.0, -1.0,  0.0)), 
    (( 2.0,  0.0,  0.0),  ( 1.0,  0.0,  1.0)), 
    (( 2.0,  0.0,  0.0),  ( 1.0,  0.0, -1.0)),
    (( 0.0,  2.0,  0.0),  ( 1.0,  1.0,  0.0)),
    (( 0.0,  2.0,  0.0),  (-1.0,  1.0,  0.0)),
    (( 0.0,  2.0,  0.0),  ( 0.0,  1.0,  1.0)),
    (( 0.0,  2.0,  0.0),  ( 0.0,  1.0, -1.0)),
    (( 0.0,  0.0,  2.0),  (-1.0,  0.0,  1.0)),    
    (( 0.0,  0.0,  2.0),  ( 1.0,  0.0,  1.0)),  
    (( 0.0,  0.0,  2.0),  ( 0.0,  1.0,  1.0)),  
    (( 0.0,  0.0,  2.0),  ( 0.0, -1.0,  1.0)),  
    ((-2.0,  0.0,  0.0),  (-1.0,  1.0,  0.0)),
    ((-2.0,  0.0,  0.0),  (-1.0, -1.0,  0.0)), 
    ((-2.0,  0.0,  0.0),  (-1.0,  0.0,  1.0)), 
    ((-2.0,  0.0,  0.0),  (-1.0,  0.0, -1.0)), 
    (( 0.0, -2.0,  0.0),  ( 1.0, -1.0,  0.0)),
    (( 0.0, -2.0,  0.0),  (-1.0, -1.0,  0.0)),
    (( 0.0, -2.0,  0.0),  ( 0.0, -1.0,  1.0)),
    (( 0.0, -2.0,  0.0),  ( 0.0, -1.0, -1.0)),
    (( 0.0,  0.0, -2.0),  ( 1.0,  0.0, -1.0)),    
    (( 0.0,  0.0, -2.0),  (-1.0,  0.0, -1.0)),  
    (( 0.0,  0.0, -2.0),  ( 0.0,  1.0, -1.0)),  
    (( 0.0,  0.0, -2.0),  ( 0.0, -1.0, -1.0)), 
    (( 1.0,  1.0,  1.0),  ( 0.0,  1.0,  1.0)),
    (( 1.0,  1.0,  1.0),  ( 1.0,  0.0,  1.0)),
    (( 1.0,  1.0,  1.0),  ( 1.0,  1.0,  0.0)),
    ((-1.0,  1.0,  1.0),  ( 0.0,  1.0,  1.0)),
    ((-1.0,  1.0,  1.0),  (-1.0,  0.0,  1.0)),
    ((-1.0,  1.0,  1.0),  (-1.0,  1.0,  0.0)),
    (( 1.0, -1.0,  1.0),  ( 0.0, -1.0,  1.0)),
    (( 1.0, -1.0,  1.0),  ( 1.0,  0.0,  1.0)), 
    (( 1.0, -1.0,  1.0),  ( 1.0, -1.0,  0.0)),
    ((-1.0, -1.0,  1.0),  ( 0.0, -1.0,  1.0)),
    ((-1.0, -1.0,  1.0),  (-1.0,  0.0,  1.0)),
    ((-1.0, -1.0,  1.0),  (-1.0, -1.0,  0.0)),
    (( 1.0,  1.0, -1.0),  ( 0.0,  1.0, -1.0)), 
    (( 1.0,  1.0, -1.0),  ( 1.0,  0.0, -1.0)), 
    (( 1.0,  1.0, -1.0),  ( 1.0,  1.0,  0.0)), 
    ((-1.0,  1.0, -1.0),  ( 0.0,  1.0, -1.0)), 
    ((-1.0,  1.0, -1.0),  (-1.0,  0.0, -1.0)), 
    ((-1.0,  1.0, -1.0),  (-1.0,  1.0,  0.0)), 
    (( 1.0, -1.0, -1.0),  ( 0.0, -1.0, -1.0)),
    (( 1.0, -1.0, -1.0),  ( 1.0,  0.0, -1.0)),
    (( 1.0, -1.0, -1.0),  ( 1.0, -1.0,  0.0)),
    ((-1.0, -1.0, -1.0),  ( 0.0, -1.0, -1.0)),
    ((-1.0, -1.0, -1.0),  (-1.0,  0.0, -1.0)), 
    ((-1.0, -1.0, -1.0),  (-1.0, -1.0,  0.0)), 
    ]


# ============ Functions

def grepNumber( command='grep "Number of occupied orbitals:" cp2k.out | tail -1' ):
    p = subprocess.Popen( command, stdout=subprocess.PIPE, shell=True)
    p.wait()
    out, err = p.communicate()
    n = int(out.split()[-1])
    return n
    
def getNBasis( fname ):
    n = 0
    with open( fname, 'r') as f:
        for line in f:
            words = line.split()
            if len(words) == 6:
                try:
                    n_ = int(words[0])
                    if n_ > n:
                        n = n_
                    else:
                        break 
                except:
                    continue
    return n

def main():
    nBasis  = getNBasis( MOfname )
    print ( "nBasis ",  nBasis )
    lower = 0
    upper = nOrb-1
    CWD = os.getcwd()
    MOfile = os.path.join(CWD,MOfname)
    print ( MOfile )
    Es, MOs = readCp2KCoeff( MOfile, nOrb, nBasis )
    print( "len(Es)", len(Es) )
    print( "Es =", np.array(Es)*27.2114 )
    
    atoms           = readXYZ( xyzfname )
    symbols         = np.array( [at.symbol for at in atoms] )
    coords_angstrom = np.array( [at.xyz    for at in atoms] )
    coords          = angstroms_to_au * coords_angstrom
    lattice_cte     = lattice_const * angstroms_to_au

    dictCGFs        = create_dict_CGFs( hdf5fname, basis_name, atoms)

    clin = np.linspace( 0.0, 1.0, nPoints )[:,None]
    kpoints = []
    klines = np.array(klines_glob)
    for kline in klines:
        kpoints.append( (1-clin)*kline[0][None,:] + clin*kline[1][None,:] )
    kpoints  = np.concatenate( kpoints )
    kpoints *= (1/lattice_cte)
    print ( "kpoints.shape ", kpoints.shape )
    
    orbitals = list(range(lower, upper + 1))
    dim_x    = len(orbitals)
    results  = np.empty((dim_x, len(klines), nPoints))
    
    print ("building basiset fourier dictionary ... ")
    chikdic = get_fourier_basis( symbols, dictCGFs, kpoints )
    print ("...fourier basis DONE !")
    for i, orb in enumerate(orbitals):
        print("Orbital: ", orb)
        mo_i      = MOs[:,orb]
        orbK      = calculate_fourier_trasform_cartesian( symbols, coords, dictCGFs, mo_i, kpoints, chikdic=chikdic )
        orbK      = np.absolute( orbK )
        orbK      = orbK.reshape( len(klines), -1 )
        print ("orbK.shape", orbK.shape)
        results[i] = orbK
    
    Es = [ Es[i]*27.2114 for i in orbitals ]
    np.save( 'klinesGX.npy', results )
    np.save( 'Es.npy', Es )
    
if __name__ == "__main__":
    #nBasis  = getNBasis( 'point_0-mo_coeff_0.out-1_0.MOLog' )
    #print( nBasis )
    main()
