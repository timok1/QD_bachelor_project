
import sys
import os
import subprocess
import argparse
import numpy as np

from functools import (partial, reduce)
from math import (pi, sqrt)
from multiprocessing import Pool
import itertools
from itertools import islice

from qmworks.common import InfoMO
from qmworks.utils import chunksOf
from qmworks.parsers.xyzParser import readXYZ
from qmworks.parsers.cp2KParser import readCp2KCoeff
from qmworks.utils import concat

sys.path.append('/home/dohnalov/git/nonAdiabaticCoupling')
import nac.integrals.realSpaceWf as rwf
from nac.schedule.components     import create_dict_CGFs
from nac import retrieve_hdf5_data

from typing import (Callable, List, Tuple)
Vector = np.ndarray
Matrix = np.ndarray

import matplotlib as mpl;  mpl.use('Agg'); print( "plot WITHOUT Xserver" ); # this makes it run without Xserver (e.g. on supercomputer) # see http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
import matplotlib.pyplot as plt



# ==== Setup
# for small cluster (approximately below 40 atoms of Si) one has to take into account that there is less than 100 orbitals (e.g. 88 for Si35)
#nCorr=13
#nOrb            = 200
#ilist = np.array(range(10)); ihomos=100-ilist-nCorr; ilumos=100+ilist+1-nCorr

#for clusters bigger than 40 atoms:
nOrb            = 200
ilist = np.array(range(10)); ihomos=100-ilist; ilumos=100+ilist+1

basis_name= "DZVP-MOLOPT-SR-GTH"
MOfname   = 'point_0-mo_coeff_0.out-1_0.MOLog'
xyzfname  = "relaxed.xyz"
hdf5fname = "/home/dohnalov/CP2K/cp2k_basis/quantum.hdf5"


# ==== Globals

rhists = { 'HOMO': [], 'LUMO': [] }
aspan  = None 

# ==== Functions

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

def wf2d( wf, axis, vmax=None ):
    #wf0     = np.zeros( wf.shape )
    #wfplus  = np.sum( np.maximum(  wf, wf0 ) , axis=axis )
    #wfminus = np.sum( np.maximum( -wf, wf0 ) , axis=axis )
    wfplus  =  np.amax(wf,axis=axis)
    wfminus = -np.amin(wf,axis=axis)
    if vmax is None:
        vmax1 = np.amax( wfplus )
        vmax2 = np.amax( wfminus )
        #print( vmax1, vmax2 )
        vmax  = np.amax( (vmax1,vmax2) )
    wfminus/=vmax; wfplus/=vmax;
    #img = np.stack( (wfminus,np.zeros(wfplus.shape),wfplus), axis=-1 )
    img = np.stack( (1-wfplus,1-0.5*(wfplus+wfminus),1-wfminus), axis=-1 )
    #print( img.shape )
    return img 
    
def radialHist( rho, XYZs, center ):
    #(nx,ny,nz)=rho.shape
    #end = (XYZs[0][-1,-1,-1],XYZs[1][-1,-1,-1],XYZs[2][-1,-1,-1])
    #R = np.sqrt( (XYZs[0]-0.5*nx)**2 + (XYZs[1]-0.5*nx)**2 + (XYZs[2]-0.5*nx)**2 )
    R = np.sqrt( (XYZs[0]-center[0])**2 + (XYZs[1]-center[0])**2 + (XYZs[2]-center[0])**2 )
    return np.histogram( R.flat, bins=100, weights=rho.flat)

def plotOrbXYZ( fname, wf, XYZs, E ):
    global aspan
    aspan = np.array( (XYZs[0][-1,-1,-1],XYZs[1][-1,-1,-1],XYZs[2][-1,-1,-1]) )
    plt.figure(figsize=(15,5))
    interpolation='nearest'
    plt.subplot(1,3,1); plt.imshow( wf2d(wf,0), interpolation=interpolation, extent=(0,aspan[1],0,aspan[2]) ); plt.title( "E=%f eV" %E )
    plt.subplot(1,3,2); plt.imshow( wf2d(wf,1), interpolation=interpolation, extent=(0,aspan[0],0,aspan[2]) )
    plt.subplot(1,3,3); plt.imshow( wf2d(wf,2), interpolation=interpolation, extent=(0,aspan[0],0,aspan[1]) )
    plt.savefig( fname+".png", bbox_inches='tight' )
    plt.close()
    hist,edges = radialHist( wf*wf, XYZs, aspan*0.5  )
    global rhists
    rhists[ fname[:4] ].append(hist)
    #plt.figure(figsize=(15,5))
    #plt.plot( edges[1:], hist )
    #plt.savefig( fname+"_Rdens.png", bbox_inches='tight' )
    #plt.close()
    
def plotRhists( rhists, prefix='rhist_HOMO' ):
    plt.figure(figsize=(5,5))
    rmax = np.sqrt( aspan[0]**2 + aspan[1]**2 + aspan[2]**2 ) * 0.5 
    rs   = np.linspace(0.0, rmax,len(rhists[0]))
    for i,rline in enumerate(rhists):
        plt.plot( rs, rline, label='%i' %i )
    plt.legend()
    plt.savefig( prefix+".png", bbox_inches='tight' )
    plt.close()

def main():
    angstroms_to_au = 1.889725989
    #nBasis = grepNumber( command='grep "Number of orbital functions:" cp2k.out | tail -1' )
    nBasis  = getNBasis( MOfname )
    print ( "nBasis ",  nBasis )

    CWD = os.getcwd()
    MOfile = os.path.join(CWD,MOfname)
    print ( MOfile )
    Es, MOs = readCp2KCoeff( MOfile, nOrb, nBasis )
    #print( "len(Es)", len(Es) )
    #print( "Es =", np.array(Es)*27.2114 )
    
    Ehomos = Es [  ihomos]*27.2114
    Elumos = Es [  ilumos]*27.2114
    homos  = np.transpose( MOs[:,ihomos] )
    lumos  = np.transpose( MOs[:,ilumos] )
    del MOs 
    print ( "Ehomos",Ehomos )
    print ( "Elumos",Elumos )
    
    atoms           = readXYZ( xyzfname )    
    dictCGFs        = create_dict_CGFs( hdf5fname, basis_name, atoms)
        
    fout = open( "dipoleTrans.dat",'w')
    Dmat = rwf.dipoleTransMatrix( homos, lumos, Ehomos, Elumos, atoms, dictCGFs, Rcut=6.0, dstep=np.array((0.5,0.5,0.5)), memLimit=14.0, fout=fout, plotFunc=plotOrbXYZ, outxsf=0 )
    fout.close()
    
    plotRhists( rhists['HOMO'], prefix='rhist_HOMO' )
    plotRhists( rhists['LUMO'], prefix='rhist_LUMO' )
    np.savetxt( 'rdens_HOMO.dat', np.transpose(np.array(rhists['HOMO'])) )
    np.savetxt( 'rdens_LUMO.dat', np.transpose(np.array(rhists['LUMO'])) )

    # Dictionary containing as key the atomic symbols and as values the set of CGFs

if __name__ == "__main__":
    main()
    
    
    
