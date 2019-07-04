
'''
usage:

'''

import sys
import numpy as np
import matplotlib as mpl;  mpl.use('Agg'); print( "plot WITHOUT Xserver" ); # this makes it run without Xserver (e.g. on supercomputer) # see http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
import matplotlib.pyplot as plt

# ============ Setup

Emin = -10.0; 
Emax = 5.0; 
dE   =  0.02;
save_bins = True

# ============ Functions

def projectionsToBins( kdata, Es, Emin=-6.0, Emax=0.0, dE=0.02 ):
    nks = kdata.shape[1]
    nE   = int((Emax-Emin)/dE) 
    bins = np.zeros( (nE, nks) )
    for i,Ei in enumerate(Es):
        iE = int((Ei-Emin)/dE)
        if (iE>0) and (iE<nE):
            bins[iE,:] = np.maximum( bins[iE,:], kdata[i,:] )
    extent=(0,1,Emin,Emax)
    return bins, extent

def print_attrs(name, obj):
    print( name )

def main():
    # ------- load data
    Es    = np.load('Es.npy') 
    kdata = np.load('klinesGL.npy')
    print ( "kdata.shape ", kdata.shape )
    kdata = np.sum( kdata, axis=1 )
    print ( "kdata.shape ", kdata.shape )
    # ------- gen bins 
    bins, extent = projectionsToBins( kdata, Es, Emin=Emin, Emax=Emax, dE=dE )
    if save_bins:
        np.savetxt('kmap_binsGL.dat', bins,  header=" nE %i nk %i dE %f Emin %f Emax %f " %(bins.shape[0],bins.shape[1],dE,Emin,Emax) )    
    # ------- plot
    plt.figure(figsize=(3,8))
    plt.imshow( np.log10(bins), interpolation='nearest', origin='image', extent=extent, cmap='jet' )
    plt.colorbar()
    plt.savefig( "kvazibandGL.png", bbox='tight' )
    plt.show()

if __name__ == "__main__":
    main()
