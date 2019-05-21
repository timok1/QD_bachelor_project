
import sys
import numpy as np
import matplotlib.pyplot as plt

#sys.path.append("/home/prokop/git/WorkScripts/PYTHONPATH") 
import pyMolecular.atomicUtils as au
import pyMolecular.plotUtils as pu
 
# ---- setup
#  use:    python replace100_file.py ../common_data/Si538H.bas ../common_data/C2H4.bas 1 14

fin          =     sys.argv[1]
fgroup       =     sys.argv[2]
typ          = int(sys.argv[3])
ofTyp        = int(sys.argv[4])

# ---- main 

atoms  = np.genfromtxt  ( fin , skip_header=1 )
group  = np.genfromtxt( fgroup, skip_header=1 )

bonds,bondVecs = au.findAllBonds( atoms, Rcut=3.0, RvdwCut=0.8 )
neighs         = au.neighs( len(atoms), bonds )

select1 = au.findTypeNeigh( atoms, neighs, ofTyp, neighTyps={ typ:(2,2)} );    # print "select1", select1

select2 = au.getAllNeighsOfSelected( select1, neighs, atoms, typs={typ} );     # print "select2", select2
select2 = list(select2.keys())

pairs   = au.findPairs_one( select2, atoms, Rcut=2.5 );   #print "pairs:", pairs 

pairs   = au.pairsNotShareNeigh( pairs, neighs );   #print "pairs:", pairs

cog     = au.findCOG( atoms[:,1:] )

atoms_  = au.replacePairs( pairs, atoms, group, up_vec=(cog,1) );   # print( "atoms_ = ",atoms_ )
au.saveAtoms( atoms_, "replaced_100.xyz", xyz=True )

# ---------- ploting

'''
rotMat = au.makeRotMat( [1.0,1.0,1.0], [0.0,1.0,0.0] )
ps     = np.dot( rotMat, np.transpose(atoms[:,1:]) )


print( ps.shape )

#plt.plot( ps[0], ps[1], 'ok' )

#pu.plotAtoms( atoms[:,0], ps[0], ps[1] )
#pu.plotBonds( [ b[0] for b in bonds], ps[0], ps[1] )
#pu.plotBonds( pairs, ps[0], ps[1] )

plt.axis('equal')
plt.show()
'''
