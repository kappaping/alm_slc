## Hartree-Fock computation of ground state.

from math import *
import numpy as np
import joblib


import sys
sys.path.append('../cmtkit/lattice')
import lattice as ltc
sys.path.append('../cmtkit/tightbinding')
import tightbinding as tb
import densitymatrix as dm
sys.path.append('../cmtkit/interaction')
import interaction as itn
sys.path.append('../cmtkit/hartreefock')
import hartreefock as hf
import slcham


# Lattice structure.
ltype='ch'
Nbl=[4,4,1]
rs,Nr=ltc.ltcsites(ltype,Nbl)
bc=1
# Lattice generation can be slow. You may want to run cmtkit/lattice/main.py to generate the lattice file for repeated use.
# To read the saved lattice, set filet to the path of the lattice file then use toread=True to read it.
filet = ""
NB,RD,RDV=ltc.ltcpairdist(ltype,rs,Nbl,bc,toread=False,filet=filet)
# Flavor and state.
Nfl=2
Nrfl=[Nr,Nfl]
Nst=tb.statenum(Nrfl)
# Filling fraction of each state.
nf=(1./2.)

sys.stdout.flush()

# File name for writing out the density matrix.
filet='data/denmat'

# Setup of initial density matrix. Ptype="rand" means start from a random density matrix. Ptype="read" reads from a file.
Ptype='rand'
fileti='data/denmat'
Pi=dm.setdenmat(Ptype,Nrfl,nf,fileti=fileti)

sys.stdout.flush()

# Sublattice-current model.
t1=1.
t2=0.1
phi2=pi/2.
tocsgns=True
csgns=[1,1] # ch
#csgns=[1,1,1] # ho, bcc
#csgns=[1,1,-1,1,-1,1] # ch3d, dia
ts=slcham.slcurrent(t1,t2,phi2,ltype,NB,RDV,rs,Nfl,tocsgns=tocsgns,csgns=csgns)
H0=tb.tbham(ts,NB,Nfl,rs)
# Interaction.
us=[4.]
UINT=itn.interaction(NB,Nrfl,us)

sys.stdout.flush()

# Algorithm: Set the parameters and run the computation.
tofile=True
optm=1
printdm=20
writedm=40
Nhf=1000000
Nhfm=10

Pf=hf.hartreefock(Pi,H0,UINT,NB,Nrfl,nf,tofile=tofile,filet=filet,optm=optm,printdm=printdm,writedm=writedm,Nhf=Nhf,Ptype=Ptype)


