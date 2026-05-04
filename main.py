## Main function

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
Nbl=[16,16,1]
rs,Nr=ltc.ltcsites(ltype,Nbl)
bc=1
filet='../../data/lattice/checkerboard/16161_bc_1'
#filet='../../data/lattice/honeycomb/18181_bc_1'
#filet='../../data/lattice/bcc0/888_bc_1'
#filet='../../data/lattice/checkerboard3d/888_bc_1'
#filet='../../data/lattice/diamond/888_bc_1'
NB,RD,RDV=ltc.ltcpairdist(ltype,rs,Nbl,bc,toread=True,filet=filet)
# Flavor and state.
Nfl=2
Nrfl=[Nr,Nfl]
Nst=tb.statenum(Nrfl)
# Filling fraction of each state.
nf=(1./2.)
# Whether to adopt the Bogoliubov-de Gennes form.
tobdg=False

sys.stdout.flush()

# File name for writing out the density matrix.
filet='../../data/test/test00'

# Setup of initial density matrix.
Ptype='rand'
#fileti='../../data/test/test00'
toptb=False
toflrot=False
Pi=dm.setdenmat(Ptype,Nrfl,nf,tobdg=tobdg,fileti=fileti,toptb=toptb,toflrot=toflrot,ltype=ltype,rs=rs,Nbl=Nbl,NB=NB,RDV=RDV,Nbli=[12,12,1])
Pi=Pi.conj()

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
# Chemical potential
mu=hf.getchempot(H0,Pi,UINT,nf,Nst,tobdg=tobdg,dnf0=1./Nst**2,toprint=True,toread=False,filet=fileti)

sys.stdout.flush()

# Algorithm: Set the parameters and run the computation.
tofile=True
optm=1
printdm=20
writedm=40
Nhf=1000000
Nhfm=10

Pf=hf.hartreefock(Pi,H0,UINT,NB,Nrfl,nf,tofile=tofile,filet=filet,optm=optm,printdm=printdm,writedm=writedm,Nhf=Nhf,Ptype=Ptype,tobdg=tobdg,mu=mu)


