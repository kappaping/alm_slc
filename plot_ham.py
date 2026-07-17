## Plotting the structure of the Hamiltonian on the lattice.

from math import *
import numpy as np
import joblib


import sys
sys.path.append('../cmtkit/lattice')
import lattice as ltc
sys.path.append('../cmtkit/tightbinding')
import tightbinding as tb
import densitymatrix as dm
import bogoliubovdegennes as bdg
sys.path.append('../cmtkit/plotlattice')
import plotlattice as pltc
import slcham

# Lattice structure.
ltype='ch'
Nbl=[4,4,1]
rs,Nr=ltc.ltcsites(ltype,Nbl)
bc=1
NB,RD,RDV=ltc.ltcpairdist(ltype,rs,Nbl,bc,toread=False)
# Flavor and state.
Nfl=2
Nrfl=[Nr,Nfl]
Nst=tb.statenum(Nrfl)
# Filling fraction of each state.
nf=(1./2.)

# Sublattice-current model.
t1=1.
t2=0.05
phi2=pi/2.
tocsgns=True
csgns=[1,1] # ch
#csgns=[1,1,1] # ho, bcc
#csgns=[1,1,-1,1,-1,1] # ch3d, dia
ts=slcham.slcurrent(t1,t2,phi2,ltype,NB,RDV,rs,Nfl,tocsgns=tocsgns,csgns=csgns)
H0=tb.tbham(ts,NB,Nfl,rs).conj()

# Plot the Hamiltonian. Select the setup based on the lattice.
rpls=[]
#rpls=[[[n0,n1,n2],sl] for n0 in [0] for n1 in [0] for n2 in [0] for sl in range(ltc.slnum(ltype))]
#rpls=[[[n0,n1,0],sl] for n0 in range(2) for n1 in range(2) for sl in range(ltc.slnum(ltype))]
rpls=[[[0,0,0],0],[[0,0,0],1],[[1,0,0],0],[[0,1,0],1]] # ch
#rpls=[[[0,1,0],0],[[0,1,0],1],[[1,1,0],0],[[1,1,0],1],[[0,2,0],0],[[1,0,0],1]] # ho
#rpls=[[[1,1,1],0],[[2,1,1],0],[[1,2,1],0],[[1,1,2],0]] # ch3d0
#rpls=[[[0,0,0],0],[[1,0,0],0],[[0,1,0],0],[[0,0,1],0]] # bcc0
#rpls=[[[0,0,0],0],[[1,0,0],0],[[0,1,0],0],[[0,0,1],0],[[0,0,0],1]] # bcc01
#rpls=[[[1,1,1],0],[[2,1,1],0],[[1,2,1],0],[[1,1,2],0]] # ch3d01
#rpls=[[[1,1,1],0],[[2,1,1],0],[[1,2,1],0],[[1,1,2],0],[[1,1,1],1],[[0,2,1],1],[[1,2,0],1],[[1,2,1],1]] # ch3d
#rpls=[[[0,0,0],0],[[1,0,0],0],[[0,1,0],0],[[0,0,1],0]] # dia
#rpls=[[[0,0,0],0],[[1,0,0],0],[[0,1,0],0],[[0,0,1],0],[[0,0,0],1],[[1,0,0],1],[[0,1,0],1],[[0,0,1],1]] # dia
Nnb=2
#res=10
res=50
show3d=False
plaz,plel,dist=0.,0.,None # 2D
#plaz,plel,dist=235.,80.,8. # bcc
#plaz,plel,dist=285.,75.,8. # ch3d
#plaz,plel,dist=260.,70.,8. # dia
filetfig='figs/fig_ham.pdf'

# Find out the n-th neighbor pairs.
nbidss=[ltc.nthneighbors(nnb,NB) for nnb in range(Nnb+1)]
# Compute the charge orders.
od4s=[[dm.orders(H0,nbidss,Nrfl,odtype='c')[0]]]
nbplidss,od4spl=pltc.collectplotelements(rs,nbidss,rpls,Nnb,od4s)

pltc.plotlattice(rs,Nnb,nbplidss,Nbl,ltype,bc,filetfig,'bislham',od4spl[0][0],res=res,show3d=show3d,plaz=plaz,plel=plel,dist=dist)

