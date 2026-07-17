## Plotting the BZ map of sublattice weight in the noninteracting models.

from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from matplotlib.patches import PathPatch
from matplotlib.text import Annotation
from matplotlib.path import Path
import matplotlib.colors as mcolors
import joblib

import sys
sys.path.append('../cmtkit/lattice')
import lattice as ltc
sys.path.append('../cmtkit/tightbinding')
import tightbinding as tb
import densitymatrix as dm
import brillouinzone as bz
sys.path.append('../cmtkit/bandtheory')
import bandtheory as bdth
sys.path.append('../cmtkit/plotlattice')
import plotband as plbd
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
t2=0.1
phi2=pi/2.
tocsgns=True
csgns=[1,1] # ch
#csgns=[1,1,1] # ho, bcc
ts=slcham.slcurrent(t1,t2,phi2,ltype,NB,RDV,rs,Nfl,tocsgns=tocsgns,csgns=csgns)
H=tb.tbham(ts,NB,Nfl,rs)

# Set the unit cell with periodicity prds.
prds=[1,1,1]
rucs,RUCRP=bdth.ftsites(ltype,rs,prds)

# Get the momentum-space Hamiltonian.
Hk=lambda k:bdth.ftham(k,H,Nrfl,RDV,rucs,RUCRP)

print('Finish model.')

todata=True
def sldiff(k):
    eigs=np.linalg.eigh(Hk(k))
    ees,eevs=eigs[0],eigs[1].conj().T
#    return sum([ees[nee]*np.linalg.multi_dot([eevs[nee],np.kron(tb.paulimat(3),tb.paulimat(0)),eevs[nee].conj().T]).real for nee in range(2)])
    return sum([np.linalg.multi_dot([eevs[nee],np.kron(tb.paulimat(3),tb.paulimat(0)),eevs[nee].conj().T]).real for nee in range(2)])

Nk=120
bzop=False
ks,dks=bz.listbz(ltype,prds,Nk,bzop)
data=[sldiff(k) for k in ks]

print('Finish data.')

filetfig='figs/fig_sublattice_weight.pdf'
tosave=True
tolabel=True
toclmax=True
cmapt,cmapdarker,cmapmax='PiYG',0.9,0.9
plbd.plotbz(ltype,prds,ks,todata=todata,data=data,ptype='gd',dks=dks,bzop=bzop,toclmax=toclmax,tolabel=tolabel,tosave=tosave,filetfig=filetfig,cmapt=cmapt,cmapdarker=cmapdarker,cmapmax=cmapmax)




