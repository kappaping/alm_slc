## Main function

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
filet='../../data/lattice/checkerboard/16161_bc_1'
NB,RD,RDV=ltc.ltcpairdist(ltype,rs,Nbl,bc,toread=False,filet=filet)
# Flavor and state.
Nfl=2
Nrfl=[Nr,Nfl]
Nst=tb.statenum(Nrfl)
# Filling fraction of each state.
nf=(1./2.)
tobdg=False

# Sublattice-current model.
t1=1.
t2=0.1
phi2=pi/2.
tocsgns=True
csgns=[1,1] # ch
#csgns=[1,1,1] # ho, bcc
#csgns=[1,1,-1,1,-1,1] # ch3d, dia
ts=slcham.slcurrent(t1,t2,phi2,ltype,NB,RDV,rs,Nfl,tocsgns=tocsgns,csgns=csgns)
H=tb.tbham(ts,NB,Nfl,rs)

# Set the unit cell with periodicity prds.
prds=[1,1,1]
rucs,RUCRP=bdth.ftsites(ltype,rs,prds)

# Get the momentum-space Hamiltonian.
Hk=lambda k:bdth.ftham(k,H,Nrfl,RDV,rucs,RUCRP,tobdg=tobdg)

print('Finish model.')


print('Finish momenta.')

todata=True
def sldiff(k):
    eigs=np.linalg.eigh(Hk(k))
    ees,eevs=eigs[0],eigs[1].T
#    return sum([ees[nee]*np.linalg.multi_dot([eevs[nee],np.kron(tb.paulimat(3),tb.paulimat(0)),eevs[nee].conj().T]).real for nee in range(2)])
    return sum([np.linalg.multi_dot([eevs[nee],np.kron(tb.paulimat(3),tb.paulimat(0)),eevs[nee].conj().T]).real for nee in range(2)])

Nk=120
bzop=False
ks,dks=bz.listbz(ltype,prds,Nk,bzop)
data=[sldiff(k) for k in ks]
'''
hsks=bz.hskpoints(ltype,prds)
kk=np.max(np.abs(np.array([hsk[1] for hsk in hsks])))
dk=kk/Nk
k0s=np.arange(-kk,kk,dk)
k1s=np.arange(-kk,kk,dk)
K0s,K1s=np.meshgrid(k0s,k1s)
Nks=len(K0s)
ks=[K0s,K1s]
data=np.array([[sldiff(np.array([K0s[nk0,nk1],K1s[nk0,nk1],0.])) for nk1 in range(Nks)] for nk0 in range(Nks)])
'''

print('Finish data.')

filetfig='../../figs/hartreefock/testbz.pdf'
tosave=True
tolabel=True
toclmax=True
cmapt,cmapdarker,cmapmax='PiYG',0.9,0.9
plbd.plotbz(ltype,prds,ks,todata=todata,data=data,ptype='gd',dks=dks,bzop=bzop,toclmax=toclmax,tolabel=tolabel,tosave=tosave,filetfig=filetfig,cmapt=cmapt,cmapdarker=cmapdarker,cmapmax=cmapmax)




