## Main function

from math import *
import numpy as np
import matplotlib.pyplot as plt
import joblib

import sys
sys.path.append('../cmtkit/lattice')
import lattice as ltc
sys.path.append('../cmtkit/tightbinding')
import tightbinding as tb
import densitymatrix as dm
import brillouinzone as bz
sys.path.append('../cmtkit/interaction')
import interaction as itn
sys.path.append('../cmtkit/bandtheory')
import bandtheory as bdth
sys.path.append('../cmtkit/plotlattice')
import plotband as plbd
sys.path.append('../cmtkit/hartreefock')
import hartreefock as hf
import slcham


# Lattice structure.
ltype='ho'
Nbl=[4,4,1]
rs,Nr=ltc.ltcsites(ltype,Nbl)
bc=1
#filet='../../data/lattice/checkerboard/16161_bc_1'
filet='../../data/lattice/honeycomb/18181_bc_1'
NB,RD,RDV=ltc.ltcpairdist(ltype,rs,Nbl,bc,toread=False,filet=filet)
# Flavor and state.
Nfl=2
Nrfl=[Nr,Nfl]
Nst=tb.statenum(Nrfl)
# Filling fraction of each state.
nf=(1./2.)
tobdg=False

# Setup of density matrix.
Ptype='copy'
filet='../../data/test/test00'
nbcpmax=2
Nbli=[18,18,1]
P=dm.setdenmat(Ptype,Nrfl,nf,fileti=filet,ltype=ltype,rs=rs,Nbl=Nbl,NB=NB,nbcpmax=nbcpmax,Nbli=Nbli)

# Sublattice-current model.
t1=1.
t2=0.1
phi2=pi/2.
tocsgns=True
#csgns=[1,1] # ch
csgns=[1,1,1] # ho, bcc
#csgns=[1,1,-1,1,-1,1] # ch3d, dia
ts=slcham.slcurrent(t1,t2,phi2,ltype,NB,RDV,rs,Nfl,tocsgns=tocsgns,csgns=csgns)
H0=tb.tbham(ts,NB,Nfl,rs)
# Interaction
us=[4.]
UINT=itn.interaction(NB,Nrfl,us)
H=hf.hfham(H0,P,UINT)
#H=H0

# Set the unit cell with periodicity prds.
prds=[1,1,1]
rucs,RUCRP=bdth.ftsites(ltype,rs,prds)

# Get the momentum-space Hamiltonian.
Hk=lambda k:bdth.ftham(k,H,Nrfl,RDV,rucs,RUCRP,tobdg=tobdg)

Nk=120

bzop=False
ks,dks=bz.listbz(ltype,prds,Nk,bzop)

#'''
todata=True
nbd,tonbd=0,False
dBfBfs=[bdth.berrycurv(k,Hk,dks,nf,tonbd=tonbd,nbd=nbd,tobdg=tobdg) for k in ks]
dBfs=[dBfBf[0] for dBfBf in dBfBfs]
Bfs=[dBfBf[1] for dBfBf in dBfBfs]
Ch=(1./(2*pi))*sum(dBfs)
print('Chern number = ',Ch)
data=[Bfs[nk] for nk in range(len(ks))]
dka=dBfBfs[0][2]
bzvol=len(ks)*dka
print('BZ volume = ',bzvol)
print('max(Bfs)=',np.max(np.array(Bfs)))
#for ndat in range(len(data)):
#    if(abs(data[ndat])>1):data[ndat]=np.sign(data[ndat])*0
'''
todata=False
bzvol=0
data=[]
'''

filetfig='../../figs/hartreefock/testbz.pdf'
tosave=True
tolabel=False
cmapt='PuOr'
cmapdarker=0.9
cmapmax=0.9
torevcmap=True
plbd.plotbz(ltype,prds,ks,todata=todata,data=data,ptype='gd',dks=dks,bzop=bzop,bzvol=bzvol,tolabel=tolabel,tosave=tosave,filetfig=filetfig,cmapt=cmapt,cmapdarker=cmapdarker,cmapmax=cmapmax,torevcmap=torevcmap)




