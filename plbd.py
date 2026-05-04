## Main function

from math import *
import numpy as np
import joblib


import sys
sys.path.append('../cmtkit/lattice')
import lattice as ltc
import brillouinzone as bz
sys.path.append('../cmtkit/tightbinding')
import tightbinding as tb
import densitymatrix as dm
import bogoliubovdegennes as bdg
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
ltype='ch'
Nbl=[4,4,1]
rs,Nr=ltc.ltcsites(ltype,Nbl)
bc=1
#filet='../../data/lattice/checkerboard/16161_bc_1'
#filet='../../data/lattice/honeycomb/18181_bc_1'
#filet='../../data/lattice/bcc0/888_bc_1'
filet='../../data/lattice/dia/888_bc_1'
NB,RD,RDV=ltc.ltcpairdist(ltype,rs,Nbl,bc,toread=False,filet=filet)
# Flavor and state.
Nfl=2
Nrfl=[Nr,Nfl]
Nst=tb.statenum(Nrfl)
# Filling fraction of each state.
nf=(1./2.)

# Setup of density matrix.
Ptype='copy'
#filet='../../data/test/test00'
#filet='../../data/hartreefock/almslc/checkerboard/t2_01_phi2_pi2_csgns_11_u0_40_16_16_1'
filet='../../data/hartreefock/almslc/honeycomb/t2_01_phi2_pi2_csgns_1n11_u0_40_18_18_1'
#filet='../../data/hartreefock/almslc/bcc/t2_01_phi2_pi2_csgns_111_u0_40_8_8_8_1'
#filet='../../data/hartreefock/almslc/checkerboard3d/t2_01_phi2_pi2_csgns_11n11n11_u0_40_8_8_8_1'
#filet='../../data/hartreefock/almslc/diamond/t2_01_phi2_pi2_csgns_11n11n11_u0_40_8_8_8_1'
nbcpmax=2
Nbli=[16,16,1]
P=dm.setdenmat(Ptype,Nrfl,nf,fileti=filet,ltype=ltype,rs=rs,Nbl=Nbl,NB=NB,nbcpmax=nbcpmax,Nbli=Nbli)

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
# Interaction
us=[4.]
UINT=itn.interaction(NB,Nrfl,us)
#H=hf.hfham(H0,P,UINT)
H=H0

# Set the unit cell with periodicity prds.
prds=[1,1,1]
rucs,RUCRP=bdth.ftsites(ltype,rs,prds)

# Get the momentum-space Hamiltonian.
Hk=lambda k:bdth.ftham(k,H,Nrfl,RDV,rucs,RUCRP)
Nk=60

filetfig='../../figs/hartreefock/testbd.pdf'
tosave=True
sn=dm.pairspin(P,0,0,Nfl)
sn=sn/np.linalg.norm(sn)
yticks=[-4,0,4]
#cmapt,cmapdarker,cmapmax='coolwarm',1.,1.
cmapt,cmapdarker,cmapmax='PiYG',0.9,0.9
plbd.plotbandcontour(Hk,ltype,prds,Nfl,Nk,nf,datatype='sl',sn=sn,cttype='pm',tosave=tosave,filetfig=filetfig,yticks=yticks,cmapt=cmapt,cmapdarker=cmapdarker,cmapmax=cmapmax)


