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
import bogoliubovdegennes as bdg
sys.path.append('../cmtkit/plotlattice')
import plotlattice as pltc
sys.path.append('../cmtkit/hartreefock')
import timedephartreefock as tdhf


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
# Filling fraction of each state.
nf=1./2.
# Whether to adopt the Bogoliubov-de Gennes form.
tobdg=False

# Setup of density matrix.
Ptype='read'
filet='../../data/test/test00'
toflrot=False
P=dm.setdenmat(Ptype,Nrfl,nf,fileti=filet,ltype=ltype,rs=rs,NB=NB,RDV=RDV,toflrot=toflrot,tobdg=tobdg)

'''
Nst=tb.statenum(Nrfl)
for stid0 in range(Nst):
    for stid1 in range(Nst):
        if((stid0//Nfl)!=(stid1//Nfl)):P[stid0,stid1]=0.
'''

sps=[dm.pairspin(P,rid,rid,Nfl) for rid in range(Nr)]
spas=[np.array([abs(sp[0]),abs(sp[1]),abs(sp[2])]) for sp in sps]
print('translation error =',sum(spas)/len(spas)-np.array([sqrt(sum([spa[nax]**2 for spa in spas])/len(spas)) for nax in range(3)]))

# Plot the orders.
rpls=[]
#rpls=[[[0,0,0],0],[[0,0,0],1],[[1,0,0],0],[[0,1,0],1]] # ch
#rpls=[[[0,1,0],0],[[0,1,0],1],[[1,1,0],0],[[1,1,0],1],[[0,2,0],0],[[1,0,0],1]] # ho
#rpls=[[[0,0,0],0],[[0,0,0],1]] # bcc
#rpls=[[[1,1,0],0],[[2,1,0],0],[[1,1,1],0],[[1,0,1],0],[[1,1,0],1],[[1,2,0],1],[[0,1,1],1],[[1,1,1],1]] # ch3d
Nnb=1
scl=1
res=50
#res=50
show3d=False
plaz,plel,dist=0.,0.,None
#plaz,plel,dist=285.,75.,8. # ch3d
#plaz,plel,dist=280.,60.,8. # bcc
#plaz,plel,dist=235.,80.,8. # bcc2
#plaz,plel,dist=260.,70.,8. # dia
filetfigc='/home/kappaping/research/figs/hartreefock/testfigc.pdf'
filetfigs='/home/kappaping/research/figs/hartreefock/testfigs.pdf'
filetfigo='/home/kappaping/research/figs/hartreefock/testfigo.pdf'
filetfigfe='/home/kappaping/research/figs/hartreefock/testfigfe.pdf'
filetfigfo='/home/kappaping/research/figs/hartreefock/testfigfo.pdf'
filetfig=[[filetfigc,filetfigs],[filetfigfe,filetfigfo]]

pltc.plotorder(P,ltype,rs,Nrfl,Nbl,bc,NB,rpls=rpls,Nnb=Nnb,scl=scl,res=res,show3d=show3d,plaz=plaz,plel=plel,dist=dist,filetfig=filetfig,tobdg=tobdg)

