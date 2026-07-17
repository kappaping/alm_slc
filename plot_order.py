## Plotting the orders on the lattice.

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


# Lattice structure.
ltype='ch'
Nbl=[16,16,1]
rs,Nr=ltc.ltcsites(ltype,Nbl)
bc=1
# Lattice generation can be slow. You may want to run cmtkit/lattice/main.py to generate the lattice file for repeated use.
# To read the saved lattice, set filet to the path of the lattice file then use toread=True to read it.
filet=''
NB,RD,RDV=ltc.ltcpairdist(ltype,rs,Nbl,bc,toread=False,filet=filet)
# Flavor and state.
Nfl=2
Nrfl=[Nr,Nfl]
# Filling fraction of each state.
nf=1./2.

# Setup of density matrix.
Ptype='read'
filet='almslc/checkerboard/t2_01_phi2_pi2_csgns_11_u0_40_16_16_1'
#filet='almslc/honeycomb/t2_01_phi2_pi2_csgns_1n11_u0_40_18_18_1'
#filet='almslc/bcc/t2_01_phi2_pi2_csgns_111_u0_40_8_8_8_1'
#filet='almslc/checkerboard3d/t2_01_phi2_pi2_csgns_11n11n11_u0_40_8_8_8_1'
#filet='almslc/diamond/t2_01_phi2_pi2_csgns_11n11n11_u0_40_8_8_8_1'
P=dm.setdenmat(Ptype,Nrfl,nf,fileti=filet,ltype=ltype)

# Consider only onsite orders.
Nst=tb.statenum(Nrfl)
for stid0 in range(Nst):
    for stid1 in range(Nst):
        if((stid0//Nfl)!=(stid1//Nfl)):P[stid0,stid1]=0.


sps=[dm.pairspin(P,rid,rid,Nfl) for rid in range(Nr)]
spas=[np.array([abs(sp[0]),abs(sp[1]),abs(sp[2])]) for sp in sps]
print('translation error =',sum(spas)/len(spas)-np.array([sqrt(sum([spa[nax]**2 for spa in spas])/len(spas)) for nax in range(3)]))

# Plot the orders. Select the setup based on the lattice.
rpls=[]
#rpls=[[[0,0,0],0],[[0,0,0],1],[[1,0,0],0],[[0,1,0],1]] # ch
#rpls=[[[0,1,0],0],[[0,1,0],1],[[1,1,0],0],[[1,1,0],1],[[0,2,0],0],[[1,0,0],1]] # ho
#rpls=[[[0,0,0],0],[[0,0,0],1]] # bcc
#rpls=[[[1,1,0],0],[[2,1,0],0],[[1,1,1],0],[[1,0,1],0],[[1,1,0],1],[[1,2,0],1],[[0,1,1],1],[[1,1,1],1]] # ch3d
Nnb=1
scl=1
res=50
show3d=False
plaz,plel,dist=0.,0.,None # 2D
#plaz,plel,dist=285.,75.,8. # ch3d
#plaz,plel,dist=280.,60.,8. # bcc
#plaz,plel,dist=235.,80.,8. # bcc2
#plaz,plel,dist=260.,70.,8. # dia
filetfigc='figs/fig_charge_order.pdf'
filetfigs='figs/fig_spin_order.pdf'
filetfig=[[filetfigc,filetfigs]]

pltc.plotorder(P,ltype,rs,Nrfl,Nbl,bc,NB,rpls=rpls,Nnb=Nnb,scl=scl,res=res,show3d=show3d,plaz=plaz,plel=plel,dist=dist,filetfig=filetfig)

