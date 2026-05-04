## Tight-binding square module

'''Tight-binding square module: Setup of tight-binding models on the square lattice'''

from math import *
import numpy as np

import sys
sys.path.append('../cmtkit/lattice')
import lattice as ltc
sys.path.append('../cmtkit/tightbinding')
import tightbinding as tb




def slcurrent(t1,t2,phi2,ltype,NB,RDV,rs,Nfl,tocsgns=False,csgns=[]):
    '''
    Cherckerboard-lattice model: Realization of trivial and Chern insulators.
    H = sum_<ij> c_i^dagger ( -t*sigmaz + i*lda*(sigma.d_{ij}) ) c_j + sum_i c_i^dagger ( m*sigmaz ) c_i
    2-flavor fermions: c=[c_1,c_2]
    t: Flavor-dependent hopping
    lda: Spin-obit-type flavor-exchange hopping with displacement d_{ij}=r_i-r_j
    m: Flavor-dependent effective mass
    For t=lda=1: Trivial if m>4 or m<-4; chern if -4<m<0 and 0<m<4; Semimetallic if m=-4,0,4
    '''
    print('Sublattice-current model: [t1,t2,phi2] =',[t1,t2,phi2],'current directions =',csgns)
    # Construct the functions for the hoppings.
    def tf(nb0,nb1,rs):
        if(nb0==nb1):return 0.*np.identity(Nfl)
        elif(NB[nb0,nb1]==1):return -t1*np.identity(Nfl)
        elif(NB[nb0,nb1]==2):
            blvs=ltc.blvecs(ltype)
            if(ltype in ['ch','ho']):blvs=blvs[0:2]
            if(ltype in ['ho','ch3d','dia']):
                blvs=list(blvs)
                blvst=[]
                for nblv0 in range(len(blvs)):
                    for nblv1 in range(len(blvs)-nblv0-1):
                        blvst+=[blvs[nblv0+nblv1+1]-blvs[nblv0]]
                blvs+=blvst
            blvs=[blv/np.linalg.norm(blv) for blv in blvs]
            if(tocsgns):blvs=[blvs[nblv]*csgns[nblv] for nblv in range(len(blvs))]
            rdv=RDV[nb0,nb1]
            rdv=rdv/np.linalg.norm(rdv)
            return -t2*sum([(abs(abs(np.dot(rdv,blv))-1.)<1e-14)*(e**(1.j*phi2*((-1)**rs[nb0][1])*np.sign(np.dot(rdv,blv)))) for blv in blvs])*np.identity(Nfl)
    return [tf,tf,tf]







