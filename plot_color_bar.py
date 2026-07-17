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
sys.path.append('../cmtkit/plotlattice')
import plotlattice as plltc

filet='figs/fig_color_bar.pdf'
#cmapt,cmapdarker,cmapmax,torevcmap='coolwarm',1.,1.,False
#cmapt,cmapdarker,cmapmax,torevcmap='PiYG',0.9,0.9,False
cmapt,cmapdarker,cmapmax,torevcmap='PuOr',0.9,0.9,True
plltc.printcbar(filet=filet,cmapt=cmapt,cmapdarker=cmapdarker,cmapmax=cmapmax,torevcmap=torevcmap)
