## Main function

from math import *
import numpy as np
import matplotlib.pyplot as plt
import joblib

import sys
sys.path.append('../../cmt_code/lattice')
import lattice as ltc
sys.path.append('../../cmt_code/tightbinding')
import tightbinding as tb
import densitymatrix as dm
import brillouinzone as bz
sys.path.append('../../cmt_code/interaction')
import interaction as itn
sys.path.append('../../cmt_code/bandtheory')
import bandtheory as bdth
sys.path.append('../../cmt_code/plotlattice')
import plotband as plbd
import plotlattice as plltc

filet='../../figs/hartreefock/testcbar.pdf'
cmapt='PuOr'
cmapdarker=0.9
cmapmax=0.9
torevcmap=True
plltc.printcbar(filet=filet,cmapt=cmapt,cmapdarker=cmapdarker,cmapmax=cmapmax,torevcmap=torevcmap)
