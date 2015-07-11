import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.misc
import math
from sympy import *
import scipy.optimize as opt
from scipy import integrate
import scipy.optimize as opt
from matplotlib.ticker import NullFormatter, MaxNLocator
from numpy import linspace
import matplotlib.image as mpimg
import matplotlib as mpl
from matplotlib import pyplot
import scipy.interpolate
import scipy.stats as st
#from sympy.stats import *
#from scipy import stats

########################################
DATA = np.zeros((5,10))
DATA[0] = 1,2,3,4,5,6,7,8,9,0
DATA[1] = 2,4,6,8,1,3,5,7,9,1
DATA[2] = 3,1,2,4,5,4,5,2,3,5
DATA[3] = 7,7,7,8,8,9,5,6,4,3
DATA[4] = 2,3,2,2,1,1,1,5,5,5
DATA = DATA.T
print DATA[0]
print DATA[1]
print DATA[2]
print DATA
print np.average(DATA[0])
print np.var(DATA[0])
