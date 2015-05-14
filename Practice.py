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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import scipy
import scipy.misc
import math
import matplotlib.image as mpimg
from numpy import linspace
from sympy import *
import scipy.optimize as opt
import scipy.interpolate
###########################################################
SIZE = 5

vect = np.linspace(0,SIZE,SIZE*4+1)
x,y = np.meshgrid(vect, vect)
R = np.zeros((2, SIZE*4,SIZE*4))
for k in range(2):
    for i in range(3,9):
        for j in range(3,9):
            R[k,i,j]= np.sqrt(((x[SIZE*2,SIZE*2]+ (1/8.)) - (x[i,j]+ (1/8.)))**2 +((y[SIZE*2,SIZE*2]+ (1/8.))- (y[i,j] + (1/8.)))**2)












