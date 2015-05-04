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
##################################
def WMAP():
    SIZE = 405
    vects = np.linspace(0,5, SIZE)
    N = np.zeros((SIZE,SIZE))
    for i in range(SIZE):
        for j in range(SIZE):
            N[i,j] = np.random.normal(0.0, 1.8)
    plt.imshow(N, origin='lower')
    return N

