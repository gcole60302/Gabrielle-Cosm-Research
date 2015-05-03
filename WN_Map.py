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

vects = np.linspace(0,5, 400)
N = np.zeros((400,400))
for i in range(400):
    for j in range(400):
        N[i,j] = np.random.normal(0.0, 1.8)

plt.imshow(N, origin='lower')

