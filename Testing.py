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

#TEST DATA AND TEST RUN W/ SAMPLE DATA (TESTED: 4/27/15)
a=200
vects = np.linspace(0, 5,400)
x,y = np.meshgrid(vects, vects)
DistR = np.zeros((400, 400))
for i in range(a-21,a+21):
    for j in range(a-21,a+21):
        DistR[i,j] = np.sqrt((x[a,a] - x[i,j])**2 +(y[a,a] - y[i,j])**2)
        
T= np.arange(0,100,2)
R= np.arange(0, 3.4, .068)
interpol = scipy.interpolate.UnivariateSpline(R,T, k=5)
#xinterpol= np.linspace(0,4,100000)
TatR= interpol(DistR)
plt.imshow(TatR, interpolation='bicubic', origin='lower')
plt.colorbar()
plt.show()
