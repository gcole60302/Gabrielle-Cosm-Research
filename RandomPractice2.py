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
plt.ion()

import numpy as np                                                                                                                                                                                                                                                             
import pylab as plt                                                                                                                                                                                                                                                            
import matplotlib.ticker as ticker                                                                                                                                                                                                                                             

# Generate data                                                                                                                                                                                                                                                                
x = np.linspace(0, 2*np.pi*1e-9)                                                                                                                                                                                                                                               
y = np.sin(x/1e-9)                                                                                                                                                                                                                                                             

# setup figures                                                                                                                                                                                                                                                                
fig = plt.figure()                                                                                                                                                                                                                                                             
ax1 = fig.add_subplot(121)                                                                                                                                                                                                                                                     
ax2 = fig.add_subplot(122)                                                                                                                                                                                                                                                     
# plot two identical plots                                                                                                                                                                                                                                                     
ax1.plot(x, y)                                                                                                                                                                                                                                                                 
ax2.plot(x, y)                                                                                                                                                                                                                                                                 

# Change only ax2                                                                                                                                                                                                                                                              
scale = 1e-9                                                                                                                                                                                                                                                                   
ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale))                                                                                                                                                                                                           
ax2.xaxis.set_major_formatter(ticks)                                                                                                                                                                                                                                           

ax1.set_xlabel("meters")                                                                                                                                                                                                                                                       
ax2.set_xlabel("nanometers")                                                                                                                                                                                                                                                   

plt.show()
