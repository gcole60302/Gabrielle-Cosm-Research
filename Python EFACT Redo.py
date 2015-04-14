##### Factor of redshift dependence for Hubble Constant

#Some imports
import matplotlib.pyplot as plt
import numpy as np
import math as m

#Some constants, n is the number of data points
OmmM = 0.27 * np.ones(n, dtype=float)
OmmK = x1 * np.ones(n, dtype=float)
OmmL = x2 * np.ones(n, dtype=float)

#Array of redshift values for clusters
z = np.arange(1, n, dtype=float)

#Define empty set for final output
EFFECT = np.zeros(n, dtype=float)

for i in range(0,n-1):
    EFFECT[i] = np.sqrt((OmmM[i])*(1+z[i])**3 + (OmmK[i])*(1 +z[i])**2 + (OmmL[i]))
#Here we used the full E(z) equation, but it may be trucated as wished    

print EFFECT

