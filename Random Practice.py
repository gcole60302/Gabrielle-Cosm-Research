##### Factor of redshift dependence for Hubble Constant

#Some imports
import matplotlib.pyplot as plt
import numpy as np
import math as m

#Some constants
OmmM = 0.27 * np.ones(100, dtype=object)
OmmK = 4 * np.ones(100, dtype=object)
OmmL = 9 * np.ones(100, dtype=object)

#Array of redshift values for clusters
z = np.arange(1, 100, dtype=object)

#Define empty set for final output
EFFECT = np.zeros(100, dtype=object)

for i in range(len(z)):
    EFFECT[i] = (OmmM[i])*(1+z[i])**3 + (OmmK[i])*(1 +z)**2 + (OmmL[i])

print EFFECT
