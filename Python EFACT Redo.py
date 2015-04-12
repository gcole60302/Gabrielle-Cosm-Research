##### Factor of redshift dependence for Hubble Constant

#Some imports
import matplotlib.pyplot as plt
import numpy as np
import math as m

#Some constants
OmmM = 0.27
OmmK =
OmmL =

#Array of redshift values for clusters
z = np.arange(1, n)

#Define empty set for final output
EFFECT = np.zeros(n)

for i in range(len(z)):
    EFFECT[i] = m.sqrt(OmmM(1+z[i])**3 + OmmK(1 +z)**2 + OmmL)

