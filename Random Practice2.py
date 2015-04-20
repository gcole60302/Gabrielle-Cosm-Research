import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.misc
import math
from sympy import *
import scipy.optimize as opt
from scipy import integrate

G= 4.30117902e-9 # gravitational constant, km^2*Mpc*s^-2

def Rho_Crit(z):
    return ((((np.sqrt((0.27)*(1. + z)**(3) + (1. - 0.27))))**(2))*(3.)*((72.)**(2.)))/((8.)*(np.pi)*(G))

def PROFILE(z, M500):
    return((M500)/((2500./3.)*(np.pi)*(Rho_Crit(z))))**(1./3.)
    


#def Plot(x,y):
 #   a = x*np.arange(20)
  #  b = y*np.arange(20)
   # plt.plot(a,b)
    #return plt.show()

#def E_FACT(z):
 #   return np.sqrt((0.27)*(1. + z)**(3) + (1. - 0.27))

#def ANG_DIAM_DIST(z):
 #   return(((3e8)/(72.))*(scipy.integrate.romberg(E_FACT, 0, z)))/(1. + z)

#def E_z(z,Omm_M):
 #   return np.sqrt((Omm_M)*(1. + z)**(3) + (1. - Omm_M))

#def Rho_Crit(z, H0, Omm_M):
 #   return (((E_z(z, Omm_M))**(2))*(3.)*((H0)**(2.)))/((8.)*(np.pi)*(G))
#y = np.arange(1,5,0.01)
#x = np.arange(2,43, 1)
#for i in range(3):
 #   print E_z(y[i], x[i])
    
#def PROFILE(z, M500):
    
    
