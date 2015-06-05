import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import scipy
import scipy.misc
import math
from sympy import *
import scipy.optimize as opt
plt.ion()


##################################################################

#Here's a whole bunch of constants
H0 = 72. #km / (s*Mpc)
OmegaM = 0.27 #matter density
hubble = H0/100.
hubble70 = H0/70.
mmw = 0.59 #mean molecular weight of gas
mu_e = 1.143 #mean molecular weight of electrons in cluster
m_p = 1.6726e-27 #proton mass in kg
m_e = 9.11e-31 #electron mass in kg
mpc = 3.08568025e22   #megaparsec in meters
Msol = 1.98892e30 #solar masses in kg
sigma_T = (6.652e-25)/((100.)**2) #thompson cross section in m^2
k_b = 1.3806503e-23 #Blotsz const in m^2*kg*s^-2*K^-1
c = 3.0e8 #speed of light in m*s^-1
q = 1.60217646e-19 #unit of charge in joules equal to electron volt
e_rm = ((m_e)*((c)**2))/((1000.)*(q)) #rest mass energy of electron in KeV
y_const = (sigma_T)/(e_rm) # compt. parameter constant
G = 4.30117902e-9 # gravitational constant, km^2*Mpc*s^-2
##############################################################

def Rho_Crit(z):
    return (((OmegaM)*(1. + z)**(3.) + (1. - OmegaM))*(3.)*((72.)**(2.)))/((8.)*(np.pi)*(G))

def E_FACT(z):
    return np.sqrt((OmegaM)*(1. + z)**(3.) + (1. - OmegaM))

###########################################
#All values choosen from cluster RXC J0003.8+0203
R500 = 0.879
P500 = 1.466
P_0 = 3.93
c500 = 1.33
alpha = 1.41
beta = 5.4905 
gamma = 0.567
index = (beta - gamma) / alpha
#all four used/calculated from arnaud paper(2010)

x = np.arange(1,1000,0.1) / (R500)
y = np.zeros(90)
def GNFW(x):
    return ((P_0)/(((c500*x)**gamma)*(1 + (c500*x)**alpha)**index))

plt.figure()
#plt.ylim(0,10e-1)
#plt.xlim(0,1000)
plt.xscale('log')
plt.yscale('log')
plt.plot(x,GNFW(x))

plt.show()

