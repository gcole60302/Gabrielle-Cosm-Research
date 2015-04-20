import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.misc
import math
from sympy import *
import scipy.optimize as opt

#Here's a whole bunch of constants
H0=72. #km / (s*Mpc)
OmegaM= 0.27 #matter density
hubble = H0/100.
hubble70 = H0/70.
mmw = 0.59 #mean molecular weight of gas
mu_e = 1.143 #mean molecular weight of electrons in cluster
m_p =1.6726e-27 #proton mass in kg
m_e =9.11e-31 #electron mass in kg
mpc = 3.08568025e22   #megaparsec in meters
Msol = 1.98892e30 #solar masses in kg
sigma_T = (6.652e-25)/(100.)**2 #thompson cross section in m^2
k_b= 1.3806503e-23 #Blotsz const in m^2*kg*s^-2*K^-1
c = 3.0e8 #speed of light in m*s^-1
q= 1.60217646e-19 #unit of charge in joules equal to electron volt
e_rm = ((m_e)*(c)**2)/((1000)*(q)) #rest mass energy of electron in KeV
y_const = (sigma_T)/(e_rm) # compt. parameter constant
G= 4.30117902e-9 # gravitational constant, km^2*Mpc*s^-2



#Arnaud profile constants
P0 = 8.403/ (hubble70)**(1.5) #Pressure nought
#PNorm
alpha = 1.0510
beta = 5.4905
gamma = 0.3081
index = (beta - gamma) / alpha
c500 = 1.177 #charcteristic concentration
alpha_p = 0.12 #constant of correction see Arnaud (2010)


#Here are some functions we need 
def Rho_Crit(z):
    return ((((np.sqrt((0.27)*(1. + z)**(3) + (1. - 0.27))))**(2))*(3.)*((72.)**(2.)))/((8.)*(np.pi)*(G))

def E_FACT(z):
    return np.sqrt((0.27)*(1. + z)**(3) + (1. - 0.27))

def ANG_DIAM_DIST(z):
    return(((3e8)/(72))*(scipy.integrate.romberg(E_FACT, 0, z)))/(1 + z)

#def ARNAUD_PROJ(x1):
 #   return (1.)/(((((c500/R500)**2)*((x1**2. + y1**2.)))**(gamma/2.))*((1 + (((c500/R500)**2.)*(x1**2. + y1**2.))**(alpha/2))**(index)))

def PROFILE(z, M500):
    R500 = ((M500)/((2500./3.)*(np.pi)*(Rho_Crit(z))))**(1./3.)
    x = np.arange(0,(100.)*(6.)*(R500)/(100.), 0.01)
    q = np.zeros(len(x))
    for i in range(len(x)):
        y1= x[i]
        r = y1
        upperlim = np.sqrt(((6.)*(R500))**(2.) - (r)**(2.))
        def ARNAUD_PROJ(x1):
            return (1.)/(((((c500/R500)**2)*((x1**2. + y1**2.)))**(gamma/2.))*((1 + (((c500/R500)**2.)*(x1**2. + y1**2.))**(alpha/2))**(index)))
        if r < (6.)*(R500):
            q[i] = scipy.integrate.romberg(ARNAUD_PROJ,0.001, upperlim)
        elif r >(6.)*(R500):
            break
    PNORM = (1.65e-3)*((E_FACT(z))**(8./3.))*((((hubble70)*(M500))/(3.0e14))**(2./3. + alpha_p))*((hubble70)**2.)*((8.403)/(hubble70)**(3./2.))*(1e6)
    dT_uK = (PNORM)*(2.73)*(y_const)*(2)*(mpc)*(1.0e6)*(q)
    r_arcmin = (((x)*(np.pi)*(60.))/(((R500)**2.)*(ANG_DIAM_DIST(z))*(180.)))
    plt.plot(r_arcmin, dT_uK)
    plt.yscale('log')
    return plt.show()
    



