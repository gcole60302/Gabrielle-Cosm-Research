#Name:
#               Arnaud Profile
#Purpose:
#               Calculate the Arnaud profile of a cluster of a given mass
#               and redshift. Calculate the SZ signal for a given cluster
#               using this profile.
#Imputs:
#   M500:       Cluster mass in units of Msun
#   Z:          Redshift of the cluster
#
#Outputs:
#   r_arcmin:   A radial array of coordinates in units of arcminutes
#   dT_uK:      An array of temperature decrements from CMB blackbody
#               in units of uK corresponding to r_arcmin
#
#Author:
#               Gabrielle
#
#Date Created:
#               Mar-2015 - June-2015
#
#Notes:
#
#
#
##################################################################
#Packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import scipy
import scipy.misc
import math
from sympy import *
import scipy.optimize as opt
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
import matplotlib as mpl
from matplotlib import pyplot
import scipy.interpolate
plt.ion()


##################################################################

#Here's a whole bunch of constants which will be used and their units where unclear
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
h= 6.626066957e-34 #planck's constant kg*m^2*s^-1

###################################################################

#Arnaud profile parameters
P0 = 8.403/ (hubble70)**(1.5) #Pressure nought
alpha = 1.0510
beta = 5.4905
gamma = 0.3081
index = (beta - gamma) / alpha
c500 = 1.177 #charcteristic concentration
alpha_p = 0.12 #constant of correction see Arnaud (2010)
cutoff  = 6. #radial range we intergrate out to 

###################################################################

#Here are some functions we need 
def Rho_Crit(z):
    return (((OmegaM)*(1. + z)**(3.) + (1. - OmegaM))*(3.)*((72.)**(2.)))/((8.)*(np.pi)*(G))

def E_FACT(z):
    return np.sqrt((OmegaM)*(1. + z)**(3.) + (1. - OmegaM))

def E_FACT1(z):
    return (1.)/(np.sqrt((OmegaM)*(1. + z)**(3.) + (1. - OmegaM)))

def ANG_DIAM_DIST(z):
    return (((c)/((72.)*(1000.)))*(scipy.integrate.romberg(E_FACT1, 0, z, divmax=20)))/(1+z)

def FREQ(a): # Temperature (mean CMB temperature), T, in Kelvin and frequency, f, in GHz 
    x = ((h)*((a)*(1.0e9)))/((k_b)*(2.73))
    return (((x)*((np.exp(x) + 1.)/(np.exp(x) - 1.)) - 4.))

###################################################################

#Here we define the output function
def PROFILE(z, M500, a):
    R500 = ((M500)/((4./3.)*(np.pi)*(500.)*(Rho_Crit(z))))**(1./3.)
    PNORM = (1.65e-3)*((E_FACT(z))**(8./3.))*((((hubble70)*(M500))/(3.0e14))**(2./3. + alpha_p))*((hubble70)**2.)*((8.403)/((hubble70)**(3./2.)))*(1.0e6)
    x = np.arange(0,(100.)*(cutoff)*(R500)/(100.), 0.01)
    q = np.zeros(len(x))
    for i in range(len(x)):
        y1= x[i]
        r = y1
        upperlim = np.sqrt(((cutoff)*(R500))**(2.) - (r)**(2.))
        def ARNAUD_PROJ(x1):
            return (1.)/(((((c500/R500)**2.)*((x1**2. + y1**2.)))**(gamma/2.))*((1. + (((c500/R500)**2.)*(x1**2. + y1**2.))**(alpha/2))**(index)))
        if r < (cutoff)*(R500):
            q[i] = scipy.integrate.romberg(ARNAUD_PROJ,0.001, upperlim, divmax=20)
        elif r >(cutoff)*(R500):
            break
    c = ((x)*(c500))/(R500)
    f = (y_const)*(q)*(PNORM)*(2.)*(mpc)
    r_over_r500 = (c)/((c500))*(R500)
    r_arcmin =(r_over_r500)/(ANG_DIAM_DIST(z))*(180.)/(np.pi)*(60.)
    dT_uK = (-1)*(f)*(1.0e6)*(2.73)*(FREQ(a))
    plt.figure()
    plt.plot(r_arcmin, dT_uK)
    #plt.text(5, 100 , r'$\mathrm{M_{500} =\/\ 4.5\ast 10^{14} \/\ M_{\odot}}$',fontsize=26)
    #plt.text(5, 25 , r'$\mathrm{Z =\/\ 0.22}$',fontsize=26)
    plt.legend((r'$\mathrm{M_{500} =\/\ 2.5\ast 10^{14} \/\ M_{\odot}}$',r'$\mathrm{M_{500} =\/\ 3.5\ast 10^{14} \/\ M_{\odot}}$',r'$\mathrm{M_{500} =\/\ 4.5\ast 10^{14} \/\ M_{\odot}}$'), 'upper right', shadow=False)
    plt.ylabel(r'$\mathrm{Comptonization,}\/\mathrm{y}\/\mathrm{(\mu K)}$', fontsize=16)
    plt.xlabel(r'$\mathrm{Radial}\/\mathrm{Distance}\/\mathrm{(Arcmin)}$', fontsize=16)
    plt.title(r'$\mathrm{Compton}\/\mathrm{Parameter}$', fontsize=18)
    plt.yscale('log')
    return 
    



