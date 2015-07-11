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
import matplotlib.image as mpimg
import matplotlib as mpl
from matplotlib import pyplot
import scipy.interpolate

###############################################
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
y_const = (sigma_T)/(e_rm) #compt. parameter constant
G= 4.30117902e-9 #gravitational constant, km^2*Mpc*s^-2
h= 6.626066957e-34 #planck's constant kg*m^2*s^-1

###################################################################

#Arnaud profile constants
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
    return ((((np.sqrt((0.27)*(1. + z)**(3) + (1. - 0.27))))**(2))*(3.)*((72.)**(2.)))/((8.)*(np.pi)*(G))

def E_FACT(z):
    return np.sqrt((0.27)*(1. + z)**(3) + (1. - 0.27))

def E_FACT1(z):
    return (1.)/(np.sqrt((0.27)*(1. + z)**(3) + (1. - 0.27)))

def ANG_DIAM_DIST(z):
    return (((c)/((72.)*(1000.)))*(scipy.integrate.romberg(E_FACT1, 0, z)))/(1+z)

def ARRAY_AVG(x):
    s = np.sum(x[x!=0])
    l = np.float64(len(x[x!=0]))
    if l == 0:
        a = 0
    else:
        a = np.sum(x[x!=0])/np.float64(len(x[x!=0]))
    return a

def ARRAY_STD(x):
    s = x[x!=0]
    l = np.float64(len(x[x!=0]))
    if l == 1:
        a = 1e-15
    elif l == 0:
        a = 0.
    else:
        a = np.std(s)
    return a

def FREQ(f): # Temperature (mean CMB temperature), T, in Kelvin and frequency, f, in GHz 
    x = ((h)*((f)*(1.0e9)))/((k_b)*(2.73))
    return (((x)*((np.exp(x) + 1.)/(np.exp(x) - 1.)) - 4.))

###################################################################
#Here we define the output function
def R500ToArcmin(z, M500):
    R500 = ((M500)/((2500./3.)*(np.pi)*(Rho_Crit(z))))**(1./3.)
    #print R500
    x = np.arange(0,(100.)*(cutoff)*(R500)/(100.), 0.01)
    convert = np.arange(0,(100.)*(1)*(R500)/(100.), 0.01)
    q = np.zeros(len(x))
    for i in range(len(x)):
        y1= x[i]
        r = y1
        upperlim = np.sqrt(((cutoff)*(R500))**(2.) - (r)**(2.))
        def ARNAUD_PROJ(x1):
            return (1.)/(((((c500/R500)**2)*((x1**2. + y1**2.)))**(gamma/2.))*((1 + (((c500/R500)**2.)*(x1**2. + y1**2.))**(alpha/2))**(index)))
        if r < (cutoff)*(R500):
            q[i] = scipy.integrate.romberg(ARNAUD_PROJ,0.001, upperlim, divmax=20)
        elif r >(cutoff)*(R500):
            break
    PNORM = (1.65e-3)*((E_FACT(z))**(8./3.))*((((hubble70)*(M500))/(3.0e14))**(2./3. + alpha_p))*((hubble70)**2.)*((8.403)/(hubble70)**(3./2.))*(1e6)

    c = ((x)*(c500))/(R500)
    convert1 = ((convert)*(c500))/(R500)

    f = (y_const)*(q)*(PNORM)*(2.)*(mpc)

    r_over_r500= (c)/((c500))*(R500)
    convert2 = (convert1)/((c500))*(R500)

    absy_150ghz = f

    r_arcmin =(r_over_r500)/(ANG_DIAM_DIST(z))*(180.)/(np.pi)*(60.)
    convert3 = (convert2)/(ANG_DIAM_DIST(z))*(180.)/(np.pi)*(60.)

    dT_uK = (absy_150ghz)*(1.0e6)*(2.73)
    return np.max(convert3)

