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
from scipy.interpolate import interp1d
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
y_const = (sigma_T)/(e_rm) # compt. parameter constant
G= 4.30117902e-9 # gravitational constant, km^2*Mpc*s^-2

###################################################################

#Arnaud profile constants
P0 = 8.403/ (hubble70)**(1.5) #Pressure nought
alpha = 1.0510
beta = 5.4905
gamma = 0.3081
index = (beta - gamma) / alpha
c500 = 1.177 #charcteristic concentration
alpha_p = 0.12 #constant of correction see Arnaud (2010)

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

##########################################    

#Here we define the output function
def MAP(z, M500):
#Create the two arrays we need (radial distance and temp)
    R500 = ((M500)/((4./3.)*(np.pi)*(500.)*(Rho_Crit(z))))**(1./3.)
    PNORM = (1.65e-3)*((E_FACT(z))**(8./3.))*((((hubble70)*(M500))/(3.0e14))**(2./3. + alpha_p))*((hubble70)**2.)*((8.403)/((hubble70)**(3./2.)))*(1.0e6)
    x = np.arange(0,(100.)*(6.)*(R500)/(100.), 0.01)
    q = np.zeros(len(x))
    for i in range(len(x)):
        y1= x[i]
        r = y1
        upperlim = np.sqrt(((6.)*(R500))**(2.) - (r)**(2.))
        def ARNAUD_PROJ(x1):
            return (1.)/(((((c500/R500)**2.)*((x1**2. + y1**2.)))**(gamma/2.))*((1. + (((c500/R500)**2.)*(x1**2. + y1**2.))**(alpha/2))**(index)))
        if r < (6.)*(R500):
            q[i] = scipy.integrate.romberg(ARNAUD_PROJ,0.001, upperlim, divmax=20)
        elif r >(6.)*(R500):
            break
    c = ((x)*(c500))/(R500)
    f = (y_const)*(q)*(PNORM)*(2.)*(mpc)
    r_over_r500= (c)/((c500))*(R500)
    r_arcmin =(r_over_r500)/(ANG_DIAM_DIST(z))*(180.)/(np.pi)*(60.)
    dT_uK = (f)*(1.0e6)*(2.73)
    R = r_arcmin
    T = dT_uK
#Round the radial distance array up to the nearest odd number
    MaxR = np.int8(np.ceil(np.max(R)))
    if MaxR %2 == 0:
        MaxR = MaxR +1
    else:
        MaxR = MaxR
    Dim = np.int8(np.ceil((MaxR)/(np.sqrt(2))))
#Create meshgrid of diagonal distance approxiamately 2*MaxR by 2*MaxR
#by quarter units
    vects = np.linspace(0,2*Dim,2*Dim*4+1)
    x,y = np.meshgrid(vects, vects)
#Create empty 2d arrays for radial distance calculation and temp population
#with same dimesions of the meshgrid
    DistR = np.zeros((2*Dim*4+1, 2*Dim*4+1))
    T_at_R = np.zeros((2*Dim*4+1, 2*Dim*4+1))
#Populate DistR with radial distances from center to all other points
    for i in range(2*Dim*4+1):
        for j in range(2*Dim*4+1):
            DistR[i,j] = np.sqrt((x[Dim*4,Dim*4] - x[i,j])**2 +(y[Dim*4,Dim*4] - y[i,j])**2)
#Populate T_at_R with temp values corresponding to radial distances
    interpol = interp1d(R,T, kind='cubic', bounds_error=False, fill_value=0) 
    for i in range(2*Dim*4+1):
        for j in range(2*Dim*4+1):
            T_at_R[i,j] = interpol(DistR[i,j])

    plt.figure() 
    plt.xlabel(r'$\mathrm{Arcmin}$',fontsize=16)
    plt.ylabel(r'$\mathrm{Arcmin}$',fontsize=16)
    plt.imshow(T_at_R, interpolation='bicubic', origin='lower',vmin=0, vmax= np.max(T))
    plt.xticks(np.linspace(0, Dim*8+0,5), np.linspace(0, Dim*2,5))
    plt.yticks(np.linspace(0, Dim*8+0,5), np.linspace(0, Dim*2,5))
    cbar = plt.colorbar()
    cbar.set_label(r'$\mathrm{Temperature}\/\mathrm{(\mu K)}$',fontsize=16)

    return Dim

##############################################################
