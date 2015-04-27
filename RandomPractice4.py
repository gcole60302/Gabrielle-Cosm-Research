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

def ANG_DIAM_DIST(z):
    return(((3e8)/(72))*(scipy.integrate.romberg(E_FACT, 0, z)))/(1 + z)

###################################################################

#Here we define the output function
def PROFILET(z, M500):
    R500 = ((M500)/((2500./3.)*(np.pi)*(Rho_Crit(z))))**(1./3.)
    x = np.arange(0,(100.)*(9.)*(R500)/(100.), 0.001)
    q = np.zeros(len(x))
    for i in range(len(x)):
        y1= x[i]
        r = y1
        upperlim = np.sqrt(((9.)*(R500))**(2.) - (r)**(2.))
        def ARNAUD_PROJ(x1):
            return (1.)/(((((c500/R500)**2)*((x1**2. + y1**2.)))**(gamma/2.))*((1 + (((c500/R500)**2.)*(x1**2. + y1**2.))**(alpha/2))**(index)))
        if r < (9.)*(R500):
            q[i] = scipy.integrate.romberg(ARNAUD_PROJ,0.001, upperlim, divmax=20)
        elif r >(9.)*(R500):
            break
    PNORM = (1.65e-3)*((E_FACT(z))**(8./3.))*((((hubble70)*(M500))/(3.0e14))**(2./3. + alpha_p))*((hubble70)**2.)*((8.403)/(hubble70)**(3./2.))*(1e6)
    c = ((x)*(c500))/(R500)
    f = (y_const)*(q)*(PNORM)*(2.)*(mpc)
    r_over_r500= (c)/((c500)*(R500))
    absy_150ghz = f
    r_arcmin =(r_over_r500)/(ANG_DIAM_DIST(z))*(180.)*(np.pi)*(60.)*(4.)
    dT_uK = (absy_150ghz)*(1.0e6)*(2.73)
    return dT_uK

##########################################    

#Here we define the output function
def PROFILER(z, M500):
    R500 = ((M500)/((2500./3.)*(np.pi)*(Rho_Crit(z))))**(1./3.)
    x = np.arange(0,(100.)*(9.)*(R500)/(100.), 0.001)
    q = np.zeros(len(x))
    for i in range(len(x)):
        y1= x[i]
        r = y1
        upperlim = np.sqrt(((9.)*(R500))**(2.) - (r)**(2.))
        def ARNAUD_PROJ(x1):
            return (1.)/(((((c500/R500)**2)*((x1**2. + y1**2.)))**(gamma/2.))*((1 + (((c500/R500)**2.)*(x1**2. + y1**2.))**(alpha/2))**(index)))
        if r < (9.)*(R500):
            q[i] = scipy.integrate.romberg(ARNAUD_PROJ,0.001, upperlim, divmax=20)
        elif r >(9.)*(R500):
            break
    PNORM = (1.65e-3)*((E_FACT(z))**(8./3.))*((((hubble70)*(M500))/(3.0e14))**(2./3. + alpha_p))*((hubble70)**2.)*((8.403)/(hubble70)**(3./2.))*(1e6)
    c = ((x)*(c500))/(R500)
    f = (y_const)*(q)*(PNORM)*(2.)*(mpc)
    r_over_r500= (c)/((c500)*(R500))
    absy_150ghz = f
    r_arcmin =(r_over_r500)/(ANG_DIAM_DIST(z))*(180.)*(np.pi)*(60.)*(4.)
    dT_uK = (absy_150ghz)*(1.0e6)*(2.73)
    return r_arcmin

######################################
#TEST DATA AND TEST RUN (WORKS 4/27/15)
vects = np.linspace(0, 5,21)
x,y = np.meshgrid(vects, vects)
DistR = np.zeros((21, 21))
for i in range(21):
    for j in range(21):
        DistR[i,j] = np.sqrt((x[10,10] - x[i,j])**2 +(y[10,10] - y[i,j])**2)
        
T= np.arange(0,100,2)
R= np.arange(0, 3.4, .068)
interpol = scipy.interpolate.UnivariateSpline(R,T, k=5)
xinterpol= np.linspace(0,4,100000)
TatR= interpol(DistR)
plt.imshow(TatR)
plt.colorbar()
######################################
#RArray = PROFILER(0.3,7.3e13)
#MaxR = np.max(RArray)
#RoundR = np.int8(np.ceil(MaxR))
#vects = np.linspace(0, RoundR,RoundR*4 +1)
#x,y = np.meshgrid(vects, vects)
#DistR = np.zeros((RoundR*4), (RoundR*4))
#for i in range(RoundR*4):
 #   for j in range(RoundR*4):
  #      DistR[i,j] = np.sqrt((x[RoundR/2,RoundR/2] - x[:j])**2 +(y[RoundR/2,RoundR/2] - y[:j])**2)
