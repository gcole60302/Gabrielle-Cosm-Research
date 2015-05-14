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

def E_FACT1(z):
    return (1.)/(np.sqrt((0.27)*(1. + z)**(3) + (1. - 0.27)))

def ANG_DIAM_DIST(z):
    return (((c)/((72.)*(1000.)))*(scipy.integrate.romberg(E_FACT1, 0, z)))/(1+z)

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
    r_arcmin =(r_over_r500)/(ANG_DIAM_DIST(z))*(180.)/(np.pi)*(60.)
    dT_uK = (absy_150ghz)*(1.0e6)*(2.73)
    return dT_uK

##########################################
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
    r_arcmin =(r_over_r500)/(ANG_DIAM_DIST(z))*(180.)/(np.pi)*(60.)
    dT_uK = (absy_150ghz)*(1.0e6)*(2.73)
    return r_arcmin
##########################################    

#Here we define the output function
def MAP(z, M500):
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
    r_arcmin =(r_over_r500)/(ANG_DIAM_DIST(z))*(180.)/(np.pi)*(60.)
    dT_uK = (absy_150ghz)*(1.0e6)*(2.73)
    R = r_arcmin
    T = dT_uK
    MaxR = np.int8(np.ceil(np.max(R)))
    if MaxR %2 == 0:
        MaxR = MaxR +1
    else:
        MaxR = MaxR
    InterR =(MaxR)*(4)+1
    CenPo = (MaxR)*(2)
    vects = np.linspace(0,MaxR,InterR)
    x,y = np.meshgrid(vects, vects)
    DistR = np.zeros((InterR, InterR))
    for i in range(InterR):
        for j in range(InterR):
            DistR[i,j] = np.sqrt((x[CenPo,CenPo] - x[i,j])**2 +(y[CenPo,CenPo] - y[i,j])**2)
    interpol = scipy.interpolate.UnivariateSpline(R,T, k=5)
    xinterpol= np.linspace(0,4,10000000)
    TatR= interpol(DistR)
    plt.imshow(TatR, interpolation='bicubic', origin='lower')
    plt.colorbar()

######################################
def NMap():
    SIZE = 405
    vects = np.linspace(0,SIZE, SIZE*4+1)
    x,y = np.meshgrid(vects, vects)
    N1 = np.zeros((SIZE*4,SIZE*4))
    for i in range(SIZE*4):
        for j in range(SIZE*4):
            N1[i,j] = np.random.normal(0.0, 5.8)
    #plt.imshow(N1, origin='lower')
    #plt.colorbar()
    #plt.show()
    return N1

def FULLMAP(n):
#Empty grid size set (SIZE x SIZE)
    SIZE = 405
    vects = np.linspace(0,SIZE, SIZE*4+1)
    x,y = np.meshgrid(vects, vects)
#Empty grids of fixed size are set   
    N1 = np.zeros((n,SIZE*4,SIZE*4))
    TatR = np.zeros((n,SIZE*4,SIZE*4))
    SUM = np.zeros((SIZE*4,SIZE*4))
    for k in range(n):
#Random location details
        X = np.arange(45, len(x)-45)
        Y = np.arange(45, len(y)-45)
        CentX = np.random.choice(X,1)[0]
        CentY = np.random.choice(Y,1)[0]
#Random Cluster details
        z = np.arange(.7,2.9,0.1)
        m500 = np.arange(1e13,10e15, 1.3e9)
        Z = np.random.choice(z,1)[0]
        M500 = np.random.choice(m500,1)[0]
        R = PROFILER(Z,M500)
        T = (1)*PROFILET(Z,M500)
        print CentX, CentY, Z, M500
#Shape of smaller plot
        MaxR = np.int8(np.ceil(np.max(R)))
        if MaxR %2 == 0:
            MaxR = MaxR +1
        else:
            MaxR = MaxR

        InterR =((MaxR)*(10))/2
    
        for i in range(CentX-InterR, CentX+InterR):
            for j in range(CentY-InterR, CentY+InterR):
                N1[k,i,j] = np.sqrt(((x[CentX,CentY] +(1/8.)) - (x[i,j] + (1/8.)))**2 +((y[CentX,CentY] +(1/8.)) - (y[i,j] + (1/8.)))**2)
        interpol = scipy.interpolate.UnivariateSpline(R,T, k=5, ext=1)
        TatR[k] = interpol(N1[k])
        tatr = interpol(0)
        for i in range(len(TatR[k])):
            for j in range(len(TatR[k])):
                if TatR[k,i,j] == tatr:
                    TatR[k,i,j] = 0
                else:
                    continue
    for k in range(n):
        SUM = SUM+ TatR[k]
    
    plt.imshow(SUM+NMap(),interpolation='bicubic', origin='lower')
    plt.colorbar()
    plt.show()
    return SUM + NMap()













    
