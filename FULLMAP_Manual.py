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

def FREQ(f, T=2.73): # Temperature (mean CMB temperature), T, in micro Kelvin and frequency, f, in GHz 
    x = ((h)*((f)*(1.0e9)))/((k_b)*((T)/(1.0e6)))
    return (T)*(((x)*((np.exp(x) + 1.)/(np.exp(x) - 1.)) - 4.))

###################################################################

#Here we define the output function
def PROFILET(z, M500):
    R500 = ((M500)/((2500./3.)*(np.pi)*(Rho_Crit(z))))**(1./3.)
    x = np.arange(0,(100.)*(cutoff)*(R500)/(100.), 0.001)
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
    f = (y_const)*(q)*(PNORM)*(2.)*(mpc)
    r_over_r500= (c)/((c500))*(R500)
    absy_150ghz = f
    r_arcmin =(r_over_r500)/(ANG_DIAM_DIST(z))*(180.)/(np.pi)*(60.)
    dT_uK = (absy_150ghz)*(1.0e6)*(2.73)
    return dT_uK

##########################################
def PROFILER(z, M500):
    R500 = ((M500)/((2500./3.)*(np.pi)*(Rho_Crit(z))))**(1./3.)
    x = np.arange(0,(100.)*(cutoff)*(R500)/(100.), 0.001)
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
    f = (y_const)*(q)*(PNORM)*(2.)*(mpc)
    r_over_r500= (c)/((c500))*(R500)
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
    r_over_r500= (c)/((c500))*(R500)
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
def NMap(var):
    SIZE = 405
    vects = np.linspace(0,SIZE, SIZE*4+1)
    x,y = np.meshgrid(vects, vects)
    N1 = np.zeros((SIZE*4,SIZE*4))
    for i in range(SIZE*4):
        for j in range(SIZE*4):
            N1[i,j] = np.random.normal(0.0, var)
    return N1

def FULLMAP(n, M500, z, x, y):
#Define some empty arrays we need for later. One for cluster information
#the other for postage stamp cluster image data and the final one for
#the arrays to make final plots of temp as function of radial
    Cluster = np.zeros((n, 4))
    TOTAL_IMG = np.zeros((n, 63,65))
    TOTAL_PLT = np.zeros((n,2, 456))
#Empty grid size set (SIZE*4 x SIZE*4 archmin, with .25archmin pixels)
    SIZE = 405
    vects = np.linspace(0,SIZE, SIZE*4+1)
    x,y = np.meshgrid(vects, vects)
#Empty grids of fixed size are set   
    SUM = np.zeros((SIZE*4,SIZE*4))
    for k in range(n):
#Define empty arrays for later
        N1 = np.zeros((SIZE*4,SIZE*4))
        T_at_R = np.zeros((SIZE*4,SIZE*4))
#Location details
        X = x[k]
        Y = y[k]
#Cluster details
        Z = z[k]
        M500 = M500[k]
#Arnad Profifle for the cluster, and info saved to array
        R = PROFILER(Z,M500)
        T = (1)*PROFILET(Z,M500)
        Cluster[k] = X, Y, Z, M500
#Shape of smaller plot
        MaxR = np.int8(np.ceil(np.max(R)))
        if MaxR %2 == 0:
            MaxR = MaxR +1
        else:
            MaxR = MaxR
        Size = 2
        InterR =((MaxR)*(Size))/2
#Create radial distance array over proper range
        for i in range(X-InterR, X+InterR):
            for j in range(Y-InterR, Y+InterR):
                N1[i,j] = np.sqrt(((x[X,Y] +(1/8.)) - (x[i,j] + (1/8.)))**2 +((y[X,Y] +(1/8.)) - (y[i,j] + (1/8.)))**2)
#Begin the extrapolation
        interpol = scipy.interpolate.UnivariateSpline(R,T, k=5, ext=1)
        for i in range(X-InterR, X+InterR):
            for j in range(Y-InterR, Y+InterR):
                T_at_R[i,j] = interpol(N1[i,j])
        SUM = SUM + T_at_R
    SUM=SUM
    plt.figure()
    plt.imshow(SUM +NMap(30.0),interpolation= 'bicubic', origin='lower')
    plt.title('var 30')
    plt.colorbar()
    plt.show()
    T = SUM +NMap()
#Make cut outs around each cluster in the map
    for q in range(n):
        A = T[Cluster[q,0]-31:Cluster[q,0]+32,Cluster[q,1]-32:Cluster[q,1]+33]
        plt.figure()
        plt.imshow(A, interpolation= 'bicubic', origin='lower')
        plt.xlabel(r'$\mathrm{Arcmin}$',fontsize=16)
        plt.ylabel(r'$\mathrm{Arcmin}$',fontsize=16)
        plt.grid()
        cbar = plt.colorbar()
        cbar.set_label(r'$\mathrm{Temperature}\/\mathrm{(\mu K)}$',fontsize=16)
        plt.show()
#Make equal size grid to calculate radial distances again in archmin with quarter archmin pixels
        SIZE1 = len(A[0])
        SIZE2 = len(A)
        vects1 = np.linspace(0,SIZE1, SIZE1*4+1)
        x1,y1 = np.meshgrid(vects1, vects1)
        DIST = np.zeros((SIZE2,SIZE1))
        for i in range(SIZE2):
            for j in range(SIZE1):
                DIST[i,j] = np.sqrt(((x1[SIZE1/2,SIZE1/2] +(1/8.)) - (x1[i,j] + (1/8.)))**2 +((y1[SIZE1/2,SIZE1/2] +(1/8.)) - (y1[i,j] + (1/8.)))**2)
#Reshape DIST and cluster array, A into a single column.
#Establish empty array for radial/temperature values, R_T
#Establish empty array for radial distances, AVG_R
#stablish empty array for temperature values, AVG_T
        LIST1 = DIST.reshape(4095,1)
        LIST2 = A.reshape(4095,1)
        R_T = np.zeros((4095,2))
        AVG_R_T =np.zeros((4095,4095))
        AVG_T = np.zeros(4095)
        AVG_R = np.zeros(4095)
#Turn DIST and A into a joint two colunm array
        for i in range(4095):
            R_T[i,0] = LIST1[i]
        for i in range(4095):
            R_T[i,1] = LIST2[i]
#Create array of all temperature values associated with a given radial range
        for i in range(4095):
            t = R_T[i,0]
            AVG_R[i] = t
            for j in range(4095):
                if R_T[j,0] == t and R_T[j,0]!= 0:
                    AVG_R_T[i,j] = R_T[j,1]
                    R_T[j,0] = 0
                else:
                    continue
#Average all values in a given radial range
        for i in range(4095):
            AVG_T[i] = ARRAY_AVG(AVG_R_T[i])
#Plot averaged temperatures as a function of radial distances
        plt.figure()
        plt.ylabel(r'$\mathrm{Temperature} \/\mathrm{Decrement}\/\mathrm{(\mu K)}$', fontsize=16)
        plt.xlabel(r'$\mathrm{Radial}\/\mathrm{Distance}\/\mathrm{(Arcmin)}$', fontsize=16)
        plt.title(r'$\mathrm{Intergrated}\/\mathrm{Compton}\/\mathrm{Parameter}\/\mathrm{Decrement}\/\mathrm{Scatter}$', fontsize=18)
        plt.scatter(AVG_R[AVG_R!=0], AVG_T[AVG_T!=0])
#Save postage clusters into one array
#Save postage clusters radials and temperatures into one array
        TOTAL_IMG[q] = A
        TOTAL_PLT[q] = AVG_R[AVG_R!=0], AVG_T[AVG_T!=0]
    return TOTAL_IMG[0], TOTAL_IMG[1]










    
