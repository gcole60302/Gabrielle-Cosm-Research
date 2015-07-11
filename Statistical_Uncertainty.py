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
def PROFILET(z, M500):
    R500 = ((M500)/((2500./3.)*(np.pi)*(Rho_Crit(z))))**(1./3.)
    #print R500
    x = np.arange(0,(100.)*(cutoff)*(R500)/(100.), 0.01)
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
    x = np.arange(0,(100.)*(cutoff)*(R500)/(100.), 0.01)
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
    x = np.arange(0,(100.)*(9.)*(R500)/(100.), 0.01)
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
#This function creates a white noise map of a specific deivation
#it can be modified to depend on the mean and variantion, but remember
#to also change it within the FULLMAP function
#Noise is per pixel and pixels are quarter archmin scale
#Noise is in units of microK
def NMap():
    SIZE = 405
    vects = np.linspace(0,SIZE, SIZE*4+1)
    x,y = np.meshgrid(vects, vects)
    N1 = np.zeros((SIZE*4,SIZE*4))
    for i in range(SIZE*4):
        for j in range(SIZE*4):
            N1[i,j] = np.random.normal(0.0, 18.0)
    #plt.imshow(N1, origin='lower')
    #plt.colorbar()
    #plt.show()
    return N1

def STAT_UNCERT(n):
#Define some empty arrays we need for later. One for cluster information
#the other for postage stamp cluster image data and the final one for
#the arrays to make final plots of temp as function of radial
    Cluster = np.zeros((n, 4))
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
#Random location details
        X = np.arange(100, len(x)-101)
        Y = np.arange(100, len(y)-101)
        CentClusIndexA = np.random.choice(X,1)[0]
        CentClusIndexB = np.random.choice(Y,1)[0]
#Random Cluster details, we use normal distributions about means 0.5 for
#the redshift and 5e14 for the mass, standard deviations are as coded
        z = np.arange(0.2,1.2,.1)
        Z = 0.5 #np.random.choice(z,1)[0]
        M500 = 3.2e14 #np.abs(np.random.normal(3.5e14, 1.0e14, 1)[0])
        R = PROFILER(Z,M500)
        T = (1)*PROFILET(Z,M500)
        Cluster[k] = CentClusIndexA, CentClusIndexB, Z, M500
        #print Z, M500
#Shape of smaller plot
        MaxR = np.int8(np.ceil(np.max(R)))
        if MaxR %2 == 0:
            MaxR = MaxR +1
        else:
            MaxR = MaxR
        Size = 2
        InterR =((MaxR)*(Size))/2
#Create radial distance array over proper range
        for i in range(CentClusIndexA-InterR, CentClusIndexA+InterR):
            for j in range(CentClusIndexB-InterR, CentClusIndexB+InterR):
                N1[i,j] = np.sqrt(((x[CentClusIndexA,CentClusIndexB] +(1/8.)) - (x[i,j] + (1/8.)))**2 +((y[CentClusIndexA,CentClusIndexB] +(1/8.)) - (y[i,j] + (1/8.)))**2)
#Begin the extrapolation
        interpol = scipy.interpolate.UnivariateSpline(R,T, k=5, ext=1)
        for i in range(CentClusIndexA-InterR, CentClusIndexA+InterR):
            for j in range(CentClusIndexB-InterR, CentClusIndexB+InterR):
                T_at_R[i,j] = interpol(N1[i,j])
        SUM = SUM + T_at_R
    SUM=SUM
    #plt.figure()
    #plt.imshow(SUM +NMap(),interpolation= 'bicubic', origin='lower')
    #plt.title('var 18')
    #plt.colorbar()
    #plt.show()
    T = SUM +NMap()
#Make cut outs around each cluster in the map
    S = T[Cluster[0,0]-31:Cluster[0,0]+32,Cluster[0,1]-32:Cluster[0,1]+33]
    SIZE1 = len(S[0])
    SIZE2 = len(S)
    TOTAL_PLT = np.zeros((n, 2, 457))#last number is the length of AVG_R[AVG_R!=0] 749 for 4, 1108 for 5 and 457 for 3 and 1530 for 6
    TOTAL_IMG = np.zeros((n, SIZE2, SIZE1))
    for q in range(n):
        A = T[Cluster[q,0]-31:Cluster[q,0]+32,Cluster[q,1]-32:Cluster[q,1]+33]

#Make equal size grid to calculate radial distances again in archmin with quarter archmin pixels
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
        LIST1 = DIST.reshape(SIZE1*SIZE2,1)
        LIST2 = A.reshape(SIZE1*SIZE2,1)
        R_T = np.zeros((SIZE1*SIZE2,2))
        AVG_R_T =np.zeros((SIZE1*SIZE2,SIZE1*SIZE2))
        AVG_T = np.zeros(SIZE1*SIZE2)
        AVG_R = np.zeros(SIZE1*SIZE2)
#Turn DIST and A into a joint two colunm array
        for i in range(SIZE1*SIZE2):
            R_T[i,0] = LIST1[i] #Radial distances
            R_T[i,1] = LIST2[i] #Temps
#Create array of all temperature values associated with a given radial range
        for i in range(SIZE1*SIZE2):
            r = R_T[i,0]
            AVG_R[i] = r
            for j in range(SIZE1*SIZE2):
                if R_T[j,0] == r and R_T[j,0]!= -1.:
                    AVG_R_T[i,j] = R_T[j,1]
                    R_T[j,0] = -1.
                else:
                    continue
#Average all values in a given radial range
        for i in range(SIZE1*SIZE2):
            AVG_T[i] = ARRAY_AVG(AVG_R_T[i])
#Plot averaged temperatures as a function of radial distances
        #plt.figure()
        #plt.ylabel(r'$\mathrm{Temperature}\/\mathrm{(\mu K)}$', fontsize=16)
        #plt.xlabel(r'$\mathrm{Radial}\/\mathrm{Distance}\/\mathrm{(Arcmin)}$', fontsize=16)
        #plt.title(str(Cluster[q,2:4]))
        ##plt.title(r'$\mathrm{Intergrated}\/\mathrm{Compton}\/\mathrm{Parameter}\/\mathrm{Decrement}\/\mathrm{Scatter}$', fontsize=18)
        #plt.scatter(AVG_R[AVG_R!=0], AVG_T[AVG_T!=0], label='Scatter')
        #plt.plot(PROFILER(Cluster[q,2], Cluster[q,3]), PROFILET(Cluster[q,2], Cluster[q,3]), 'r', label='Arnaud Profile')
        #plt.legend(loc='upper right', shadow=False)
#Save the postage size image data for each reapective cluster into one array
#Save postage clusters radial arrays and temperatures arrays into one array for each respective cluster
        TOTAL_IMG[q] = A
        TOTAL_PLT[q] = AVG_T[AVG_T !=0.], AVG_T[AVG_T !=0.]
#################################################
#Create some empty arrays to begin the binning process
    BIN1 = np.zeros((SIZE2,SIZE1))
    BIN2 = np.zeros((SIZE2,SIZE1))
    BIN3 = np.zeros((SIZE2,SIZE1))
    BIN4 = np.zeros((SIZE2,SIZE1))
    BIN5 = np.zeros((SIZE2,SIZE1))
    ALL_BINS = np.zeros((5,SIZE2,SIZE1))  
#We bin the data according to various features (Mass/Redshift, etc) to make the bined heat maps
#2 is for redshift, 3 is for mass
    for i in range(n):
        if Cluster[i,3] <= 5e14:
            BIN1 += TOTAL_IMG[i]
        if Cluster[i,3] <= 5e14 and  Cluster[i,3] >= 4e14:
            BIN2 += TOTAL_IMG[i]
        if Cluster[i,3] <= 4e14 and  Cluster[i,3] >= 3e14:
            BIN3 += TOTAL_IMG[i]
        if Cluster[i,3] <= 3e14 and  Cluster[i,3] >= 2e14:
            BIN4 += TOTAL_IMG[i]
        if Cluster[i,3] <= 2e14:
            BIN5 += TOTAL_IMG[i]
    ALL_BINS[0] = BIN1
    ALL_BINS[1] = BIN2
    ALL_BINS[2] = BIN3
    ALL_BINS[3] = BIN4
    ALL_BINS[4] = BIN5
    for k in range(1):
        A = ALL_BINS[k]
        vects1 = np.linspace(0,SIZE1, SIZE1*4+1)
        x1,y1 = np.meshgrid(vects1, vects1)
        DIST = np.zeros((SIZE2,SIZE1))
        for i in range(SIZE2):
            for j in range(SIZE1):
                DIST[i,j] = np.sqrt(((x1[SIZE1/2,SIZE1/2] +(1/8.)) - (x1[i,j] + (1/8.)))**2 +((y1[SIZE1/2,SIZE1/2] +(1/8.)) - (y1[i,j] + (1/8.)))**2)
        LIST1 = DIST.reshape(SIZE1*SIZE2,1)
        LIST2 = A.reshape(SIZE1*SIZE2,1)
        R_T = np.zeros((SIZE1*SIZE2,2))
        BINS = 20
        BINT = np.zeros((BINS, SIZE1*SIZE2))
        BINR = np.zeros((BINS, SIZE1*SIZE2))
        for i in range(SIZE1*SIZE2):
            R_T[i,0] = LIST1[i] #Radial distances
            R_T[i,1] = LIST2[i] #Temps
        for i in range(SIZE1*SIZE2):
            if R_T[i,0] <= .50:
                BINT[0,i] = R_T[i,1]
                BINR[0,i] = R_T[i,0]
            elif R_T[i,0] > .50 and R_T[i,0] <= 0.75:
                BINT[1,i] = R_T[i,1]
                BINR[1,i] = R_T[i,0]
            elif R_T[i,0] > 0.75 and R_T[i,0] <= 1.0:
                BINT[2,i] = R_T[i,1]
                BINR[2,i] = R_T[i,0]
            elif R_T[i,0] > 1.0 and R_T[i,0] <= 1.5:
                BINT[3,i] = R_T[i,1]
                BINR[3,i] = R_T[i,0]
            elif R_T[i,0] > 1.5 and R_T[i,0] <= 2.0:
                BINT[4,i] = R_T[i,1]
                BINR[4,i] = R_T[i,0]
            elif R_T[i,0] > 2.0 and R_T[i,0] <= 3.0:
                BINT[5,i] = R_T[i,1]
                BINR[5,i] = R_T[i,0]
            elif R_T[i,0] > 3.0 and R_T[i,0] <= 4.0:
                BINT[6,i] = R_T[i,1]
                BINR[6,i] = R_T[i,0]
            elif R_T[i,0] > 4.0 and R_T[i,0] <= 5.0:
                BINT[7,i] = R_T[i,1]
                BINR[7,i] = R_T[i,0]
            elif R_T[i,0] > 5.0 and R_T[i,0] <= 6.0:
                BINT[8,i] = R_T[i,1]
                BINR[8,i] = R_T[i,0]
            elif R_T[i,0] > 6.0 and R_T[i,0] <= 6.5:
                BINT[9,i] = R_T[i,1]
                BINR[9,i] = R_T[i,0]
            elif R_T[i,0] > 6.5 and R_T[i,0] <= 7.0:
                BINT[10,i] = R_T[i,1]
                BINR[10,i] = R_T[i,0]
            elif R_T[i,0] > 7.0 and R_T[i,0] <= 7.5:
                BINT[11,i] = R_T[i,1]
                BINR[11,i] = R_T[i,0]
            elif R_T[i,0] > 7.5 and R_T[i,0] <= 8.0:
                BINT[12,i] = R_T[i,1]
                BINR[12,i] = R_T[i,0]
            elif R_T[i,0] > 8.0 and R_T[i,0] <= 8.5:
                BINT[13,i] = R_T[i,1]
                BINR[13,i] = R_T[i,0]
            elif R_T[i,0] > 8.5 and R_T[i,0] <= 9.0:
                BINT[14,i] = R_T[i,1]
                BINR[14,i] = R_T[i,0]
            elif R_T[i,0] > 9.0 and R_T[i,0] <= 9.5:
                BINT[15,i] = R_T[i,1]
                BINR[15,i] = R_T[i,0]
            elif R_T[i,0] > 9.5 and R_T[i,0] <= 10.0:
                BINT[16,i] = R_T[i,1]
                BINR[16,i] = R_T[i,0]
            elif R_T[i,0] > 10.0 and R_T[i,0] <= 10.5:
                BINT[17,i] = R_T[i,1]
                BINR[17,i] = R_T[i,0]
            elif R_T[i,0] > 10.5 and R_T[i,0] <= 11.0:
                BINT[18,i] = R_T[i,1]
                BINR[18,i] = R_T[i,0]
            elif R_T[i,0] > 11.0 and R_T[i,0] <= 11.5:
                BINT[19,i] = R_T[i,1]
                BINR[19,i] = R_T[i,0]  
            elif R_T[i,0] > 11.5 and R_T[i,0] <= 12.0:
                BINT[20,i] = R_T[i,1]
                BINR[20,i] = R_T[i,0]
            elif R_T[i,0] > 12.0:
                BINT[21,i] = R_T[i,1]
                BINR[21,i] = R_T[i,0]

        AVG_T = np.zeros(BINS)
        AVG_R = np.zeros(BINS)
        VAR_T = np.zeros(BINS)
        
        for i in range(BINS):
            AVG_T[i] = ARRAY_AVG(BINT[i])
            AVG_R[i] = ARRAY_AVG(BINR[i])
            VAR_T[i] = ARRAY_STD(BINT[i])
        if n == 1:
            R_Stand = AVG_R
            T_Avg = AVG_T
            Var = VAR_T

            Var_Avg = np.mean(Var)
            #print Var_Avg
            #plt.figure()
            #plt.errorbar(R_Stand, T_Avg, yerr = Var, fmt = 'o') 
            #plt.title('Mean as a Function of Radius (1 Cluster, Variance as Errorbars)', fontsize=19)
            #plt.xlabel(r'$\mathrm{Radial}\/\mathrm{Distance}\/\mathrm{(Arcmin)}$', fontsize=19)
            #plt.ylabel(r'$\mathrm{Comptonization,}\/\mathrm{y}\/\mathrm{(\mu K)}$',fontsize=19)

            #plt.figure()
            #plt.scatter(R_Stand, T_Avg)
            #plt.title('Comptonization Factor as a Function of Radius (1 Clusters)', fontsize=19)
            #plt.xlabel(r'$\mathrm{Radial}\/\mathrm{Distance}\/\mathrm{(Arcmin)}$', fontsize=19)
            #plt.ylabel(r'$\mathrm{Comptonization,}\/\mathrm{y}\/\mathrm{(\mu K)}$',fontsize=19)

            #plt.figure()
            #plt.scatter(R_Stand, T_Avg/Var)
            #plt.title('Signal Over Variance as a Function of Radius (1 Cluster)', fontsize=19)
            #plt.xlabel(r'$\mathrm{Radial}\/\mathrm{Distance}\/\mathrm{(Arcmin)}$', fontsize=19)
            #plt.ylabel(r'$\mathrm{Comptonization,}\/\mathrm{y}\/\mathrm{(\mu K)}$',fontsize=19)
        elif n > 1:
            R_Stand = AVG_R
            T_Avg = (AVG_T) / (n)
            Var = (VAR_T) /(np.sqrt(n))

            Var_Avg = np.mean(Var)
            #print Var_Avg
            #plt.figure()
            #plt.errorbar(R_Stand, T_Avg, yerr = Var, fmt = 'o') 
            #plt.title('Mean as a Function of Radius (10 Clusters, Variance as Errorbars)', fontsize=19)
            #plt.xlabel(r'$\mathrm{Radial}\/\mathrm{Distance}\/\mathrm{(Arcmin)}$', fontsize=19)
            #plt.ylabel(r'$\mathrm{Comptonization,}\/\mathrm{y}\/\mathrm{(\mu K)}$',fontsize=19)

            #plt.figure()
            #plt.scatter(R_Stand, T_Avg)
            #plt.title('Comptonization Factor as a Function of Radius (1 Clusters)', fontsize=19)
            #plt.xlabel(r'$\mathrm{Radial}\/\mathrm{Distance}\/\mathrm{(Arcmin)}$', fontsize=19)
            #plt.ylabel(r'$\mathrm{Comptonization,}\/\mathrm{y}\/\mathrm{(\mu K)}$',fontsize=19)
    

            #plt.figure()
            #plt.scatter(R_Stand, T_Avg/Var)
            #plt.title('Signal Over Variance as a Function of Radius (10 Clusters)', fontsize=19)
            #plt.xlabel(r'$\mathrm{Radial}\/\mathrm{Distance}\/\mathrm{(Arcmin)}$', fontsize=19)
            #plt.ylabel(r'$\mathrm{Comptonization,}\/\mathrm{y}\/\mathrm{(\mu K)}$',fontsize=19)
    return Var_Avg #R_Stand, T_Avg, Var

#With SN_ALL_BIN we can run multiple simulations and then plot them together
#n is the number of simulations and c is the number of clusters in each sismulation
#Bins must be factors of 456, bins specify the number of data points we will plot
def STAT_UNCER_ALL(num_sim,num_clus):
    R = STAT_UNCERT(num_clus)[0]
    DATA = np.zeros((num_sim, len(R)))
    for i in range(num_sim):
        DATA[i] = STAT_UNCERT(num_clus)[1]

    print DATA
    MEAN = DATA.mean(axis=0)
    VAR = DATA.std(axis=0)
    print VAR
    VAR_Avg = np.mean(VAR)
    
    #plt.figure()
    #plt.scatter(R, MEAN)
    #plt.title('Comptonization as a Function of Radius (5 Simulations of 5 Clusters Each)', fontsize=19)
    #plt.xlabel(r'$\mathrm{Radial}\/\mathrm{Distance}\/\mathrm{(Arcmin)}$', fontsize=19)
    #plt.ylabel(r'$\mathrm{Comptonization,}\/\mathrm{y}\/\mathrm{(\mu K)}$',fontsize=19)
    
    plt.figure()
    plt.errorbar(R, MEAN, yerr = VAR, fmt = 'o') 
    plt.title('Comptonization as a Function of Radius (10 Simulations of 15 Clusters Each, Variance as Errorbars)', fontsize=19)
    plt.xlabel(r'$\mathrm{Radial}\/\mathrm{Distance}\/\mathrm{(Arcmin)}$', fontsize=19)
    plt.ylabel(r'$\mathrm{Comptonization,}\/\mathrm{y}\/\mathrm{(\mu K)}$',fontsize=19)
    
    plt.figure()
    plt.scatter(R, (MEAN)/(VAR))
    plt.title('Mean Over Variance as a Function of Radius (10 Simulations of 15 Clusters Each)', fontsize=19)
    plt.xlabel(r'$\mathrm{Radial}\/\mathrm{Distance}\/\mathrm{(Arcmin)}$', fontsize=19)
    plt.ylabel(r'$\mathrm{Comptonization,}\/\mathrm{y}\/\mathrm{(\mu K)}$',fontsize=19)
     
    return 


