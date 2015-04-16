import matplotlib.pyplot as plt
import numpy as np
from pylab import *

#Here's a whole bunch of constants
H0=72. #km / (s*Mpc)
OmegaM= 0.27 #matter density
hubble = H0/100.
hubble70 = H0/70.
mmw = 0.59 #mean molecular weight of gas
mu_e = 1.143 #mean molecular weight of electrons in cluster
mpc = 3.08568025e22   #megaparsec in meters
Msol = 1.98892e30 #solar masses in kg
#Arnaud profile constant
P0 = 8.403/ (hubble70)^1.5
alpha = 1.0510
beta = 5.4905
gamma = 0.3081
index = (beta - gamma) / alpha
c500 = 1.177




P_0 = 8.93 

c500 = 1.33 
#used that from arnaud code

h_0 = 0.675
#used h_0 as mean value bwt 0.5 and 0.85

rho_crit = (1.878e-29)*((h_0)**2)
#used rho_crit as 1.878e-29 g/cm^3, could have used 1.0553e-5 GeV/cm^3
#also in later models, rho-crit is function of the hubble constant which is
#a function o red shift

R500 = (0.879)#*(3.08568025e22)
#value 1.010 was choosen randomly from arnaud (2010) cluster RXC J1516+0005

M500 =(((4.)*(np.pi))/(3.))*(((R500)**3)*(500)*(rho_crit))
#used equation from arnaud paper(2010)

x = np.arange(1,10,0.1) / (R500)
y = np.zeros(90)
for i in range(len(x)):
    y[i]=((P_0)/(((c500*x[i])**gamma)*(1 + (c500*x[i])**alpha)**index))/ (1.466)#e-3)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
line = ax.plot(y, color='blue')
ax.set_xscale('log')
ax.set_yscale('log')
plt.plot(x,y)

plt.show()

