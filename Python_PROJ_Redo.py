import matplotlib.pyplot as plt
import numpy as np
from pylab import *

#All values choosen from cluster RXC J0003.8+0203
alpha =1.41
beta = 5.4905 
gamma =0.567
index = (beta - gamma) / alpha
#all four used/calculated from arnaud paper(2010)

h_70 = (72./70.)
#used that from arnaud code

P_0 = 3.93 #8.403*(h_70)**(-3./2.): from arnaud paper(2010)

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

q=10e3
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

