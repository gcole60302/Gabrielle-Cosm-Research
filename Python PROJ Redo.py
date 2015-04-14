import matplotlib.pyplot as plt
import numpy as np


alpha = 1.0510 
beta = 5.4905 
gamma = 0.3081 
index = (beta - gamma) / alpha
#all four used/calculated from arnaud paper(2010)

h_70 = (72./70.)
#used that from arnaud code

P_0 = 8.403*(h_70)**(-3./2.)
#used that from arnaud paper(2010)

c500 = 1.177 
#used that from arnaud code

h_0 = 0.675
#used h_0 as mean value bwt 0.5 and 0.85

rho_crit = (1.8782e-29)*((h_0)**2)
#used rho_crit as 1.878e-29 g/cm^3, could have used 1.0553e-5 GeV/cm^3
#also in later models, rho-crit is function of the hubble constant which is
#a function of red shift

R500 = (1.010)*(3.08568025e22)
#value 1.010 was choosen randomly from arnaud (2010) cluster RXC J1516+0005

M500 =(((4.)*(np.pi))/(3.))*(((R500)**3)*(500)*(rho_crit))
#used equation from arnaud paper(2010)


x = np.arange(1,100) / (R500)
y = np.zeros(99)
for i in range(len(x)):
    y[i]=(P_0)/(((c500*x[i])**gamma)*(1 + (c500*x[i])**alpha)**index)

plt.plot(x,y)
plt.show()

