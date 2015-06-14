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
plt.ion()
import numpy as np                                                                                                                                                                                                                                                             
import pylab as plt                                                                                                                                                                                                                                                            
import matplotlib.ticker as ticker                                                                                                                                                                                                                                             
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy import interpolate

F_BIN1 = scipy.interpolate.interp1d(T_R_BIN1, T_T_BIN1, kind='cubic')
    F_BIN2 = scipy.interpolate.interp1d(T_R_BIN2, T_T_BIN2, kind='cubic')
    F_BIN3 = scipy.interpolate.interp1d(T_R_BIN3, T_T_BIN3, kind='cubic')
    F_BIN4 = scipy.interpolate.interp1d(T_R_BIN4, T_T_BIN4, kind='cubic')
    F_BIN5 = scipy.interpolate.interp1d(T_R_BIN5, T_T_BIN5, kind='cubic')

    SN_BIN1 = np.zeros((5, len(AVG_R[AVG_R!=0])))
    R_Stand = AVG_R[AVG_R!=0]

    for i in range(len(R_Stand)):
        SN_BIN1[i] = (F_BIN1(R_Stand[i]))/ (np.sqrt(BIN11[i] - F_BIN1(R_Stand[i])))

    for i in range(len(R_Stand)):
        SN_BIN1[i] = (F_BIN2(R_Stand[i]))/ (np.sqrt(BIN21[i] - F_BIN1(R_Stand[i])))

    for i in range(len(R_Stand)):
        SN_BIN1[i] = (F_BIN2(R_Stand[i]))/ (np.sqrt(BIN31[i] - F_BIN1(R_Stand[i])))

