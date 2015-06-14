    Radial_1 = np.zeros((n, 7900))
    Temp_1 = np.zeros((n, 7900))
    Radial_2 = np.zeros((n, 7900))
    Temp_2 = np.zeros((n, 7900))
    Radial_3 = np.zeros((n, 7900))
    Temp_3 = np.zeros((n, 7900))
    Radial_4 = np.zeros((n, 7900))
    Temp_4 = np.zeros((n, 7900))
    Radial_5 = np.zeros((n, 7900))
    Temp_5 = np.zeros((n, 7900))
####################
    for i in range(n):
        if Cluster[i,3] >= 5e14:
            a = len(PROFILER(Cluster[i,2], Cluster[i,3]))
            Radial_1[i,0:a] = PROFILER(Cluster[i,2], Cluster[i,3])
            Temp_1[i,0:a] = PROFILET(Cluster[i,2], Cluster[i,3])
        if Cluster[i,3] <= 5e14 and  Cluster[i,3] >= 4e14:
            a = len(PROFILER(Cluster[i,2], Cluster[i,3]))
            Radial_2[i,0:a] = PROFILER(Cluster[i,2], Cluster[i,3])
            Temp_2[i,0:a] = PROFILET(Cluster[i,2], Cluster[i,3])
        if Cluster[i,3] <= 4e14 and  Cluster[i,3] >= 3e14:
            a = len(PROFILER(Cluster[i,2], Cluster[i,3]))
            Radial_3[i,0:a] = PROFILER(Cluster[i,2], Cluster[i,3])
            Temp_3[i,0:a] = PROFILET(Cluster[i,2], Cluster[i,3])
        if Cluster[i,3] <= 3e14 and  Cluster[i,3] >= 2e14:
            a = len(PROFILER(Cluster[i,2], Cluster[i,3]))
            Radial_4[i,0:a] = PROFILER(Cluster[i,2], Cluster[i,3])
            Temp_4[i,0:a] = PROFILET(Cluster[i,2], Cluster[i,3])
        if Cluster[i,3] <= 2e14:
            a = len(PROFILER(Cluster[i,2], Cluster[i,3]))
            Radial_5[i,0:a] = PROFILER(Cluster[i,2], Cluster[i,3])
            Temp_5[i,0:a] = PROFILET(Cluster[i,2], Cluster[i,3])

#############################
    T_R_BIN1 = np.trim_zeros(Radial_1.max(axis=0), 'b')
    T_T_BIN1 = np.zeros((n, len(T_R_BIN1)))

    T_R_BIN2 = np.trim_zeros(Radial_2.max(axis=0), 'b')
    T_T_BIN2 = np.zeros((n, len(T_R_BIN2)))

    T_R_BIN3 = np.trim_zeros(Radial_3.max(axis=0), 'b')
    T_T_BIN3 = np.zeros((n, len(T_R_BIN3)))

    T_R_BIN4 = np.trim_zeros(Radial_4.max(axis=0), 'b')
    T_T_BIN4 = np.zeros((n, len(T_R_BIN4)))

    T_R_BIN5 = np.trim_zeros(Radial_5.max(axis=0), 'b')
    T_T_BIN5 = np.zeros((n, len(T_R_BIN5)))
################################# 
    for i in range(n):
        if Cluster[i,3] >= 5e14:
            interpol = scipy.interpolate.interp1d(PROFILER(Cluster[i,2], Cluster[i,3]), PROFILET(Cluster[i,2], Cluster[i,3]),kind='cubic', bounds_error=False, fill_value=0.)
            T_T_BIN1[i,1:] = interpol(T_R_BIN1[1:])
        if Cluster[i,3] <= 5e14 and  Cluster[i,3] >= 4e14:
            interpol = scipy.interpolate.interp1d(PROFILER(Cluster[i,2], Cluster[i,3]), PROFILET(Cluster[i,2], Cluster[i,3]),kind='cubic', bounds_error=False, fill_value=0.)
            T_T_BIN2[i,1:] = interpol(T_R_BIN2[1:])
        if Cluster[i,3] <= 4e14 and  Cluster[i,3] >= 3e14:
            interpol = scipy.interpolate.interp1d(PROFILER(Cluster[i,2], Cluster[i,3]), PROFILET(Cluster[i,2], Cluster[i,3]),kind='cubic', bounds_error=False, fill_value=0.)
            T_T_BIN3[i,1:] = interpol(T_R_BIN3[1:])
        if Cluster[i,3] <= 3e14 and  Cluster[i,3] >= 2e14:
            interpol = scipy.interpolate.interp1d(PROFILER(Cluster[i,2], Cluster[i,3]), PROFILET(Cluster[i,2], Cluster[i,3]),kind='cubic', bounds_error=False, fill_value=0.)
            T_T_BIN4[i,1:] = interpol(T_R_BIN4[1:])
        if Cluster[i,3] <= 2e14:
            interpol = scipy.interpolate.interp1d(PROFILER(Cluster[i,2], Cluster[i,3]), PROFILET(Cluster[i,2], Cluster[i,3]),kind='cubic', bounds_error=False, fill_value=0.)
            T_T_BIN5[i,1:] = interpol(T_R_BIN5[1:])
##################################
    T_T_BIN1 = T_T_BIN1.sum(axis=0)
    T_T_BIN2 = T_T_BIN2.sum(axis=0)
    T_T_BIN3 = T_T_BIN3.sum(axis=0)
    T_T_BIN4 = T_T_BIN4.sum(axis=0)
    T_T_BIN5 = T_T_BIN5.sum(axis=0)
#################################
    
    F_BIN1 = scipy.interpolate.interp1d(T_R_BIN1, T_T_BIN1, kind= 'cubic', bounds_error=False, fill_value=0.)
    F_BIN2 = scipy.interpolate.interp1d(T_R_BIN2, T_T_BIN2, kind= 'cubic', bounds_error=False, fill_value=0.)
    F_BIN3 = scipy.interpolate.interp1d(T_R_BIN3, T_T_BIN3, kind= 'cubic', bounds_error=False, fill_value=0.)
    F_BIN4 = scipy.interpolate.interp1d(T_R_BIN4, T_T_BIN4, kind= 'cubic', bounds_error=False, fill_value=0.)
    F_BIN5 = scipy.interpolate.interp1d(T_R_BIN5, T_T_BIN5, kind= 'cubic', bounds_error=False, fill_value=0.)

    SN_ALL_BIN = np.zeros((5, len(AVG_R[AVG_R!=0])))
    R_Stand = AVG_R[AVG_R!=0]

#SN_ALL_BIN is an array of the S/N ratio of the clusters in a specific bin
    for i in range(len(R_Stand)):
        SN_ALL_BIN[0,i] = (F_BIN1(R_Stand[i]))/ (np.sqrt((BIN11[i] - F_BIN1(R_Stand[i]))**2.))
        SN_ALL_BIN[1,i] = (F_BIN2(R_Stand[i]))/ (np.sqrt((BIN21[i] - F_BIN2(R_Stand[i]))**2.))
        SN_ALL_BIN[2,i] = (F_BIN3(R_Stand[i]))/ (np.sqrt((BIN31[i] - F_BIN3(R_Stand[i]))**2.))
        SN_ALL_BIN[3,i] = (F_BIN4(R_Stand[i]))/ (np.sqrt((BIN41[i] - F_BIN4(R_Stand[i]))**2.))
        SN_ALL_BIN[4,i] = (F_BIN5(R_Stand[i]))/ (np.sqrt((BIN51[i] - F_BIN5(R_Stand[i]))**2.))

    plt.figure()
    plt.title('Bin 1')
    plt.scatter(R_Stand,SN_ALL_BIN[0])

    plt.figure()
    plt.title('Bin 2')
    plt.scatter(R_Stand,SN_ALL_BIN[1])

    plt.figure()
    plt.title('Bin 3')
    plt.scatter(R_Stand,SN_ALL_BIN[2])

    plt.figure()
    plt.title('Bin 4')
    plt.scatter(R_Stand,SN_ALL_BIN[3])

    plt.figure()
    plt.title('Bin 5')
    plt.scatter(R_Stand,SN_ALL_BIN[4])
