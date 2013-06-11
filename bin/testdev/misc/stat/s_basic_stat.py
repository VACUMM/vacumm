
# Module de calcul de statistiques de validation simples (RMS, ...)
#**********************************************************************************
# MOD : s_basic_stat.py
# OBJ : Statistiques
#
# CRE : G. Charria (Ifremer)
#       S. Theetten (Ifremer)
# VER : 1.0 (12/2009)
#**********************************************************************************
def spatial_average(lon,lat,data):
    # Function to compute the spatial average of data at each time step
    # -----------------------------------------------------------------------------
    import numpy
   
    res_timeserie=numpy.empty((len(data)))

    for itime in range(0, len(data)):
       res_timeserie[itime]=numpy.mean(data[itime])
    
    # Controle
    #print res_timeserie


    # Figure de controle
    #import matplotlib
    #import matplotlib.pyplot as plt
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(res_timeserie)
    #plt.title('Mean Obs')
    #plt.show()




    return res_timeserie
    #-- Fin du calcul
    #------------------------------------------------------------------------------
def spatial_std(lon,lat,data):
    # Function to compute the standard deviation of data at each time step
    # -----------------------------------------------------------------------------
    import numpy
   
    res_timeserie=numpy.empty((len(data)))

    for itime in range(0, len(data)):
       res_timeserie[itime]=numpy.std(data[itime])
    
    # Controle
    #print res_timeserie


    # Figure de controle
    #import matplotlib
    #import matplotlib.pyplot as plt
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(res_timeserie)
    #plt.title('Std Obs')
    #plt.show()




    return res_timeserie
    #-- Fin du calcul
    #------------------------------------------------------------------------------
def rmsdiff_time(lon,lat,data):
    # Function to compute the RMS at each time step
    # -----------------------------------------------------------------------------
    import numpy
   
    res_timeserie=numpy.empty((len(data)))

    for itime in range(0, len(data)):
       res_timeserie[itime]=numpy.sqrt(numpy.mean(numpy.multiply(data[itime],data[itime])))
    
    # Controle
    #print res_timeserie


    # Figure de controle
    import matplotlib
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(res_timeserie)
    plt.title('Std Obs')
    plt.show()




    return res_timeserie
    #-- Fin du calcul
    #------------------------------------------------------------------------------



