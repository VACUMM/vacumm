# -*- coding: utf8 -*-
# Copyright or © or Copr. Actimar/IFREMER (2010-2015)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
# Classe valid: Validations 1D/2D/3D
#**********************************************************************************
# CRE : G. Charria (Ifremer)
#
# VER : 1.0 (04/2010)
#
#**********************************************************************************
#from valid import *
from vacumm.validator.valid import *

class ValidXYT(Valid):
    """ Validation spatiale XY """
    def __init__(self,model,obs):
        Valid.__init__(self,model,obs)

    def spatial_mean(self):
        """ Moyenne spatiale en chaque pas de temps pour les champs modelises et observes """
        import cdutil,  cdms2,  MV2,  numpy
        from vacumm.misc.grid import meshweights

        #self.modelspatialmean = MV2.average(MV2.average(self.model, axis=-1), axis=-1)
        # --> equivalent a nanmean(nanmean(sst')) sous matlab
        #self.modelspatialmean = MV2.average(MV2.average(self.model, axis=1), axis=1)
        # --> equivalent a nanmean(nanmean(sst)) sous matlab
        #print self.modelspatialmean


        self.model.spa_mean = list() #Initialise une liste
        self.obs.spa_mean = list() #Initialise une liste

        # Estimation de la moyenne a chaque pas de temps
        tps = self.obs.getTime()
        for i,  t in enumerate(tps):
            if self.obs[i, :, :].count()==0:
                #self.obs.spa_mean.append(self.obs.fill_value) # voir pour mettre un masque dans la variable a tracer ... ou pour choisir une autre valeur
                self.obs.spa_mean.append(0)
                self.model.spa_mean.append(0)
            else:
                #modif J.Gatti- pour que moyenne spatiale fontionne avec une grille curvilineaire et poids soient en metres

                lonmod=self.model.getLongitude()
                latmod=self.model.getLatitude()
                # lon et lat doivent etre 2D pour que meshweigths calcule des poids en m2
                if lonmod.rank()==1:
                    lonmod,latmod=numpy.meshgrid(lonmod,latmod)

                wtsmod=meshweights(lonmod, latmod, proj=True)
                # la matrice de poids doit avoir le meme type que la variable moyennee (MV dans notre cas)
                wtsmod = MV2.masked_array(wtsmod, MV2.getmask(self.model[i, :, :]), axes=self.model[i, :, :].getAxisList())

                self.model.spa_mean.append(cdutil.averager(self.model[i, :, :], axis='yx', weights=wtsmod))

                lonobs=self.obs.getLongitude()
                latobs=self.obs.getLatitude()

                if lonobs.rank()==1:
                    lonobs,latobs=numpy.meshgrid(lonobs,latobs)

                wtsobs=meshweights(lonobs, latobs, proj=True)
                wtsobs = MV2.masked_array(wtsobs, MV2.getmask(self.obs[i, :, :]), axes=self.obs[i, :, :].getAxisList())

                self.obs.spa_mean.append(cdutil.averager(self.obs[i, :, :], axis='yx', weights=wtsobs))

                #self.model.spa_mean.append(cdutil.averager(self.model[i, :, :],axis='yx'))
                #self.obs.spa_mean.append(cdutil.averager(self.obs[i, :, :],axis='yx'))
                #fin modif J.Gatti


        # Transformation en variable cdms2
        self.model.spa_mean = cdms2.createVariable(self.model.spa_mean,typecode='f',id='model_spa_mean', attributes=dict(units=self.model.units))
        self.obs.spa_mean = cdms2.createVariable(self.obs.spa_mean,typecode='f',id='obs_spa_mean', attributes=dict(units=self.obs.units))


        # Ajout de l'axe des temps correspondant a self...
        self.model.spa_mean.setAxis(0, tps)
        self.obs.spa_mean.setAxis(0, tps)

        # Les deux methodes suivantes donne le meme resultat mais la seconde est plus longue:
        #toto=cdutil.area_weights(self.model)
        #print toto
        #test=cdutil.averager(self.model, axis = 'yx', weights = ['weighted', 'weighted'], combinewts=1)
        #test=cdutil.averager(self.model, axis = 'xy')

    def temporal_mean(self):
        """ Moyennes temporelles pour les champs modelises et observes """
        import MV2, cdutil
        import numpy as np

        # --------- !!!!!!!  -------- Refelchir a Mettre plutot des return pour les fonctions ------- !!!!!!! ---------
        if np.ndim(self.model) < 3 :
            self.model.temp_mean = self.model
            self.obs.temp_mean = self.obs
        else :
            # Calcul de la moyenne temporelle ... Resultat => Map
            self.model.temp_mean = MV2.average(self.model, axis=0)
            self.obs.temp_mean = MV2.average(self.obs, axis=0)

    def bias_extrema(self):
        """ Carte des biais moyens extremes (min,max) entre les champs modelises et observes """
        import numpy as np
        import MV2
        # modif J.Gatti
        biais=self.obs-self.model

        # index des min et max
        imin=np.ma.argmin(abs(biais), axis=0)
        imax=np.ma.argmax(abs(biais), axis=0)

        #signe des min et max
        biasminsign=np.zeros((biais.shape[1], biais.shape[2]))
        biasmaxsign=np.zeros((biais.shape[1], biais.shape[2]))
        ind1=range(biais.shape[1])
        ind2=range(biais.shape[2])
        for i in ind1:
            for j in ind2:
                biasminsign[i, j]=np.sign(biais.data[imin[i, j], i, j])
                biasmaxsign[i, j]=np.sign(biais.data[imax[i, j], i, j])

        # min et max
        self.biasmin=MV2.min(abs(biais), axis=0) *biasminsign
        self.biasmax=MV2.max(abs(biais), axis=0)*biasmaxsign

        del biais
        #fin modif
# cartes des biais a la date ou le biais spatial moyen est minimum/maximum
#        self.spatial_mean()
#
#        biais = self.obs.spa_mean.data-self.model.spa_mean.data
#        imin = np.nonzero(abs(biais) == abs(biais).min())
#        imax = np.nonzero(abs(biais) == abs(biais).max())
#
#        self.biasmin = self.obs[imin[0][0], :, :]-self.model[imin[0][0], :, :]
#        self.biasmax = self.obs[imax[0][0], :, :]-self.model[imax[0][0],  :, :]

    def systematic_bias(self):
        """ Biais systematique (+1: sous-estimation du modele, -1: sur-estimation du modele) """
        import numpy as np
        import cdms2
        from vacumm.misc.grid import get_grid,  set_grid

        map_bias = list() #Initialise une liste

        # Estimation de la moyenne a chaque pas de temps
        tps = self.obs.getTime()
        for i,  t in enumerate(tps):
            map_bias.append(np.sign(self.obs[i, :, :]-self.model[i, :, :] ))

        map_bias = cdms2.createVariable(map_bias,typecode='f',id='computation')

        res = np.sum(map_bias,  axis=0)
        res
        # Trouve valeurs egale a la longueur de la serie temporelle (<=> uniquement valeurs positives)
        one = (res==len(map_bias)).nonzero()

        # Trouve valeurs egales a moins la longueur de la serie temporelle (<=> uniquement valeurs negatives)
        mone = (res==-len(map_bias)).nonzero()

        self.biassyst = np.zeros(res.shape)
        self.biassyst[one]=1.
        self.biassyst[mone]=-1.

        self.biassyst = cdms2.createVariable(self.biassyst, typecode='f',id='syst_bias')
        ggm = get_grid(self.model)
        set_grid(self.biassyst, ggm, axes=True)




    def temporal_std(self):
        """ Ecart-type en chaque point geographique """
        from genutil import statistics
        from vacumm.misc.grid import get_grid,  set_grid
        # centered and biased std (cf. http://www2-pcmdi.llnl.gov/cdat/manuals/cdutil/cdat_utilities-2.html)
        self.model.temp_std = () #Initialise un tuple
        self.obs.temp_std = () #Initialise un tuple

        self.model.temp_std = statistics.std(self.model, axis = 0)
        self.obs.temp_std = statistics.std(self.obs, axis = 0)

        ggm = get_grid(self.model)
        set_grid(self.model.temp_std, ggm, axes=True)

        ggo = get_grid(self.obs)
        set_grid(self.obs.temp_std, ggm, axes=True)


    def spatial_std(self):
        """ Ecart-type pour chaque pas de temps """
        from genutil import statistics
        import cdms2
        # centered and biased std (cf. http://www2-pcmdi.llnl.gov/cdat/manuals/cdutil/cdat_utilities-2.html)
        self.model.spa_std = list() #Initialise une liste
        self.obs.spa_std = list() #Initialise une liste

        # Estimation de la std a chaque pas de temps
        tps = self.obs.getTime()
        for i,  t in enumerate(tps):
            self.model.spa_std.append(statistics.std(self.model[i, :, :],axis='yx'))
            self.obs.spa_std.append(statistics.std(self.obs[i, :, :],axis='yx'))

        # Transformation en variable cdms2
        self.model.spa_std = cdms2.createVariable(self.model.spa_std,typecode='f',id='model_spa_std', attributes=dict(units=self.model.units))
        self.obs.spa_std = cdms2.createVariable(self.obs.spa_std,typecode='f',id='obs_spa_std', attributes=dict(units=self.obs.units))


        # Ajout de l'axe des temps correspondant a self...
        self.model.spa_std.setAxis(0, tps)
        self.obs.spa_std.setAxis(0, tps)



    def extrema(self):
        """ Extraction des valeurs minimum et maximum dans les champs modelises et observes """
        from genutil import minmax

        self.model.extrema=minmax(self.model)
        self.obs.extrema=minmax(self.obs)

        # Affichage du resultat
        m=str(self.model.extrema)
        o=str(self.obs.extrema)

        #print(" => Minimum et maximum simules: "+ m)
        #print(" => Minimum et maximum observes: "+ o)

    def temporal_rms(self):
        """ RMS entre modele et observations - calculee en chaque point sur la dimension
            temporelle / Resultat => Carte de RMS uncentered and biased """
        from genutil import statistics
        from vacumm.misc.grid import get_grid,  set_grid


        self.temp_rms = () #Initialise un tuple

        self.temp_rms = statistics.rms(self.model, self.obs, axis = 0)


        gg = get_grid(self.model)
        set_grid(self.temp_rms, gg, axes=True)

    def temporal_rmsc(self):
        """ RMS entre modele et observations - calculee en chaque point sur la dimension
            temporelle / Resultat => Carte de RMS centered and biased """
        from genutil import statistics
        from vacumm.misc.grid import get_grid,  set_grid


        self.temp_rmsc = () #Initialise un tuple

        self.temp_rmsc = statistics.rms(self.model, self.obs, axis = 0,  centered = 1)


        gg = get_grid(self.model)
        set_grid(self.temp_rmsc, gg, axes=True)

    def spatial_rms(self):
        """ RMS pour chaque pas de temps """
        from genutil import statistics
        import cdms2
        import numpy as np

        # centered and biased RMS
        self.spa_rms = list() #Initialise une liste

        # Estimation de la std a chaque pas de temps
        tps = self.obs.getTime()
        for i,  t in enumerate(tps):
            minterm = np.reshape(self.model[i, :, :], (self.model.shape[-1]*self.model.shape[-2]))
            ointerm = np.reshape(self.obs[i, :, :], (self.obs.shape[-1]*self.obs.shape[-2]))
            #modif J.Gatti - pour MFS
            if not minterm._grid_  is None:
                minterm._grid_=None
            if not ointerm._grid_  is None:
                ointerm._grid_=None
            #fin modif J.Gatti
            self.spa_rms.append(statistics.rms(minterm, ointerm))

        # Transformation en variable cdms2
        self.spa_rms = cdms2.createVariable(self.spa_rms,typecode='f',id='spa_rms', attributes=dict(units=self.model.units))


        # Ajout de l'axe des temps correspondant a self...
        self.spa_rms.setAxis(0, tps)

    def spatial_rmsc(self):
        """ RMS centree pour chaque pas de temps """
        from genutil import statistics
        import cdms2
        import numpy as np
        # centered and biased RMS
        self.spa_rmsc = list() #Initialise une liste

        # Estimation de la std a chaque pas de temps
        tps = self.obs.getTime()
        for i,  t in enumerate(tps):
            minterm = np.reshape(self.model[i, :, :], (self.model.shape[-1]*self.model.shape[-2]))
            ointerm = np.reshape(self.obs[i, :, :], (self.obs.shape[-1]*self.obs.shape[-2]))
            #modif J.Gatti - pour MFS
            if not minterm._grid_  is None:
                minterm._grid_=None
            if not ointerm._grid_  is None:
                ointerm._grid_=None
            #fin modif J.Gatti
            self.spa_rmsc.append(statistics.rms(minterm, ointerm, centered=1))

        # Transformation en variable cdms2
        self.spa_rmsc = cdms2.createVariable(self.spa_rmsc,typecode='f',id='spa_rmsc', attributes=dict(units=self.model.units))


        # Ajout de l'axe des temps correspondant a self...
        self.spa_rmsc.setAxis(0, tps)

    def temporal_corr(self):
        """ Correlation entre modele et observations - calculee en chaque point sur la dimension
            temporelle / Resultat => Carte de Correlation uncentered and biased """
        from genutil import statistics
        from vacumm.misc.grid import get_grid,  set_grid


        self.temp_corr = () #Initialise un tuple

        self.temp_corr = statistics.correlation(self.model, self.obs, axis = 0)

        gg = get_grid(self.model)
        set_grid(self.temp_corr, gg, axes=True)


    def obs_coverage(self):
        """ Calcul du taux de couverture des observations en chaque point sur la periode """
        import numpy as np
        from vacumm.misc.grid import get_grid,  set_grid
        import cdms2

        self.obs_cov = () #Initialise un tuple

        self.obs_cov = cdms2.createVariable(self.obs.count(axis=0),typecode='f',id='obs_cov', attributes=dict(units='%'))
        self.obs_cov = self.obs_cov / len(self.obs) * 100.

        gg = get_grid(self.obs)
        set_grid(self.obs_cov, gg, axes=True)

    def spatial_obs_coverage(self):
        """ Calcul du taux de couverture spatiale des observations pour chaque pas de temps """
        import numpy as np
        from vacumm.misc.grid import get_grid,  set_grid
        import cdms2

        self.obs_spacov = list() #Initialise une liste

        # Estimation du taux de couverture a chaque pas de temps
        tps = self.obs.getTime()
        for i,  t in enumerate(tps):
            self.obs_spacov.append(self.obs[i, :, :].count(axis=None))

        # Transformation en variable cdms2
        self.obs_spacov = cdms2.createVariable(self.obs_spacov,typecode='f',id='obs_spacov', attributes=dict(units='%'))

        ## Passage en pourcentages
        #self.obs_spacov = self.obs_spacov / (self.obs.shape[1]*self.obs.shape[2]) * 100.

        # Ajout de l'axe des temps correspondant a self.obs
        self.obs_spacov.setAxis(0, tps)







