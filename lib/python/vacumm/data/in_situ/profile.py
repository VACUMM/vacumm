# -*- coding: utf8 -*-
# Copyright or © or Copr. Actimar/IFREMER (2013-2015)
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

class Profile() :

    glider_begin = "EXGL"

    glider_type = "GLIDER"
    recopeca_type = "RECOPESCA"

    #cle des differentes tables
    time_key = "time"
    lon_key = "lon"
    lat_key = "lat"
    depth_key = "depth"
    temp_key = "temp"
    sal_key = "sal"


    """ Profile vertical """
    def __init__(self, time, platform_name, lon, lat):

        self.name = "Vertical Profile"
        self.shortname = "VertProf"

        self.time = time
        self.lon = lon
        self.lat = lat

        #tableau numpy
        self.depth = []
        self.temp = []
        self.sal = []
        #...
        # AJOUTER NOUVELLES VARIABLES ICI
        #...

        self.platform_name = platform_name
        if self.platform_name.startswith(self.glider_begin):
            self.platform_type = self.glider_type
        else:
            self.platform_type = self.recopeca_type


    # initialisation : creation des list du profil
    def add_time_step(self, depth=None, temp=None, sal=None):
        self.depth.append(depth)
        self.temp.append(temp)
        self.sal.append(sal)
        #...
        # AJOUTER NOUVELLES VARIABLES ICI
        #...


    #fin initialistion :creation table numpy
    def convert_tables(self):
        from numpy import array
        self.depth = array(self.depth, dtype='float')
        self.temp = array(self.temp, dtype='float')
        self.sal = array(self.sal, dtype='float')
        #...
        # AJOUTER NOUVELLES VARIABLES ICI
        #...


    # creation de tables pour tracer les profils
    # retourne une map de ndarray
    def create_tabs(self):
        from numpy import ones_like
        from matplotlib.dates import date2num

        #duplication sur le nombre de profondeur (utile pour tracer le profil)
        time_array = ones_like(self.depth) * date2num(self.time)
        lon_array = ones_like(self.depth) * float(self.lon)
        lat_array = ones_like(self.depth) * float(self.lat)

        map_return = {self.time_key : time_array,
                      self.lon_key :lon_array,
                      self.lat_key :lat_array
                      }

        return map_return


    # retourne les infos date,lon et lat du profil
    def get_position(self):
        return [self.time, self.lon, self.lat]

