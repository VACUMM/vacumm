# -*- coding: utf8 -*-
"""Conventions about data formats and names"""
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
#

from warnings import warn
from collections import OrderedDict
import string

import cdms2, MV2, re
from vacumm import VACUMMError
from vacumm.misc import kwfilter, dict_merge
from vacumm.misc.axes import create as create_axis, isaxis
from vacumm.misc.grid import create_axes2d
from vacumm.misc.io import ncmatch_obj
from vacumm.data.misc.arakawa import ARAKAWA_LOCATIONS

__all__ = ['VAR_SPECS', 'AXIS_SPECS',
    'GENERIC_AXIS_NAMES', 'GENERIC_VAR_NAMES', 'GENERIC_NAMES',
    'format_var', 'format_axis', 'format_grid', 'match_obj',
    'cf2atts', 'cf2search', 'cp_suffix', 'specs_def_loc', 'get_loc',
    'change_loc', 'change_loc_single', 'dupl_loc_specs', 'no_loc_single',
    'change_loc_specs', 'squeeze_loc_single', 'squeeze_loc', 'get_physloc',
    'HIDDEN_CF_ATTS', 'set_loc'
]

ARAKAWA_SUFFIXES = [('_'+p) for p in ARAKAWA_LOCATIONS]
arakawa_suffixes = ARAKAWA_SUFFIXES # compat
arakawa_locations = ARAKAWA_LOCATIONS

#: Specifications for variables
VAR_SPECS = OrderedDict(

    # Thermodynamics
    temp = dict(
        names = ['temp', 'temperature', 'TEMP'],
        standard_names = ['sea_water_temperature', 'sea_water_potential_temperature'],
        long_names = 'Temperature',
        units =  'degrees_celsius',
#        axes = {'t':'time', 'x':'lon', 'y':'lat'},
    ),
    ptemp = dict(
        names = ['ptemp'],
        standard_names = ['sea_water_potential_temperature'],
        long_names = 'Potential temperature',
        units =  'degrees_celsius',
    ),
    sal = dict(
        names=['sal', 'psal', 'salinity', 'SAL'],
        standard_names='sea_water_salinity',
        long_names = 'Salinity',
        units = 'PSU',
    ),
    sst = dict(
        names = ['sst'],
        standard_names = 'sea_surface_temperature',
        long_names = 'Sea surface temperature',
        units = 'degrees_celsius',
    ),
    sss = dict(
        names = ['sss'],
        standard_names = 'sea_surface_salinity',
        long_names = 'Sea surface salinity',
        units = 'PSU',
    ),
    dens = dict(
        names = ['dens'],
        standard_names = 'sea_water_density',
        long_names = 'Sea water density',
        units = 'kg m-3',
    ),
    sigmat = dict(
        names = ['sigmatheta'],
        standard_names = 'sea_water_sigma_t',
        long_names = 'Sea water density minus 1000',
        units = 'kg m-3',
    ),
    ndens = dict(
        names = ['ndens'],
        standard_names = 'sea_water_neutral_density',
        long_names = 'Sea water neutral density',
        units = 'kg m-3',
    ),
    sigmatheta = dict(
        names = ['sigmatheta'],
        standard_names = 'sea_water_sigma_theta',
        long_names = 'Sea water potential density minus 1000',
        units = 'kg m-3',
    ),
    pdens = dict(
        names = ['pdens', 'sigma0'],
        standard_names = 'sea_water_potential_density',
        long_names = 'Sea water potential density',
        units = 'kg m-3',
    ),
    sigma0 = dict(
        inherit = 'pdens'
    ),
    sigma1 = dict(
        names = ['sigma1'],
        standard_names = 'sea_water_potential_density',
        long_names = 'Sea water potential density with ref at 1000 dbar',
        units = 'kg m-3',
    ),
    sigma2 = dict(
        names = ['sigma2'],
        standard_names = 'sea_water_potential_density',
        long_names = 'Sea water potential density with ref at 2000 dbar',
        units = 'kg m-3',
    ),
    sigma3 = dict(
        names = ['sigma3'],
        standard_names = 'sea_water_potential_density',
        long_names = 'Sea water potential density with ref at 3000 dbar',
        units = 'kg m-3',
    ),
    sigma4 = dict(
        names = ['sigma4'],
        standard_names = 'sea_water_potential_density',
        long_names = 'Sea water potential density with ref at 4000 dbar',
        units = 'kg m-3',
    ),
    ssd = dict(
        names = ['ssd'],
        standard_names = 'sea_surface_density',
        long_names = 'Sea surface density',
        units = 'PSU',
    ),
    conduct = dict(
        standard_names = 'sea_water_electrical_conductivity',
        long_names = 'Sea water electrial conductivity',
        units = 'S m-1',
    ),
    sndspd = dict(
        standard_names = 'speed_of_sound_in_sea_water',
        long_names = 'Speed of sound in water',
        units = 'm s-1',
    ),
    mld = dict(
        standard_names = 'mixed_layer_depth',
        long_names = "Mixed layer depth",
        units = "m",
        physloc = 't',
    ),
    ped = dict(
        standard_names = 'potential_energy_deficit',
        long_names = "Potential energy deficit",
        units = "J m-2",
        physloc = 't',
    ),
    ohc = dict(
        standard_names = 'ocean_heat_content',
        long_names = 'Ocean heat content',
        units = 'J',
        physloc = 't',
    ),
    osc = dict(
        standard_names = 'ocean_salt_content',
        long_names = 'Ocean salt content',
        units = 'kg',
        physloc = 't',
    ),
    cp = dict(
        standard_names = 'specific_heat_capacity',
        long_names = 'Specific heat capacity',
        units = 'J K-1',
        physloc = 't',
    ),



    # Dynamics
    ssh = dict(
        names=['ssh', 'xe'],
        standard_names=['sea_surface_height_above_sea_level','sea_surface_height_above_geoid' ],
        long_names = 'Sea surface height',
        units = 'm',
    ),
    u3d = dict(
        names=['uz', 'u3d'],
        standard_names=['sea_water_x_velocity', #'sea_water_x_velocity_at_u_location',
            'eastward_sea_water_velocity'
            ],
        long_names = "Sea water velocity along X",
        units = "m s-1",
        axes = dict(x=['lon_u'], y=['lat_u']),
        atlocs = ['t', 'u', 'v'],
        physloc = 'u',
    ),
    v3d = dict(
        names=['vz', 'v3d'],
        standard_names=['sea_water_y_velocity', #'sea_water_y_velocity_at_v_location',
            'northward_sea_water_velocity'],
        long_names = "Sea water velocity along Y",
        units = "m s-1",
        axes = dict(x=['lon_v'], y=['lat_v']),
        atlocs = ['t', 'u', 'v'],
        physloc = 'v',
    ),
    u = dict(
        names=['uz', 'u3d'],
        standard_names=['sea_water_x_velocity', #'sea_water_x_velocity_at_u_location',
            'eastward_sea_water_velocity'
            ],
        long_names = "Sea water velocity along X",
        units = "m s-1",
        axes = dict(x=['lon_u'], y=['lat_u']),
        atlocs = ['t', 'u', 'v'],
        physloc = 'u',
    ),
    v = dict(
        names=['vz', 'v3d'],
        standard_names=['sea_water_y_velocity', #'sea_water_y_velocity_at_v_location',
            'northward_sea_water_velocity'],
        long_names = "Sea water velocity along Y",
        units = "m s-1",
        axes = dict(x=['lon_v'], y=['lat_v']),
        atlocs = ['t', 'u', 'v'],
        physloc = 'v',
    ),
    w3d = dict(
        names=['wz', 'w3d'],
        standard_names=['sea_water_z_velocity_at_w_location', 'sea_water_z_velocity'],
        long_names = "Sea water velocity along Z at W location",
        units = "m s-1",
        physloc = 'w',
    ),
    ubt = dict(
        names=['ubt', 'u2d', 'u'],
        standard_names=['barotropic_sea_water_x_velocity'],
        long_names = "Sea water barotropic velocity along X",
        units = "m s-1",
        axes = dict(x=['lon_u'], y=['lat_u']),
        atlocs = ['t', 'u', 'v'],
        physloc = 'u',
    ),
    vbt = dict(
        names=['vbt', 'v2d', 'v'],
        standard_names=['barotropic_sea_water_y_velocity'],
        long_names = "Sea water barotropic velocity along Y",
        units = "m s-1",
        axes = dict(x=['lon_v'], y=['lat_v']),
        atlocs = ['t', 'u', 'v'],
        physloc = 'v',
    ),
    ubc = dict(
        names=['ubc', 'u'],
        standard_names=['baroclinic_sea_water_x_velocity'],
        long_names = "Sea water baroclinic velocity along X",
        units = "m s-1",
        axes = dict(x=['lon_u'], y=['lat_u']),
        atlocs = ['t', 'u', 'v'],
        physloc = 'u',
    ),
    vbc = dict(
        names=['vbc', 'v'],
        standard_names=['baroclinic_sea_water_y_velocity'],
        long_names = "Sea water baroclinic velocity along Y",
        units = "m s-1",
        axes = dict(x=['lon_v'], y=['lat_v']),
        atlocs = ['t', 'u', 'v'],
        physloc = 'v',
    ),
    usurf = dict(
        names = ['usurf'],
        standard_names=['sea_surface_x_velocity'],
        long_names = ["Sea surface velocity along X"],
        units = "m s-1",
        axes = dict(y=['lat'], x=['lon']),
        atlocs = ['t', 'u', 'v'],
        physloc = 'u',
    ),
    vsurf = dict(
        names = ['vsurf'],
        standard_names=['sea_surface_y_velocity', 'sea_surface_y_velocity_at_y_location'],
        long_names = ["Sea surface velocity along Y", "Sea surface velocity along Y at V location"],
        units = "m s-1",
        axes = dict(x=['lon_v'], y=['lat_v']),
        atlocs = ['t', 'u', 'v'],
        physloc = 'v',
    ),
    ugbt = dict(
        names=['ugbt'],
        standard_names=['barotropic_sea_water_x_geostrophic_velocity',
            'eastward_geostrophic_current_velocity'],
        long_names = "Sea water barotropic geostrophic velocity along X",
        atlocs = ['t', 'u', 'v'],
        physloc = 'u',
        units = "m s-1",
    ),
    vgbt = dict(
        names=['vgbt'],
        standard_names=['barotropic_sea_water_y_geostrophic_velocity',
            'northward_geostrophic_current_velocity'],
        long_names = "Sea water barotropic geostrophic velocity along Y",
        atlocs = ['t', 'u', 'v'],
        physloc = 'v',
        units = "m s-1",
    ),
    speed = dict(
        names = ['speed'],
        standard_names=['sea_water_speed'],
        long_names = ["Sea water speed"],
        units = "m s-1",
        axes = dict(y=['lat'], x=['lon']),
        atlocs = ['t', 'u', 'v'],
        physloc = 't',
    ),
    cdir = dict(
        names = ['cdir'],
        standard_names=['direction_of_sea_water_velocity'],
        long_names = ["Direction of sea water velocity"],
        units = "degrees",
        axes = dict(y=['lat'], x=['lon']),
        atlocs = ['t', 'u', 'v'],
        physloc = 't',
    ),
    ke = dict(
        names = ['ke'],
        standard_names ='kinetic_energy',
        long_names = 'Kinetic energy',
        units = "m2 s-2",
    ),
    eke = dict(
        names = ['eke'],
        standard_names ='eddy_kinetic_energy',
        long_names = 'Eddy kinetic energy',
        units = "m2 s-2",
    ),
    tke = dict(
        names = ['tke'],
        standard_names ='turbulent_kinetic_energy',
        long_names = 'Turbulent kinetic energy',
        units = "m2 s-2",
    ),
    mke = dict(
        names = ['mke'],
        standard_names ='mean_kinetic_energy',
        long_names = 'Mean kinetic energy',
        units = "m2 s-2",
    ),
    kz = dict(
        names = ['kz', 'kzm'],
        standard_names = 'average_ocean_vertical_tracer_diffusivity',
        long_names = "Vertical diffusivity",
        units = "m2 s-1",
    ),


    # Bathymetry
    bathy = dict(
        names=['bathy', 'h0'],
        standard_names=['model_sea_floor_depth_below_sea_level', 'model_sea_floor_depth_below_geoid', "sea_floor_depth_below_geoid" ],
        long_names = 'Bathymetry',
        units = 'm',
        atlocs = 't',
    ),
    bathy_u = dict(
        names = ['hx'],
        standard_names = [],
        long_names =['bathymetry at u-location'],
    ),
    bathy_v = dict(
        names = ['hy'],
        standard_names = [],
        long_names =['bathymetry at v-location'],
    ),

    depth = dict(
        names = ['depth', 'dep', 'deptht', 'depthu', 'depthv'],
        standard_names = ['ocean_layer_depth'],
        long_names =  'Depth',
        units = 'm',
        atlocs = ['t', 'u', 'v', 'w'],
    ),

    # Cell sizes
    dx = dict(
        names = ['dx'], #, 'dx_u', 'dx_v', 'dx_f'],
        standard_names = 'cell_x_size',
        long_names = "Mesh size along x",
        units = 'm',
        atlocs = ['t', 'u', 'v', 'f'],
    ),
    dy = dict(
        names = ['dy'], #, 'dy_u', 'dy_v', 'dy_f'],
        standard_names = 'cell_y_size',
        long_names = "Mesh size along y",
        units = 'm',
        atlocs = ['t', 'u', 'v', 'f'],
    ),
    dz = dict(
        names = ['dz'],
        standard_names = 'ocean_layer_thickness',
        long_names = "Ocean layer thickness",
        units = "m",
        atlocs = ['t', 'u', 'v', 'w'],
    ),
    dlon = dict(
        names = ['dlon'], #, 'dlon_u', 'dlon_v', 'dlon_f'],
        standard_names = 'cell_x_size',
        long_names = "Mesh size along x",
        units = 'degrees',
        atlocs = ['t', 'u', 'v', 'f'],
    ),
    dlat = dict(
        names = ['dlat'], #, 'dlat_u', 'dlat_v', 'dlat_f'],
        standard_names = 'cell_y_step',
        long_names = "Mesh step along y",
        units = 'degrees',
        atlocs = ['t', 'u', 'v', 'f'],
    ),

    # Cell volume
    cvol = dict(
        names = ['cvol'],
        standard_names = 'cell_volume',
        long_names = "Volume of the cell",
        units = "m3",
        atlocs = ['t', 'u', 'v', 'w'],
    ),

    # Standard volume
    vol = dict(
        names = ['vol'],
        standard_names = 'seawater_volume',
        long_names = "Volume of the sea water",
        units = "m3",
        atlocs = ['t', 'u', 'v', 'w'],
    ),

    # Coriolis
    corio = dict(
        names = ['corio', 'f0'],
        standard_names = "coriolis_parameter",
        long_names = "Coriolis parameter",
        units = "s-1",
        atlocs = ['t', 'u', 'v', 'f'],
    ),
    beta = dict(
        names = [],
        standard_names = "meridional_derivative_of_coriolis_parameter",
        long_names = "Meridional derivative of coriolis parameter",
        units = "m-1 s-1",
    ),

    # Atmosphere
    lwhf = dict(
        names = [],
        standard_names = ["surface_net_downward_longwave_flux"],
        long_names = "Net longwave radiation (positive when directed downward)",
        units = "W.m-2",
    ),
    swhf = dict(
        names = [],
        standard_names = ["surface_net_downward_shortwave_flux"],
        long_names = "Net shortwave radiation (positive when directed downward)",
        units = "W.m-2",
    ),
    lathf = dict(
        names = [],
        standard_names = ["surface_downward_latent_heat_flux"],
        long_names = "latent heat flux (positive when directed downward)",
        units = "W.m-2",
    ),
    senhf = dict(
        names = [],
        standard_names = ["surface_downward_sensible_heat_flux"],
        long_names = "sensible heat flux (positive when directed downward)",
        units = "W.m-2",
    ),
    evap = dict(
        names = [],
        standard_names = ["lwe_thickness_of_water_evaporation_amount"],
        long_names = "evaporation (positive when directed downward)",
        units = "m",
    ),
    rain = dict(
        names = [],
        standard_names = ["lwe_thickness_of_precipitation_amount"],
        long_names = "precipitation (positive when directed downward)",
        units = "m",
    ),
    wspd = dict(
        names=[],
        standard_names = 'wind_speed',
        long_names = "Wind speed",
        units = "m s-1",
    ),

    wdir = dict(names=[],
        standard_names = ['wind_to_direction', 'wind_from_direction'],
        long_names = "Wind direction",
        units = "degrees",
    ),
    wfdir = dict(names=[],
        standard_names = ['wind_from_direction'],
        long_names = "Wind from direction",
        units = "degrees",
    ),
    wtdir = dict(names=[],
        standard_names = ['wind_to_direction'],
        long_names = "Wind to direction",
        units = "degrees",
    ),
    ua = dict(
        names = [],
        standard_names = ["eastward_wind", "x_wind"],
        long_names = "Zonal wind speed (westerly)",
        units = "m s-1",
    ),
    va = dict(
        names = [],
        standard_names = ["northward_wind", "y_wind"],
        long_names = "Meridional wind speed (northerly)",
        units = "m s-1",
    ),

    # Ocean Atmosphere interface
    u10m = dict(
        names = [],
        standard_names = ["x_wind_at_10m","x_wind","x_wind_at_u_location",
            "x_wind_at_10m_at_u_location","eastward_wind"],
        long_names = "10-m zonal wind speed (westerly)",
        units = "m s-1",
    ),
    v10m = dict(
        names = [],
        standard_names = ["y_wind_at_10m","y_wind","y_wind_at_v_location",
            "y_wind_at_10m_at_v_location","northward_wind"],
        long_names = "10-m meridional wind speed (northerly)",
        units = "m s-1",
    ),
    ux10m = dict(
        names = [],
        standard_names = ["x_wind_at_10m", "x_wind", "grid_eastward_wind",
                          "x_wind_at_u_location",
                          "x_wind_at_10m_at_u_location"],
        long_names = "10-m wind speed along X",
        units = "m s-1",
    ),
    vy10m = dict(
        names = [],
        standard_names = ["y_wind_at_10m", "y_wind", "grid_northward_wind",
                          "y_wind_at_v_location",
                          "y_wind_at_10m_at_v_location"],
        long_names = "10-m wind speed along Y",
        units = "m s-1",
    ),

    tauu = dict(
        names = [],
        standard_names = ["surface_downward_eastward_stress",
            "surface_eastward_stress"],
        long_names = "Surface eastward wind stress",
        units = "N m-2",
        axes = dict(x=['lon_u'], y=['lat_u']),
        physloc = 'u',
    ),
    tauv = dict(
        names = [],
        standard_names = ["surface_downward_northward_stress",
            "surface_northward_stress"],
        long_names = "Surface northward wind stress",
        units = "N m-2",
        axes = dict(x=['lon_v'], y=['lat_v']),
        physloc = 'v',
    ),
    taux = dict(
        names = ['ustress'],
        standard_names = ["surface_downward_x_stress", "surface_x_stress",
                          "surface_downward_x_stress_at_u_location"],
        long_names = "Surface wind stress along X",
        units = "N m-2",
        axes = dict(x=['lon_u'], y=['lat_u']),
        physloc = 'u',
    ),
    tauy = dict(
        names = ['vstress'],
        standard_names = ["surface_downward_y_stress", "surface_y_stress",
                          "surface_downward_y_stress_at_v_location"],
        long_names = "Surface wind stress along Y",
        units = "N m-2",
        axes = dict(x=['lon_v'], y=['lat_v']),
        physloc = 'v',
    ),

    # Surfaces waves
    hs = dict(
        names = ['hs'],
        standard_names = ["significant_height_of_wind_and_swell_waves"],
        long_names = "Significant height of wind and swell waves",
        units = "m",
    ),
    fp = dict(
        names = ['fp'],
        standard_names = [],
        long_names = "Frequency of wind and swell waves at spectral peak",
        units = "s-1",
    ),
    th1p = dict(
        names = ['th1p'],
        standard_names = ["sea_surface_wave_from_direction"],
        long_names = "Mean direction of wind and swell waves at spectral peak",
        units = "degree",
    ),


)
var_specs = VAR_SPECS # compat

#: Specifications for axes
AXIS_SPECS = OrderedDict(

    # Axes
    time = dict(
        names = ['time'],
        standard_names = ['time'],
        long_names = 'Time',
        axis = 'T',
    ),
    lon = dict(
        names = ['lon', 'longitude', 'nav_lon'],
        standard_names = ['longitude'],
        long_names = 'Longitude',
        units = ['degrees_east', 'degree_east', 'degree_e', 'degrees_e', 'degreee', 'degreese'],
        axis = 'X',
        iaxis = 'ni',
        jaxis = 'nj',
        atlocs=['t', 'u', 'v', 'f'],
    ),
    lat = dict(
        names = ['lat', 'latitude', 'nav_lat'],
        standard_names = ['latitude'],
        long_names = 'Latitude',
        units = ['degrees_north', 'degree_north', 'degree_n', 'degrees_n', 'degreen', 'degreesn'],
        axis = 'Y',
        iaxis = 'ni',
        jaxis = 'nj',
        atlocs=['t', 'u', 'v', 'f'],
    ),
    depth = dict(
        inherit = 'depth',
        axis = 'Z',
    ),
    depth_t = dict(
        inherit = 'depth_t',
        axis = 'Z',
    ),
    depth_u = dict(
        inherit = 'depth_u',
        axis = 'Z',
    ),
    depth_v = dict(
        inherit = 'depth_v',
        axis = 'Z',
    ),
    depth_w = dict(
        inherit = 'depth_w',
        axis = 'Z',
    ),
    level = dict(
        names = ['level'],
        standard_names = ["model_level_number", "ocean_sigma_coordinate", "ocean_s_coordinate", "ocean_sigma_coordinate_at_w_location", "ocean_s_coordinate_at_w_location"],
        long_names = ['Model level number',  'Sigma level', 'Sigma level at W location'],
        axis='Z',
    ),
    level_w = dict(
        names = ['level_w'],
        standard_names = ["ocean_sigma_coordinate_at_w_location", "ocean_s_coordinate_at_w_location"],
        long_names = ['Sigma level', 'Sigma level at W location'],
        axis='Z',
    ),

    # Subaxes
    ni = dict(
        names = ['ni'],
        standard_names =  ['x_grid_index'],
        long_names = ["x-dimension of the grid"],
        axis = 'X',
    ),
    nj = dict(names = ['nj'],
        standard_names =  ['y_grid_index'],
        long_names = ["y-dimension of the grid"],
        axis = 'Y',
    ),
    ni_u = dict(
        names = ['ni_u'],
        standard_names =  ['x_grid_index_at_u_location'],
        long_names = ["x-dimension of the grid at U location"],
        axis = 'X',
    ),
    nj_u = dict(names = ['nj_u'],
        standard_names =  ['y_grid_index_at_u_location'],
        long_names = ["y-dimension of the grid at U location"],
        axis = 'Y',
    ),
    ni_v = dict(names = ['ni_v'],
        standard_names =  ['x_grid_index_at_v_location'],
        long_names = ["x-dimension of the grid at V location"],
        axis = 'X',
    ),
    nj_v = dict(names = ['nj_v'],
        standard_names =  ['y_grid_index_at_v_location'],
        long_names = ["y-dimension of the grid at V location"],
        axis = 'Y',
    ),



)
axis_specs = AXIS_SPECS # compat

#: Specifications for grid formating
GRID_SPECS = {
    't': dict(lon='lon', lat='lat', level='depth_t'),
    'u': dict(lon='lon_u', lat='lat_u', level='depth_t'),
    'v': dict(lon='lon_v', lat='lat_v', level='depth_t'),
    'f': dict(lon='lon_f', lat='lat_f', level='depth_t'),
    'w': dict(lon='lon', lat='lat', level='depth_w'),
}
#TODO: 't' must use lon_t and lat_t
GRID_SPECS['r'] = GRID_SPECS['t']
grid_specs = GRID_SPECS # compat

_reloc =  OrderedDict(
    [('name',  re.compile('(_([a-z]))?$', re.I).search),
    ('standard_name', re.compile('(_at_([a-z])_location)?$', re.I).search),
    ('long_name',  re.compile('( at ([a-z]) location)?$', re.I).search)],
)

#: Default location of gridded variables
DEFAULT_LOCATION = 't'
default_location = DEFAULT_LOCATION # compat

def get_loc(var, stype=None, mode='loc', default=None):
    """Get the location testing name, standard_name or long_name and other attributes

    :Params:

        - **var**: Generic name of the variable or axis OR an object CDAT object
          (variable or axis). If a CDAT object, it first searches for the
          :attr:`_vacumm_cf_location` attributes.
        - **stype**, optional: One of 'name', 'standard_name' or 'long_name'.
          If ``None``, they are all tested.
          If ``var`` is a string, location is searched interpreting it as of type ``'stype``.
          Else if ``var`` is an object ``stype`` specifies what attribute to test.
        - **mode**, optional: What you get if something found.

            - ``"loc"``: Get the location lowercase letter (like "u").
            - ``"ext"``: Same as ``"loc"`` but searches for the :attr:`position` first
              if case of a CDAT variable,
              then the :attr:`_vacumm_cf_physloc` attributes or the "physloc" specification
              if ``var`` has an entry in :attr:`VAR_SPECS`, else defaults to
              :attr:`DEFAULT_LOCATION`.
            - Else, the matched string.

    :Example:

        >>> loc = get_loc('toto_u')
        >>> suffix = get_loc('toto_at_u_location', stype='standard_name', mode='full')
        >>> loc = get_loc(myvar, stype=['name','long_name'])

    :Return:

        - ``None`` if no location spec found,
        - if ``mode=="loc"``, one of possible physical
          :attr:`~vacumm.data.misc.arakawa.ARAKAWA_LOCATIONS`,
        - else, the complete matching string, like "_at_u_location".
    """
    # What to test
    if stype is None: stype = _reloc.keys() # all by default
    stypes = stype if isinstance(stype, list) else [stype]

    # CDAT object (not a string)
    if isinstance(var, list):
        if len(var)==0:
            return default
        var = var[0]
    if not isinstance(var, basestring):

        # First: _vacumm_cf_location
        if getattr(var, '_vacumm_cf_location', None):
            return var._vacumm_cf_location.lower()

        # Location from names
        for stype in stypes:
            if stype=='name':
                loc = get_loc(var.id, stype, mode=mode)
            elif hasattr(var, stype):
                loc = get_loc(getattr(var, stype), stype, mode)
            else:
                continue
            break
        if loc is not None:
            return loc

        # Ext mode
        if mode=='ext':

            # Position attribute
            if hasattr(var, 'position'): return var.position.lower()

            # Physloc
            gname = var.id if not hasattr(var, '_vacumm_cf_name') else var._vacumm_cf_name
            gname = no_loc_single(gname, 'name')
            if gname in VAR_SPECS:

                # From _vacumm_cf_physloc attribute
                if hasattr(var, '_vacumm_cf_physloc'):
                    return var._vacumm_cf_physloc

                # From specs
                if VAR_SPECS[gname].get('physloc'):
                    return VAR_SPECS[gname]['physloc']

        return default

    # String
    outloc = mode=='loc' or mode=='ext'
    for stype in stypes:

        # From string
        group = 2 if outloc else 1
        if stype not in _reloc: stype = 'name'
        loc = _reloc[stype](var).group(group)
        if loc and outloc: loc = loc.lower()
        if loc: return loc

        # Extended mode
        if mode=='ext':

            # From physloc
            if stype=='name' and var in VAR_SPECS and VAR_SPECS[var].get('physloc'):
                return VAR_SPECS[var].get('physloc')

#            # Default location
#            return DEFAULT_LOCATION
    return default


def no_loc_single(name, stype):
    """Remove location specification

    :Params:

        - **name**: Generic ame of the variable or axis.
        - **stype**: One of 'name', 'standard_name' or 'long_name'.
    """
    loc = get_loc(name, stype, mode='full')
    if loc is not None:
        return name[:-len(loc)]
    return name

def _loc_(loc=None):
    """Get location as a single char or empty lowercase string"""
    if loc is None: loc = ''
    if loc.startswith('_'): loc = loc[1:]
    loc = loc[:1]
    if loc not in string.ascii_letters:
        raise TypeError('Wrong location: '+loc)
    return loc.lower()

def change_loc_single(name, stype, loc, squeeze=None):
    """Change location specification

    :Params:

        - **name**: String to change
        - **stype**: String type: name, standard_name or long_name.
        - **loc**: Location as a letter or None.
        - **squeeze**, optional: If specified and equal to ``loc``, location is removed
          instead of being changed.

    :Example:

        >>> change_loc_single('usurf_u', 'name', 't')
        'usurf_t'
        >>> change_loc_single('sst_at_u_location', 'standard_name', 'v')
        'sst_at_v_location'
        >>> change_loc_single('usurf_t', 'name', 'u', squeeze='u')
        'usurf'
    """
    basename = no_loc_single(name, stype)
    loc = _loc_(loc)
    if loc and (not squeeze or squeeze!=loc):
        if stype=='standard_name':
            return basename+'_at_%s_location'%loc
        if stype=='long_name':
            return basename+' at %s location'%loc.upper()
        return basename+'_'+loc
    return basename

def change_loc_specs(loc, names=None, standard_names=None, long_names=None, axes=None,
    iaxis=None, jaxis=None, squeeze=None, **kwargs):
    """Change location specification in names, standard names, long names and axes names

    :Example:

        >>> specs = change_loc_specs('v', names = ['usurf_u', 'u_u'])

    :Return: A dictionary
    """
    specs = kwargs.copy()
    if 'atlocs' in specs: del specs['atlocs']

    # Attributes
    for stype in 'name', 'standard_name', 'long_name', 'iaxis', 'jaxis':
        sname = stype if stype in ('iaxis', 'jaxis') else (stype+'s')
        values = eval(sname)
        if values is None: continue
        if isinstance(values, list):
            values = list(values)
            tmp = [change_loc_single(value, stype, loc, squeeze=squeeze) for value in values]
            values = []
            for value in tmp: # unique
                if value not in values:
                    values.append(value)

        else:
            values = change_loc_single(values, stype, loc, squeeze=squeeze)
        specs[sname] = values

    # Axes
    if axes is not None:
        axes = axes.copy()
        for l in axes.keys(): # loop on x, y
            if l=='t': continue # skip time
            laxes = axes[l]
            single = not isinstance(laxes, list)
            if single: laxes = [laxes] # work on lists
            laxes = [change_loc_single(axis, 'name', loc, squeeze=DEFAULT_LOCATION) for axis in laxes]+laxes # duplicate
            lnewaxes = []
            for axis in laxes: # unique
                if axis not in lnewaxes:
                    lnewaxes.append(axis)
            if single and len(lnewaxes)==1:
                lnewaxes = lnewaxes[0]
            axes[l] = lnewaxes
        specs['axes'] = axes

    return specs

def change_loc(var, toloc, axes=True, squeeze=True):
    """Change location specifications of a variable and its axes

    It affects the id, the standard_name and the long_name when they are defined.

    :Params:

        -

    :Example:

        >>> change_loc(sst, 'u').id
        'sst_u'
        >>> change_loc(sst, 't', squeeze=True).id
        'sst'

    """

    # Squeeze physloc
    if not isinstance(squeeze, basestring):
        if squeeze:
            squeeze = get_physloc(var)
        else:
            squeeze = None

    # Change attributes
    # - get
    specs = change_loc_specs(toloc, squeeze=squeeze,
        names=var.id,
        standard_names = getattr(var, 'standard_name', None),
        long_names = getattr(var, 'long_name', None),
        )
    # - set
    var.id = specs['names']
    for att in 'standard_name', 'long_name':
        if att+'s' in specs:
            setattr(var, att, specs[att+'s'])
    # - cf name
    if hasattr(var, '_vacumm_cf_name'):
        var._vacumm_cf_name = change_loc_single(var._vacumm_cf_name, 'name', toloc,
            squeeze=squeeze)

    # Change axes and grid attributes
    if axes and cdms2.isVariable(var) or cdms2.isGrid(var):

        # Axes
        # - usual
        for meth in 'Longitude', 'Latitude', 'Time', 'Level':
            if hasattr(var, 'get'+meth):
                axismet = getattr(var, 'get'+meth) # var.getLongitude
                if axismet() is not None:
                    change_loc(axismet(), toloc, squeeze=squeeze)
        # - 2d axes
        if cdms2.isVariable(var) and isaxis(var) and var.ndim>2:
            for i in -1, -2:
                change_loc(var.getAxis(i), toloc, squeeze=squeeze)

        # Grid
        if cdms2.isVariable(var) and var.getGrid() is not None:
            change_loc(var.getGrid(), toloc, squeeze=squeeze)

    # Reference attribute
    set_loc(var, toloc)

    return var

def set_loc(var, loc, addpos=False):
    """Define (or remove) the location of a variable with :attr:`_vacumm_cf_location` attribute

    It also makes sure that the :attr:`position` attribute is the same if present.

    :Params:

        - **var**: CDAT variable.
        - **loc**: Location. If empty, remove :attr:`_vacumm_cf_location` and
          :attr:`position` attributes.
        - **addpos**, optional: Set also the :attr:`position` attribute.

    """
    # Remove
    if not loc:
        for att in '_vacumm_cf_location', 'position':
            if hasattr(var, att):
                delattr(var, att)

    # Set
    loc = loc.lower()
    var._vacumm_cf_location = loc
    if addpos:
        var.position = loc
    elif hasattr(var, 'position') and var.position.lower()!=loc:
        var.position = loc
    return var


def get_physloc(var):
    """Get the physical location of a variable

    .. note:: The physical location may be different from the current location.

    It defaults to :attr:`DEFAULT_LOCATION`
    """
    # Direct
    for att in ['_vacumm_cf_physloc']:
        if hasattr(var, att):
            return getattr(var, att)

    # Entry name in VAR_SPECS (string)
    if isinstance(var, basestring):
        if var in VAR_SPECS:
            var = VAR_SPECS[var]
        else:
            return DEFAULT_LOCATION

    # Dict of specs
    if isinstance(var, dict):
        return var.get('physloc', DEFAULT_LOCATION)

    # Using VAR_SPECS
    if hasattr(var, '_vacumm_cf_name'):
        name = var._vacumm_cf_name
    elif var.id in VAR_SPECS:
        name = var.id
    else:
        return DEFAULT_LOCATION
    return VAR_SPECS[name].get('physloc', DEFAULT_LOCATION)

def squeeze_loc(var):
    """Squeeze out location specification of a variable if the location is
    physloc or :attr:`DEFAULT_LOCATION`

    :Params:

        - **var**: CDAT variable.
        - **physloc**, optional: ``physloc`` specs given by the entry in :attr:`VAR_SPECS`
          or :attr:`DEFAULT_LOCATION`.

    :Return: The same variable
    """
    # Physloc
    physloc = get_physloc(var)

    # Remove loc only for physloc
    change_loc(var, None, axes=True, squeeze=physloc)

def squeeze_loc_single(name, stype, physloc=None):
    """Squeeze location specs if it matches physical location"""
    if physloc is None and stype=='name' and name in VAR_SPECS:
        physloc = VAR_SPECS[name].get('physloc')
    if physloc is None:
        physloc = DEFAULT_LOCATION
    loc = get_loc(name, stype, mode='loc')
    if loc==physloc:
        name = change_loc_single(name, stype, None)
    return name

def dupl_loc_specs(all_specs, fromname, toloc):
    """Duplicate the specification for a variable or an axis to another or several locations

    The following rules apply:

        - If the original specifications are from a name without a specific location
          (generic), new specs (located) are appended (merged) with original specs.
        - If the specifications at target location already exist,
          it merges new specs with old specs.
        - Generic specification are systematically created or updated by merging of
          specialized ones.

    :Example:

        >>> dupl_loc_specs(VAR_SPECS, 'corio', 'u') # Create the 'corio_u' entry in VAR_SPECS and update 'corio'
        >>> dupl_loc_specs(VAR_SPECS, 'corio_t', 'u') # Create 'corio_u' and 'corio_t'

    """
    if not fromname in all_specs:
        raise KeyError('No such entry in specifications: '+fromname)
    single = not isinstance(toloc, (list, tuple))
    if single:
        toloc = [toloc]
    tonames = []
    tomerge = []
    fromnoloc = no_loc_single(fromname, 'name')==fromname
    for loc in toloc:

        # New name (id)
        toname = change_loc_single(fromname, 'name', loc)
        tonames.append(toname)

        # New specs
        tospecs = change_loc_specs(loc, **all_specs[fromname])

        # Add a version of standard_name and long_name for T loc without location spec
        if loc=='t': # 'toto_at_t_location' -> ['toto_at_t_location','toto']
            for spname in 'standard_names', 'long_names':
                if spname not in tospecs: continue
                addspecs = change_loc_specs(None, **{spname:tospecs[spname]})
                for dd in tospecs, addspecs: # to lists
                    if not isinstance(dd[spname], list):
                        dd[spname] = [dd[spname]]
                tospecs[spname].extend(addspecs[spname])

        # For merging specs back with old specs if old has no location
#        if fromnoloc: # generic
        tomerge.append(tospecs)

        # Merge with old spec if existing
        if toname in all_specs:
            tospecs = dict_merge(tospecs, all_specs[toname], mergelists=True)

#        # Default location -> add its generic version ('_t' -> + '')
#        if loc==DEFAULT_LOCATION:
#            genspecs = change_loc_specs(None, **tospecs)
#            tospecs = dict_merge(tospecs, genspecs, mergelists=True)

        # Remove atlocs attribute
        if 'atlocs' in tospecs:
            tospecs.pop('atlocs')

        # Store it
        all_specs[toname] = tospecs

    # Make a generic entry that merges all specialized ones
    if fromnoloc:
        genspecs = all_specs[fromname] # existing generic specs
    else:
        genspecs = change_loc_specs(None, **all_specs[fromname]) # remove location info
    tomerge.insert(0, genspecs)
    kw = dict(mergelists=True)
    genname = genspecs['names'][0]
    all_specs[genname] = dict_merge(*tomerge, **kw)

    if single: return tonames[0]
    return tonames

#def specs_def_loc(all_specs, names, suffixes=ARAKAWA_SUFFIXES):
#    """Make specs of a variable at a special location the specs for default location
#
#    :Example:
#
#        >>> specs_def_loc(all_specs, ['usurf_u', 'ssh_f']) # specs of 'usurf_u' is copied to 'usurf'
#    """
#    if not isinstance(names, (list, tuple)): names = [names]
#    if isinstance(suffixes, basestring): suffixes = [suffixes]
#    for name in names:
#        m = re.match('.+(%s)$'%'|'.join(suffixes), name)
#        if m is None: continue
#        dupl_loc_specs(all_specs, name, '', suffixes=ARAKAWA_SUFFIXES)


def cp_suffix(idref, id, suffixes=ARAKAWA_SUFFIXES):
    """Copy a suffix if found in an id (name) to another id"""
    if isinstance(suffixes, basestring): suffixes = [suffixes]
    m = re.match('.+(%s)$'%'|'.join(suffixes), idref)
    if m is None: return id
    return id+m.group(1)

# Format specifications
# - Makes sure to have lists, except for 'axis' and 'inherit'
# - Check geo axes
# - Check inheritance
# - Makes sure that axes specs have no 'axes' key
# - Makes sure that specs key is the first entry of 'names'
# - add standard_name to list of names
# - Check duplication to other locations ('toto' -> 'toto_u')
for all_specs in VAR_SPECS, AXIS_SPECS:

    from_atlocs = []
    for name, specs in all_specs.items():

        # Entry already generated with the atlocs key
        if name in from_atlocs: continue

        # Always lists (except for dict and some keys)
        for key, value in specs.items():
            if isinstance(value, dict): continue
            if key in ['axis', 'inherit', 'physloc']: continue
            if isinstance(value, tuple):
                specs[key] = list(specs[key])
            elif not isinstance(value, list):
                specs[key] = [value]

        # Physloc must be the first atlocs
        if 'physloc' in specs and 'atlocs' in specs:
            p = specs['physloc']
            if p in specs['atlocs']:
                specs['atlocs'].remove(p)
            specs['atlocs'].insert(0, p)

        # Geo axes (variables only)
        if all_specs is VAR_SPECS:
            if 'axes' not in specs or not specs['axes']:
                specs['axes'] = {}
            specs['axes'].setdefault('t', 'time')
            suffixes = [('_'+s) for s in 'rftuv']
            specs['axes'].setdefault('y', cp_suffix(name, 'lat', suffixes=suffixes))
            specs['axes'].setdefault('x', cp_suffix(name, 'lon', suffixes=suffixes))
            if 'atlocs' in specs: # need list for later duplication of location
                for l in 'yx':
                    if not isinstance(specs['axes'][l], list):
                        specs['axes'][l] = [specs['axes'][l]]
            for l, n in ('t', 'time'), ('y', 'lat'), ('x', 'lon'):
                if isinstance(specs['axes'][l], list) and n not in specs['axes'][l]:
                    specs['axes'][l].append(n)

        # Inherits from other specs (merge specs with dict_merge)
        if 'inherit' in specs:
            objname = specs['inherit']
            obj_specs = None
            if name==objname: # same name
                if all_specs is VAR_SPECS and objname in AXIS_SPECS:
                    obj_specs = AXIS_SPECS
                elif all_specs is AXIS_SPECS and objname in VAR_SPECS:
                    obj_specs = VAR_SPECS
            elif objname in VAR_SPECS: # check in VAR_SPECS first
                obj_specs = VAR_SPECS
            elif objname in AXIS_SPECS: # then in AXIS_SPECS
                obj_specs = AXIS_SPECS
            if obj_specs is not None:
                all_specs[name] = specs = dict_merge(specs, obj_specs[objname])

        # No geo axes for axes!
        if all_specs is AXIS_SPECS and 'axes' in specs:
            del specs['axes']

        # Key = first name (id)
        if "names" not in specs:
            specs['names'] = [name]
        else:
            if name in specs['names']:
                specs['names'].remove(name)
            specs['names'].insert(0, name)

        # Standard_names in names
        if ('standard_names' in specs and specs['standard_names']):
            for standard_name in specs['standard_names']:
                if standard_name not in specs['names']:
                    specs['names'].append(standard_name)

        # Duplicate at other locations
        if 'atlocs' in specs:
            tonames =  dupl_loc_specs(all_specs, name, specs['atlocs'])
            from_atlocs.extend(tonames)

del specs

#: List of generic variable names
GENERIC_VAR_NAMES = VAR_SPECS.keys()
generic_var_names = GENERIC_VAR_NAMES # compat

#: List of generic axis names
GENERIC_AXIS_NAMES = AXIS_SPECS.keys()
generic_axis_names = GENERIC_AXIS_NAMES # compat

#: List of all generic names (axes and variables)
GENERIC_NAMES= GENERIC_VAR_NAMES + GENERIC_AXIS_NAMES
generic_names = GENERIC_NAMES # compat

def cf2search(name, mode=None, raiseerr=True, **kwargs):
    """Extract specs from :attr:`AXIS_SPECS` or :attr:`VAR_SPECS` to form a search dictionary

    :Params:

        - **name**: Generic name of an axis or a variable.
        - **mode**, optional: Search mode [default: None->``"ns"``].
          A string containg one or more of the following letters:

            - ``"n"``: Search using names (ids).
            - ``"s"``: Search using standard_name attribute.
            - ``"a"``: Search using axis attribute.
            - ``"l"``: Search using long_name attribute.
            - ``"u"``: Search using units attribute.

          The order is important.

    :Return: An :class:`colections.OrderedDict`

    :Example:

        >>> cf2search('sst', mode='nsu')
        {'names':['sst'], 'standard_names':['sea_surface_temperature'], 'units':['degrees_celsius']}
    """
    # Get specs
    if name in VAR_SPECS:
        specs = VAR_SPECS[name]
    elif name in AXIS_SPECS:
        specs = AXIS_SPECS[name]
    else:
        if raiseerr:
            raise VACUMMError("Wrong generic name. It should be one of: "+' '.join(GENERIC_AXIS_NAMES+GENERIC_VAR_NAMES))
        else:
            return

    # Form search dict
    if not isinstance(mode, basestring): mode = 'nsa'
    keys = []
    for m in mode:
        for key in ['names', 'standard_names', 'axis', 'long_names', 'units']:
            if key.startswith(m):
                keys.append(key)
                break
    return OrderedDict([(k, specs[k]) for k in keys if k in specs])


_attnames_plurals = ['standard_name', 'long_name']
_attnames_exclude = ['names', 'atlocs', 'inherit', 'axes', 'physloc', 'iaxis', 'jaxis', 'select']
_attnames_firsts = ['standard_name', 'long_name', 'units', 'axis', 'valid_min', 'valid_max']
def cf2atts(name, select=None, exclude=None, ordered=True, **extra):
    """Extract specs from :attr:`AXIS_SPECS` or :attr:`VAR_SPECS` to form
    a dictionary of attributes (units and long_name)"""
    # Get specs
    if isinstance(name, dict):
        specs = name.copy()
    elif name in VAR_SPECS:
        specs = VAR_SPECS[name]
    elif name in AXIS_SPECS:
        specs = AXIS_SPECS[name]
    else:
        raise VACUMMError("Wrong generic name: %s. It should be one of: "%name+' '.join(GENERIC_AXIS_NAMES+GENERIC_VAR_NAMES))

    # Which attributes
    atts = OrderedDict() if ordered else {}
    if exclude is None: exclude = []
    elif isinstance(exclude, basestring): exclude = [exclude]
    exclude.extend(_attnames_exclude)
    for attname in _attnames_firsts+specs.keys():

        # Distinction between specs key and attribute name
        if attname in _attnames_plurals:
            key = attname+'s'
        elif attname[:-1] in _attnames_plurals:
            key = attname
            attname = attname[:-1]
        else:
            key = attname

        # Skip some attributes
        if key not in specs or attname in exclude or attname in atts or \
            (select is not None and attname not in select):
            continue

        # No lists or tuples
        value = specs[key]
        if isinstance(value, (list, tuple)):
            if len(value)==0:
                continue
            value = value[0]

        # Store it
        atts[attname] = value

    # Extra
    for att, val in extra.items():
        atts[att] = val

    return atts



# Format a variable
def format_var(var, name, force=True, format_axes=True, order=None, nodef=True,
        mode='warn', **kwargs):
    """Format a MV2 variable according to its generic name


    :Params:

        - **var**: A :mod:`numpy` or :mod:`MV2` variabe.
        - **name**: Generic name of variable. It should be one of
          those listed by :attr:`GENERIC_VAR_NAMES`.
        - **force**, optional: Overwrite attributes in all cases.
        - **format_axes**, optional: Also format axes.
        - **nodef**, optional: Remove location specification when it refers to the
          default location (:attr:`DEFAULT_LOCATION`).
        - **mode**: "silent", "warn" or "raise".
        - Other parameters are passed as attributes, except those:

            - present in specifications to allow overriding defaults,
            - starting with 'x', 'y', or 't' which are passed
              to :func:`format_axis`.

    :Examples:

        >>> var = format_var(myarray, 'sst', valid_min=-2, valid_max=100)

    """
    # Filter keywords for axis formating
    axismeths = {'t':'getTime', 'y':'getLatitude', 'x':'getLongitude'}
    kwaxes = {}
    for k in axismeths.keys():
        kwaxes[k] = kwfilter(kwargs, k+'_')

    # Always a MV2 array
    if not cdms2.isVariable(var):
        var = MV2.asarray(var)

    # Check specs
    if name not in GENERIC_NAMES:
        if var.id in GENERIC_NAMES:
            name = var.id
        elif mode=='warn':
            warn("Generic var name not found '%s'."%name)
            return var
        elif mode=='silent':
            return var
        else:
            raise KeyError("Generic var name not found '%s'. Please choose one of: %s"%(
                name, ', '.join(GENERIC_NAMES)))
    isaxis = name in GENERIC_AXIS_NAMES
    if isaxis:
        specs = AXIS_SPECS[name].copy()
        if 'axis' in specs:
            del specs['axis']
    else:
        specs = VAR_SPECS[name].copy()
    # - merge kwargs and specs
    for key, val in kwargs.items():
        if val is None or key not in specs: continue
        # Check type
        if not isinstance(val, list) and isinstance(specs[key], list):
            val = [val]
        # Set
        specs[key] = val
        del kwargs[key]
    # - remove default location
    if nodef:
        refloc = specs.get('physloc', DEFAULT_LOCATION)
        for stype in 'name', 'long_name', 'standard_name':
            sname = stype+'s'
            if sname not in specs: continue
            if get_loc(specs[sname], stype)==refloc:
                specs[sname] = [no_loc_single(specs[sname][0], stype)]
        name = specs['names'][0]
    # - id
    if (force or var.id.startswith('variable_') or
            (isaxis and var.id.startswith('axis_'))): # FIXME: use regexp
        var.id = name
    # - attributes
    for att, val in cf2atts(specs, **kwargs).items():
        if force or not getattr(var, att, ''):
            setattr(var, att, val)
    # - physical location
    loc = get_loc(var, mode='ext')
    if not loc and 'physloc' in specs:
        loc = specs['physloc']
    if loc:
        if 'physloc' in specs and loc in specs['physloc']:
            var._vacumm_cf_physloc = loc.lower()
        set_loc(var, loc)
    # - store cf name
    var._vacumm_cf_name = name

    # Axes
    if format_axes:

        # Order
        order = var.getOrder() if not isinstance(order, basestring) else order
        if order is not None:
            if not re.match('^[xyzt-]+$', order):
                raise VACUMMError("Wrong cdms order type: "+order)
            if len(order)!=var.ndim:
                raise VACUMMError("Cdms order should be of length %s instead of %s"%(var.ndim, len(order)))

        # First check
        if 'axes' in specs:
            axspecs = specs['axes']
            formatted = []
            for key, meth in axismeths.items():
                axis = getattr(var, meth)()
                if order is not None: order.replace(key, '-')
                if axis is not None:
                    format_axis(axis, axspecs[key], **kwaxes[key])
                    formatted.append(key)

        # Check remaining simple axes (DOES NOT WORK FOR 2D AXES)
        if order is not None and order!='-'*len(order):
            for key in axismeths.keys():
                if key in order and key not in formatted:
                    axis = var.getAxis(order.index(key))
                    format_axis(axis, axspecs[key], **kwaxes[key])


    return var

# Format an axis
def format_axis(axis, name, force=True, recreate=False, format_subaxes=True,
    nodef=True, **kwargs):
    """Format a MV2 axis according to its generic name


    :Params:

        - **var**: A :mod:`numpy` or :mod:`MV2` variabe.
        - **name**: Single or list of generic axis names. It should be one of
          those listed by :attr:`GENERIC_AXIS_NAMES`.
        - **force**, optional: Overwrite attributes in all cases.
        - **nodef**, optional: Remove location specification when it refers to the
          default location (:attr:`DEFAULT_LOCATION`).
        - **axes2d_<param>**, optional: <param> is passed to
          :func:`~vacumm.misc.grid.misc.create_axes2d` for 2D axes.
        - **recreate**, optional: Force recreation of axes using either

        - Other parameters are passed as attributes.

    :Examples:

        >>> axis1d = format_axis(array1d, 'lon')
        >>> axis2d = format_axis(array2d, 'lon_u', axes2d_iid = 'xi',  axes2d_jid = 'xi', recreate=True)


    """

    # Guess best generic name from a list
    if name is None: return axis
    if isinstance(name, (list, tuple)):
        for nm in name:
            if match_obj(axis, nm):
                name = nm
                break

    # Check specs
    if name not in GENERIC_AXIS_NAMES:
        raise KeyError("Generic axis name not found '%s'. Please choose one of: %s"%(
            name, ', '.join(GENERIC_AXIS_NAMES)))
    specs = AXIS_SPECS[name]
    # - merge kwargs and specs
    for key, val in kwargs.items():
        if val is None or key not in specs: continue
        # Check type
        if not isinstance(val, list) and isinstance(specs[key], list):
            val = [val]
        # Set
        specs[key] = val
        del kwargs[key]

    # Always a MV2 axis (1D or 2D)
    axis._oldid = axis.id
    kwaxed2d = kwfilter(kwargs, 'axes2d_')
    if not isaxis(axis) or recreate:
        if len(axis.shape)==1:
            axis = cdms2.createAxis(axis)
        else:
            xy = specs['axis'].lower()
            kwaxed2d[xy] = axis
            axis = create_axes2d(**kwaxed2d)
            return axis
    axis2d = len(axis.shape)==2

    # Apply specs
    # - merge kwargs and specs
    for key, val in kwargs.items():
        if val is None or key not in specs: continue
        # Check type
        if not isinstance(val, list) and isinstance(specs[key], list):
            val = [val]
        # Set
        specs[key] = val
        del kwargs[key]
    # - remove default location
    if nodef:
        for stype in 'name', 'long_name', 'standard_name':
            sname = stype+'s'
            if sname not in specs: continue
            if get_loc(specs[sname], stype)==DEFAULT_LOCATION:
                specs[sname] = [no_loc_single(specs[sname][0], stype)]
    # - id
    if force or axis.id.startswith('variable_') or axis.id.startswith('axis_'):
        axis.id = specs['names'][0]
    # - attributes
    for att, val in cf2atts(specs, exclude=['axis'] if axis2d else None, **kwargs).items():
        if force or not getattr(axis, att, ''):
            setattr(axis, att, val)
    # - store cf name
    axis._vacumm_cf_name = axis.id

    # Sub-axes for 2D axes
    if axis2d and format_subaxes:
        format_axis(axis.getAxis(-1), specs['iaxis'][0])
        format_axis(axis.getAxis(-2), specs['jaxis'][0])

    return axis

def format_grid(grid, pt, **kwargs):
    """Format a grid and its axes"""
    if cdms2.isVariable(grid): grid = grid.getGrid()
    if grid is None: return
    gs = GRID_SPECS[pt]
    lon = grid.getLongitude()
    lat = grid.getLatitude()
    format_axis(lon, gs['lon'])
    format_axis(lat, gs['lat'])


#: Hidden attributes of variable useful of this module
HIDDEN_CF_ATTS = ['_vacumm_cf_name', '_vacumm_cf_physloc', '_vacumm_cf_location']
hidden_cf_atts = HIDDEN_CF_ATTS # compat

def match_obj(obj, name, searchmode=None, **kwargs):
    """Check if a variable or an axis matches generic specifications

    :Params:

        - **obj**: a numpy or cdms2 axis or variable.
        - **name**: A generic names.
        - **searchmode**, optional: Passed to :func:`~vacumm.misc.io.ncmatch_obj`.
    """
    search = cf2search(name, mode=searchmode, raiseerr=False)
    if search is None: return False
    search.update(kwargs)
    return ncmatch_obj(obj, searchmode=searchmode, **search)

match_var = match_obj
