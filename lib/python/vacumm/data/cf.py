"""Conventions about data formats and names"""

import cdms2, MV2, re
from collections import OrderedDict

from vacumm import VACUMMError
from vacumm.misc import kwfilter, dict_merge
from vacumm.misc.axes import create as create_axis, isaxis
from vacumm.misc.grid import create_axes2d
from vacumm.misc.io import ncmatch_obj

__all__ = ['var_specs', 'axis_specs', 
    'generic_axis_names', 'generic_var_names', 'generic_names', 
    'format_var', 'format_axis', 'format_grid', 'match_var', 
    'cf2atts', 'cf2search', 'cp_suffix', 
]

#: Specifications for variables
var_specs = dict(

    # Thermodynamics
    temp = dict(
        names = ['temp', 'temperature'],
        standard_names = ['sea_water_temperature', 'sea_water_potential_temperature'],
        long_names = 'Temperature',
        units =  'degrees_celsius', 
#        axes = {'t':'time', 'x':'lon', 'y':'lat'}, 
    ),
    sal = dict(
        names=['sal', 'psal', 'salinity'],
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
    pdens = dict(
        names = ['pdens'], 
        standard_names = 'sea_water_potential_density', 
        long_names = 'Sea water potential density', 
        units = 'kg m-3',
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
        standard_names=['sea_water_x_velocity_at_u_location', 'sea_water_x_velocity', 
#            'eastward_sea_water_velocity'
            ],
        long_names = "Sea water velocity along X at U location", 
        units = "m s-1", 
    ),
    v3d = dict(
        names=['vz', 'v3d'],
        standard_names=['sea_water_y_velocity_at_v_location', 'sea_water_y_velocity', 
            'northward_sea_water_velocity'],
        long_names = "Sea water velocity along Y at V location", 
        units = "m s-1", 
    ),
    w3d = dict(
        names=['wz', 'w3d'],
        standard_names=['sea_water_z_velocity_at_w_location', 'sea_water_z_velocity'],
        long_names = "Sea water velocity along Z at W location", 
        units = "m s-1", 
    ),
    ubt = dict(
        names=['u','ubt', 'u2d'],
        standard_names=['barotropic_sea_water_x_velocity_at_u_location'],
        long_names = "Sea water barotropic velocity along X at U location", 
        units = "m s-1", 
    ),
    vbt = dict(
        names=['v','vbt', 'v2d'],
        standard_names=['barotropic_sea_water_y_velocity_at_v_location'],
        long_names = "Sea water barotropic velocity along Y at V location", 
        units = "m s-1", 
    ),
    usurf = dict(
        names = ['usurf'], 
        standard_names=['sea_surface_x_velocity_at_u_location', 'sea_surface_x_velocity'],
        long_names = "Sea surface velocity along X at U location", 
        units = "m s-1", 
    ), 
    vsurf = dict(
        names = ['vsurf'], 
        standard_names=['sea_surface_y_velocity_at_y_location', 'sea_surface_y_velocity'],
        long_names = "Sea surface velocity along Y at V location", 
        units = "m s-1", 
    ), 
    ugbt = dict(
        names=['ugbt'],
        standard_names=['barotropic_sea_water_x_geostrophic_velocity_at_u_location'],
        long_names = "Sea water barotropic geostrophic velocity along X at U location", 
        units = "m s-1", 
    ),
    vgbt = dict(
        names=['vgbt'],
        standard_names=['barotropic_sea_water_y_geostrophic_velocity_at_v_location'],
        long_names = "Sea water barotropic geostrophic velocity along Y at V location", 
        units = "m s-1", 
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

    # Thermodynamics
    mld = dict(
        names = [], 
        standard_names = 'mixed_layer_depth', 
        long_names = "Mixed layer depth", 
        units = "m", 
    ), 
    ped = dict(
        names = [], 
        standard_names = 'potential_energy_deficit', 
        long_names = "Potential energy deficit", 
        units = "J m-2", 
    ), 


    # Bathymetry
    bathy = dict(
        names=['bathy', 'h0'], 
        standard_names=['model_sea_floor_depth_below_sea_level', 'model_sea_floor_depth_below_geoid', "sea_floor_depth_below_geoid" ], 
        long_names = 'Bathymetry', 
        units = 'm', 
    ), 
    bathy_u = dict(
        names=['bathy_u', 'hx'], 
        standard_names=['model_sea_floor_depth_below_sea_level_at_u_location', 'model_sea_floor_depth_below_geoid_at_u_location', "sea_floor_depth_below_geoid_at_u_location" ], 
        long_names = 'Bathymetry at U location', 
        units = 'm', 
        
    ), 
    bathy_v = dict(
        names=['bathy_v', 'hy'], 
        standard_names=['model_sea_floor_depth_below_sea_level_at_v_location', 'model_sea_floor_depth_below_geoid_at_v_location' "sea_floor_depth_below_geoid_at_v_location" ], 
        long_names = 'Bathymetry at V location', 
        units = 'm', 
    ), 
    depth = dict(
        names = ['depth', 'dep', 'depthu', 'depthv'], 
        standard_names = ['ocean_layer_depth', 'ocean_layer_depth_at_t_location', 'ocean_layer_depth_at_u_location', 
            'ocean_layer_depth_at_v_location'], 
        long_names =  'Depth', 
        units = 'm', 
    ), 
    depth_u = dict(
        names = ['depth_u', 'depthu'], 
        standard_names = 'ocean_layer_depth_at_u_location', 
        long_names =  'Depth at U location', 
        units = 'm', 
    ), 
    depth_v = dict(
        names = ['depthv', 'depth_v'], 
        standard_names = 'ocean_layer_depth_at_v_location', 
        long_names =  'Depth at V location', 
        units = 'm', 
    ), 
    depth_w = dict(
        names = ['depthw', 'depth_w'], 
        standard_names = 'ocean_interface_depth_at_w_location', 
        long_names =  'Depth of interfaces at W location', 
        units = 'm', 
    ), 
        
    # Cell sizes
    dx = dict(
        names = ['dx', 'dx_u', 'dx_v', 'dx_f'],
        standard_names = 'cell_x_size',
        long_names = "Mesh size along x",
        units = 'm',
    ),
    dx_u = dict(
        names = ['dx_u'],
        standard_names = 'cell_x_size_at_u_location',
        long_names = "Mesh size along x at u location",
        units = 'm',
    ),
    dx_v = dict(
        names = ['dx_v'],
        standard_names = 'cell_x_size_at_v_location',
        long_names = "Mesh size along x at v location",
        units = 'm',
    ),
    dx_f = dict(
        names = ['dx_f'],
        standard_names = 'cell_x_size_at_f_location',
        long_names = "Mesh size along x at f location",
        units = 'm',
    ),
    dy = dict(
        names = ['dy', 'dy_u', 'dy_v', 'dy_f'],
        standard_names = 'cell_y_size',
        long_names = "Mesh size along y",
        units = 'm',
    ),
    dy_u = dict(
        names = ['dy_u'],
        standard_names = 'cell_y_size_at_u_location',
        long_names = "Mesh size along y at u location",
        units = 'm',
    ),
    dy_v = dict(
        names = ['dy_v'],
        standard_names = 'cell_y_size_at_v_location',
        long_names = "Mesh size along y at v location",
        units = 'm',
    ),
    dy_f = dict(
        names = ['dy_f'],
        standard_names = 'cell_y_size_at_f_location',
        long_names = "Mesh size along y at f location",
        units = 'm',
    ),
    dz = dict(
        names = ['dz'], 
        standard_names = 'ocean_layer_thickness', 
        long_names = "Ocean layer thickness", 
        units = "m", 
    ), 
    dz_u = dict(
        names = ['dz_u'], 
        standard_names = 'ocean_layer_thickness_at_u_location', 
        long_names = "Ocean layer thickness at U location", 
        units = "m", 
    ), 
    dz_v = dict(
        names = ['dz_v'], 
        standard_names = 'ocean_layer_thickness_at_v_location', 
        long_names = "Ocean layer thickness at Vlocation", 
        units = "m", 
    ), 
    dz_w = dict(
        names = ['dz_w'], 
        standard_names = 'ocean_layer_thickness_at_w_location', 
        long_names = "Ocean layer thickness at W location", 
        units = "m", 
    ), 
    dlon = dict(
        names = ['dlon', 'dlon_u', 'dlon_v', 'dlon_f'],
        standard_names = 'cell_x_size',
        long_names = "Mesh size along x",
        units = 'degrees',
    ),
    dlon_u = dict(
        names = ['dlon_u'],
        standard_names = 'cell_x_step_at_u_location',
        long_names = "Mesh step along x at u location",
        units = 'degrees',
    ),
    dlon_v = dict(
        names = ['dlon_v'],
        standard_names = 'cell_x_step_at_v_location',
        long_names = "Mesh step along x at v location",
        units = 'degrees',
    ),
    dlon_f = dict(
        names = ['dlon_f'],
        standard_names = 'cell_x_step_at_f_location',
        long_names = "Mesh step along x at f location",
        units = 'degrees',
    ),
    dlat = dict(
        names = ['dlat', 'dlat_u', 'dlat_v', 'dlat_f'],
        standard_names = 'cell_y_step',
        long_names = "Mesh step along y",
        units = 'degrees',
    ),
    dlat_u = dict(
        names = ['dlat_u'],
        standard_names = 'cell_y_step_at_u_location',
        long_names = "Mesh step along y at u location",
        units = 'degrees',
    ),
    dlat_v = dict(
        names = ['dlat_v'],
        standard_names = 'cell_y_step_at_v_location',
        long_names = "Mesh step along y at v location",
        units = 'degrees',
    ),
    dlat_f = dict(
        names = ['dlat_f'],
        standard_names = 'cell_y_step_at_f_location',
        long_names = "Mesh step along y at f location",
        units = 'degrees',
    ),
   
    # Cell volumes
    vol = dict(
        names = ['vol'], 
        standard_names = 'cell_volume', 
        long_names = "Volume of the cell", 
        units = "m3", 
    ), 
    vol_u = dict(
        names = ['vol_u'], 
        standard_names = 'cell_volume_at_u_location', 
        long_names = "Volume of the cell at U location", 
        units = "m3", 
    ), 
    vol_v = dict(
        names = ['vol_v'], 
        standard_names = 'cell_volume_at_v_location', 
        long_names = "Volume of the cell at Vlocation", 
        units = "m3", 
    ), 
    vol_w = dict(
        names = ['vol_w'], 
        standard_names = 'cell_volume_at_w_location', 
        long_names = "Volume of the cell at W location", 
        units = "m3", 
    ), 
   
    # Coriolis
    f0 = dict(
        names = [], 
        standard_names = "coriolis_parameter", 
        long_names = "Coriolis parameter", 
        units = "s-1", 
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
    u10m = dict(
        names = [], 
        standard_names = ["eastward_wind","x_wind","x_wind_at_10m","x_wind_at_u_location","x_wind_at_10m_at_u_location"], 
        long_names = "10-m zonal wind speed (westerly)", 
        units = "m s-1", 
    ), 
    v10m = dict(
        names = [], 
        standard_names = ["northward_wind","y_wind","y_wind_at_10m","y_wind_at_v_location","y_wind_at_10m_at_v_location"], 
        long_names = "10-m meridional wind speed (northerly)", 
        units = "m s-1", 
    ), 
    ux10m = dict(
        names = [], 
        standard_names = ["x_wind", "grid_eastward_wind",
                          "x_wind_at_10m","x_wind_at_u_location","x_wind_at_10m_at_u_location"], 
        long_names = "10-m wind speed along X", 
        units = "m s-1", 
    ), 
    vy10m = dict(
        names = [], 
        standard_names = ["y_wind", "grid_northward_wind",
                          "y_wind_at_10m","y_wind_at_v_location","y_wind_at_10m_at_v_location"], 
        long_names = "10-m wind speed along Y", 
        units = "m s-1", 
    ), 

   # Ocean Atmosphere interface
   tauu = dict(
        names = [], 
        standard_names = ["surface_downward_eastward_stress",  
            "surface_eastward_stress"], 
        long_names = "Surface eastward wind stress", 
        units = "N m-2", 
    ), 
    tauv = dict(
        names = [], 
        standard_names = ["surface_downward_northward_stress",
            "surface_northward_stress"], 
        long_names = "Surface northward wind stress", 
        units = "N m-2", 
    ), 
    taux = dict(
        names = ['ustress'], 
        standard_names = ["surface_downward_x_stress", "surface_x_stress",
                          "surface_downward_x_stress_at_u_location"], 
        long_names = "Surface wind stress along X", 
        units = "N m-2", 
    ), 
    tauy = dict(
        names = ['vstress'], 
        standard_names = ["surface_downward_y_stress", "surface_y_stress",
                          "surface_downward_y_stress_at_v_location"], 
        long_names = "Surface wind stress along Y", 
        units = "N m-2", 
    ), 
    
    # Surfaces waves
    hs = dict(
        names = ['hs'], 
        standard_names = ["significant_height_of_wind_and_swell_waves"], 
        long_names = "Significant height of wind and swell waves", 
        units = "m", 
    ), 
    

)

#: List of generic variable names
generic_var_names = var_specs.keys()


#: Specifications for axes
axis_specs = dict(

    # Axes
    time = dict(
        names = ['time'], 
        standard_names = ['time'], 
        long_names = 'Time', 
        axis = 'T', 
    ),
    lon = dict(
        names = ['lon', 'longitude'], 
        standard_names = ['longitude'], 
        long_names = 'Longitude', 
        units = ['degrees_east', 'degree_east', 'degree_e', 'degrees_e', 'degreee', 'degreese'], 
        axis = 'X', 
        iaxis = 'ni', 
        jaxis = 'nj', 
    ), 
    lat = dict(
        names = ['lat', 'latitude'], 
        standard_names = ['latitude'], 
        long_names = 'Latitude', 
        units = ['degrees_north', 'degree_north', 'degree_n', 'degrees_n', 'degreen', 'degreesn'], 
        axis = 'Y', 
        iaxis = 'ni', 
        jaxis = 'nj', 
    ), 
    lon_u = dict(
        names = ['lon_u',], 
        standard_names = ['longitude_at_u_location'], 
        long_names = 'Longitude at U location', 
        units = ['degrees_east', 'degree_east', 'degree_e', 'degrees_e', 'degreee', 'degreese'], 
        axis = 'X', 
        iaxis = 'ni_u', 
        jaxis = 'nj_u', 
    ), 
    lat_u = dict(
        names = ['lat_v',], 
        standard_names = ['latitude_at_u_location'], 
        long_names = 'Latitude at V location', 
        units = ['degrees_north', 'degree_north', 'degree_n', 'degrees_n', 'degreen', 'degreesn'], 
        axis = 'Y', 
        iaxis = 'ni_u', 
        jaxis = 'nj_u', 
    ), 
    lon_v = dict(
        names = ['lon_v',], 
        standard_names = ['longitude_at_v_location'], 
        long_names = 'Longitude at V location', 
        units = ['degrees_east', 'degree_east', 'degree_e', 'degrees_e', 'degreee', 'degreese'], 
        axis = 'X', 
        iaxis = 'ni_v', 
        jaxis = 'nj_v', 
    ), 
    lat_v = dict(
        names = ['lat_v',], 
        standard_names = ['latitude_at_v_location'], 
        long_names = 'Latitude at V location', 
        units = ['degrees_north', 'degree_north', 'degree_n', 'degrees_n', 'degreen', 'degreesn'], 
        axis = 'Y', 
        iaxis = 'ni_v', 
        jaxis = 'nj_v', 
    ), 
    depth = dict(
        inherit = 'depth', 
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
#: List of generic axis names
generic_axis_names = axis_specs.keys()

#: List of all generic names (axes and variables)
generic_names= generic_var_names+generic_axis_names

#: Specifications for grid formating
grid_specs = {
    't': dict(lon='lon', lat='lat'), 
    'u': dict(lon='lon_u', lat='lat_u'), 
    'v': dict(lon='lon_v', lat='lat_v'), 
    'f': dict(lon='lon_f', lat='lat_f'), 
}
grid_specs['r'] = grid_specs['t']

def cp_suffix(idref, id, suffixes=['_u', '_v', '_w', '_z']):
    """Copy a suffix if found in an id to another id"""
    if isinstance(suffixes, basestring): suffixes = [suffixes]
    m = re.match('.+(%s)$'%'|'.join(suffixes), idref)
    if m is None: return id
    return id+m.group(1)

# Format specifications
for all_specs in var_specs, axis_specs:
    
    for name, specs in all_specs.items():
        
        # Always lists (except for dict and some keys)
        for key, value in specs.items():
            if isinstance(value, dict): continue
            if key in ['axis', 'inherit']: continue
            if isinstance(value, tuple):
                specs[key] = list(specs[key])
            elif not isinstance(value, list):
                specs[key] = [value]
                
        # Geo axes (variables only)
        if all_specs is var_specs:
            if 'axes' not in specs or not specs['axes']:
                specs['axes'] = {}
            specs['axes'].setdefault('t', 'time')
            suffixes = [('_'+s) for s in 'rftuv']
            specs['axes'].setdefault('y', cp_suffix(name, 'lat', suffixes=suffixes))
            specs['axes'].setdefault('x', cp_suffix(name, 'lon', suffixes=suffixes))
                    
        # Inherits from other specs
        if 'inherit' in specs:
            objname = specs['inherit']
            obj_specs = None
            if key==objname: # same name
                if key in generic_var_names and objname in generic_axis_names:
                    obj_specs = generic_axis_names
                elif key in generic_axis_names and objname in generic_var_names:
                    obj_specs = generic_var_names
            elif objname in generic_var_names:
                obj_specs = var_specs 
            else:
                obj_specs = axis_specs
            if obj_specs is not None:
                all_specs[name] = specs = dict_merge(specs, var_specs[objname])
                
        # No geo axes for axes!
        if all_specs is axis_specs and 'axes' in specs:
            del specs['axes']
        
        # Key = first id
        if name in specs['names']:
            specs['names'].remove(name)
        specs['names'].insert(0, name)

del specs

def cf2search(name, mode=None):
    """Extract specs from :attr:`axis_specs` or :attr:`var_specs` to form a search dictionary
    
    :Params:
    
        - **name**: Generic name of an axis or a variable.
        - **mode**, optional: Search mode [default: None->``"ns"``].
          A string containg one or more of the follwing letters:
        
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
    if name in var_specs:
        specs = var_specs[name]
    elif name in axis_specs:
        specs = axis_specs[name]
    else:
        VACUMMError("Wrong generic name. It should be one of: "+' '.join(generic_axis_names+generic_var_names))
    
    # Form search dict
    if not isinstance(mode, basestring): mode = 'nsa'
    keys = []
    for m in mode:
        for key in ['names', 'standard_names', 'axis', 'long_names', 'units']:
            if key.startswith(m):
                keys.append(key)
                break
    return OrderedDict((k, specs[k]) for k in keys if k in specs)

def cf2atts(name):
    """Extract specs from :attr:`axis_specs` or :attr:`var_specs` to form
    a dictionary of atrributes (units and long_name)"""
    # Get specs
    if name in var_specs:
        specs = var_specs[name]
    elif name in axis_specs:
        specs = axis_specs[name]
    else:
        VACUMMError("Wrong generic name. It should be one of: "+' '.join(generic_axis_names+generic_var_names))
    
    # Form search dict
    return dict((k, specs[k]) for k in ['long_names', 'units', 'axis'] if k in specs)



# Format a variable
def format_var(var, name, force=True, format_axes=True, order=None, **kwargs):
    """Format a MV2 variable according to its generic name
    
    
    :Params:
    
        - **var**: A :mod:`numpy` or :mod:`MV2` variabe.
        - **name**: Generic name of variable. It should be one of
          those listed by :attr:`generic_var_names`.
        - **force**: Overwrite attributes in all cases.
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
    var = MV2.asarray(var)
    
    # Check specs
    if name not in generic_var_names:
        raise KeyError("Generic var name not found '%s'. Please choose one of: "%(
            name, ', '.join(generic_var_names)))
    specs = var_specs[name]
    # - merge kwargs and specs
    for key, val in kwargs.items():
        if val is None or key not in specs: continue
        # Check type
        if not isinstance(val, list) and isinstance(specs[key], list): 
            val = [val]
        # Set
        specs[key] = val
        del kwargs[key]
    # - id
    if force or var.id.startswith('variable_'):
        var.id = specs['names'][0]
    # - other specs
    for att, val in specs.items():
        if att in ['long_names', 'standard_names']:
            att = att[:-1]
        elif att in ['names', 'axes']: continue
        if force or not getattr(var, att, ''):
            setattr(var, att, val[0])
    # - user attributes
    for att, val in kwargs.items():
        if force or not getattr(var, att, ''):
            setattr(var, att, val)
    
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
        axspecs = specs['axes']
        for key, meth in axismeths.items():
            axis = getattr(var, meth)()
            if order is not None: order.replace(key, '-')
            if axis is not None:
                format_axis(axis, axspecs[key], **kwaxes[key])
        
        # Check remaining simple axes (DOES NOT WORK FOR 2D AXES)
        if order is not None and order!='-'*len(order):
            for key in axismeths.keys():
                if key in order:
                    format_axis(var.getAxis(order.index(key)), axspecs[key], **kwaxes[key])
            
    
    return var

# Format an axis
def format_axis(axis, name, force=True, recreate=False, format_subaxes=True, **kwargs):
    """Format a MV2 axis according to its generic name
    
    
    :Params:
    
        - **var**: A :mod:`numpy` or :mod:`MV2` variabe.
        - **name**: Generic name of variable. It should be one of
          those listed by :attr:`generic_var_names`.
        - **force**, optional: Overwrite attributes in all cases.
        - **axes2d_<param>**, optional: <param> is passed to 
          :func:`~vacumm.misc.grid.misc.create_axes2d` for 2D axes.
        - **recreate**, optional: Force recreation of axes using either
          
        - Other parameters are passed as attributes.
    
    :Examples:
    
        >>> axis1d = format_axis(array1d, 'lon')
        >>> axis2d = format_axis(array2d, 'lonu', axes2d_iid = 'xi',  axes2d_jid = 'xi', recreate=True)
        
        
    """
    
    # Check specs
    if name not in generic_axis_names:
        raise KeyError("Generic axis name not found '%s'. Please choose one of: %s"%(
            name, ', '.join(generic_axis_names)))
    specs = axis_specs[name]
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
            axis = cdms2.create_axis(axis)
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
    # - id
    if force or axis.id.startswith('variable_') or axis.id.startswith('axis_'):
        axis.id = specs['names'][0]
    # - other specs
    for att, val in specs.items():
        if att in ['long_names', 'standard_names']:
            att = att[:-1]
        elif att in ['names', 'iaxis', 'jaxis'] or (att=='axis' and axis2d): continue
        if force or not getattr(axis, att, ''):
            setattr(axis, att, val[0])
    # - user attributes
    for att, val in kwargs.items():
        if force or not getattr(axis, att, ''):
            if isinstance(val, list): val = val[0]
            setattr(axis, att, val)
    
    # Sub-axes for 2D axes
    if axis2d and format_subaxes:
        format_axis(axis.getAxis(-1), specs['iaxis'][0])
        format_axis(axis.getAxis(-2), specs['jaxis'][0])
    
    return axis
    
def format_grid(grid, pt, **kwargs):
    """Format a grid and its axes"""
    if cdms2.isVariable(grid): grid = grid.getGrid()
    if grid is None: return
    gs = grid_specs[pt]
    lon = grid.getLongitude()
    lat = grid.getLatitude()
    format_axis(lon, gs['lon'])
    format_axis(lat, gs['lat'])

def match_var(var, name, mode=None, **kwargs):
    """Check if a variable match some specifications"""
    search = cf2search(name, mode=mode)
    search.update(kwargs)
    return ncmatch_obj(var, **search)

