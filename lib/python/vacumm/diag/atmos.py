# -*- coding: utf8 -*-
import numpy as N, MV2
from vacumm.data.cf import format_var

def wind_stress(u, v, rhoa=1.25, cd=0.016, format_axes=False, alongxy=None):
    """Compute the sea surface zonal and meridional wind stress from 10m wind components
    
    Output variables are formatted using :func:`~vacumm.data.cf.format_var`.
    
    :Formula:
    
        .. math::
        
            \\tau_x = \\rho_a C_d (U_{10m}^2+V_{10m}^2)^{\\frac{1}{2}}U_{10m}
            
            \\tau_y = \\rho_a C_d (U_{10m}^2+V_{10m}^2)^{\\frac{1}{2}}V_{10m}
    
    :Params:
    
        - **u/v**: Wind at 10m above sea surface.
        - **rhoa**, optional: Air density (in kg.m-3).
        - **cd**,  optional: Drag coefficient.
        - **format_axes**, optional: Also format axes using 
          :func:`~vacumm.data.cf.format_axis`.
        - **alongxy**, optional: Format variables considering components are along X/Y
          direction and not along zonal/meridional direction.
    
    :Return: ``us, vs``
    """
    # Init and format variables
    us = MV2.asarray(u).clone()
    vs = MV2.asarray(v).clone()
    if alongxy is None:
        alongxy =  hasattr(u, 'long_name') and 'x' in u.long_name.lower()        
    if alongxy:
        usname, vsname = 'taux', 'tauy'
    else:
        usname, vsname = 'tauu', 'tauv'
    format_var(us, usname, format_axes=format_axes)
    format_var(vs, vsname, format_axes=format_axes)
    
    # Compute
    uvmod = N.ma.sqrt(u**2+v**2)
    us.assignValue(cd*rhoa*uvmod*u)
    vs.assignValue(cd*rhoa*uvmod*v)
    del uvmod
    return us, vs
    
def ws2w(us, vs, rhoa=1.25, cd=0.016, format_axes=False, alongxy=None):
    """Convert from wind stress to 10m wind components
    
    This function is the reverse one of :func:`wind_stress`.
    Output variables are formatted using :func:`~vacumm.data.cf.format_var`.
    
    :Formula:
    
        .. math::
        
            U_{10m} = (\\rho_a C_d)^{-\\frac{1}{2}} 
            (\\tau_x^2+\\tau_y^2)^{-\\frac{1}{4}} \\tau_x
            
            V_{10m} = (\\rho_a C_d)^{-\\frac{1}{2}} 
            (\\tau_x^2+\\tau_y^2)^{-\\frac{1}{4}} \\tau_y
            
    :Params:
    
        - **us/vs**: Wind stress components.
        - **rhoa**, optional: Air density (in kg.m-3).
        - **cd**,  optional: Drag coefficient.
        - **format_axes**, optional: Also format axes using 
          :func:`~vacumm.data.cf.format_axis`.
        - **alongxy**, optional: Format variables considering components are along X/Y
          direction and not along zonal/meridional direction.
    
    :Return: ``u, v``
    """
   # Init and format variables
    u = MV2.asarray(us).clone()
    v = MV2.asarray(vs).clone()
    if alongxy is None:
        alongxy =  hasattr(u, 'long_name') and 'x' in us.long_name.lower()        
    if alongxy:
        uname, vname = 'ux10m', 'vy10m'
    else:
        uname, vname = 'u10m', 'v10m'
    format_var(u, uname, format_axes=format_axes)
    format_var(v, vname, format_axes=format_axes)
    
    # Compute
    zero =  us.filled(1)==0.
    zero &= vs.filled(1)==0.
    uvsmod = (us**2+vs**2)**-0.25
    uvsmod /= N.sqrt(rhoa*cd)
    u.assignValue(uvsmod*us)
    v.assignValue(uvsmod*vs)
    u.assignValue(MV2.where(zero, 0., u))
    v.assignValue(MV2.where(zero, 0., v))
    del uvsmod
    return u, v
    
    
if __name__=='__main__':
    import MV2
    u = MV2.arange(3.)
    v = 2*u
    us, vs = wind_stress(u, v)
    U, V = ws2w(us, vs)
    print u, v
    print U, V
    
