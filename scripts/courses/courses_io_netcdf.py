#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Testing to read/modidy/(re)write NEtCDF files without loosing attributes"""

# Inits
ncfile = "nobounds.nc"
outfile = "result.nc"

# Imports
from vcmq import *
from vacumm.misc.misc import get_atts

# Warnings off
warnings.simplefilter('ignore')

# Open input file
f = cdms2.open(data_sample(ncfile))

# Open output file
fw = cdms2.open(outfile,'w')

# Cancel the default behaviour (<=> generating bounds variables)
cdms2.setAutoBounds('off')

# -- Read and write in another file variables
for var in f.listvariable():
    print '>> '+var+ ' <<'
    data=f(var)
    try: # - Variable is an n-dimension array
        data.getValue()
        # - Example of moficiation - mask an area
        if var == 'TEMP':
            data[:,:,:,80:200]=MV2.masked
        # => Practice: Try to modify the variable and its mask before writing.
        # => Practice: Change attributes of the variable.
    except: # - Variable is a single float
        interm = f[var]
        # For single float, need to re-allocate attributes
        data=cdms2.createVariable(data,typecode='f',id=interm.id,attributes=get_atts(interm))
    data.foo = 'bar'
        
    # Write variable
    fw.write(data)
    
# -- Global attributes
# (add automatically Conventions = "CF-1.0" as global attribute)
ga = f.listglobal()
for att in ga:
    setattr(fw,att,f.getglobal(att))
fw.toto = 'tutu'
# => Practice: Add new global attributes to the file.

# Close files
f.close()
fw.close()

# Test variables in output file
fn = cdms2.open(outfile)
fn.showvariable()
fn.close

# -----------------------------------------------------------
# ---- Read large datasets ----
# - Example with previmer files
print 10*'-'+' How to use ncread_files ...'
rep='/home/oo4/oo/modeles_previmer/f1_e2500/best_estimate/2013'
xe = ncread_files(os.path.join(rep,"PREVIMER_F1-MARS3D-MANGAE2500_%Y%m%dT%H00Z.nc"), 'xe',  
    ('2013-01-05', '2013-01-07', 'co'))    
print xe.shape
# => Practice: Try different parameters in ncread_files.
    
