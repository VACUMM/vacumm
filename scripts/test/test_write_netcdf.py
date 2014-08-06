"""Testing to read/modidy/(re)write NEtCDF files without loosing attributes"""

# Inits
ncfile = "nobounds.nc"
outfile = "result.nc"

# Imports
from vcmq import cdms2, warnings, data_sample, code_dir_name, os
from vacumm.misc.misc import get_atts

# Warnings off
warnings.simplefilter('ignore')

# Open input file
f = cdms2.open(data_sample(ncfile))

# Open output file
outfile = os.path.join(code_dir_name(), outfile)
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
            data[:,:,:,80:200].mask=True
        
        
    except: # - Variable is a single float
        interm = f[var]
        # For single float, need to re-allocate attributes
        data=cdms2.createVariable(data,typecode='f',id=interm.id,attributes=get_atts(interm))
        
    # Write variable
    fw.write(data)
    
# -- Global attributes
# (add automatically Conventions = "CF-1.0" as global attribute)
ga = f.listglobal()
for att in ga:
    setattr(fw,att,f.getglobal(att))

# Close files
f.close()
fw.close()

# Test variables in output file
fn = cdms2.open(outfile)
fn.showvariable()
fn.close
    
