
# Module de colocalisation des sorties de modele avec les observations
#**********************************************************************************
# MOD : c_coloc_obs_mod.py
# OBJ : Module de colocalisation base sur le module SCRIP
#
# CRE : G. Charria (Ifremer)
#       S. Theetten (Ifremer)
# VER : 1.0 (01/2010)
#**********************************************************************************
def lonlatobs_to_lonlatmodel(WORKDIR,DIR_REGRID,lon_obs,lat_obs,lon_model,lat_model,sst_nar):
    #******************************************************************************
    # Interpolation horizontale des observations sur la grille du modele
    # lon_obs,lat_obs => lon_model,lat_model
    #
    # Utilisation du module SCRIP
    #******************************************************************************
    import os
    import cdms, regrid
    os.chdir(DIR_REGRID)  # on se place dans le repertoire de travail

    scrip_exec='/export/home/logiciels/CDAT/installation/SCRIP/scrip'

    scrip_in="""
    &remap_inputs
        num_maps = 2
        grid1_file = '%s'
        grid2_file = '%s'
        interp_file1 = '%s'
        interp_file2 = '%s'
        map1_name = '%s Conservative Mapping'
        map2_name = '%s Conservative Mapping'
        map_method = 'conservative'
        normalize_opt = 'frac'
        output_opt = 'scrip'
        restrict_type = 'latitude'
        num_srch_bins = 90
        luse_grid1_area = .false.
        luse_grid2_area = .false.
    /
    
    """

    # ------------------------------------------------------------------------------------------
    def toScripFormat(slab,file):
        """ Dumps information about a slab into a file that can be used by the SCRIP regridder
        Usage
        toSCRIP(slab,file)
        where:
        file: is the output file to which information in SCRIP style will be dumped
        slab: is the slab from which informations are gathered
        """
        
        ## Get definitions of grid
        g=slab.getGrid()
        print 'lat',g.getMesh()[0,0]
        print 'lon',g.getMesh()[0,1]
        m=slab.mask()
        if m is not None:
            print m.shape
            g.setMask(m)
        print g.shape
        if g is None:
            raise 'Error, no grid defined'
        f=cdms.Cdunif.CdunifFile(file,'w')
        g.writeScrip(f,'test_scrip')
        f.close()
        return
    # -------------------------------------------------------------------------------------------



    #-------------------------------------------------------
    print ''
    print '---------- INTERPOLATION HORIZONTALE ----------'
    print ''
    #-------------------------------------------------------

    # Open the SCRIP remapping file and data file
    direc = ''
    fremap = cdms.open('rmp_NAR_to_MARS_conserv.nc')

    # Input data array: sst_nar

    infile=open(sys.argv[1])

    for l in infile.xreadlines(): ## Loop thru input file
        sp=l.split()
        f=cdms.open(sp[0]) # Obs file
        var=sp[1]
        acc=sp[2] # Accronym for obs
        finobs=acc+'_obs_SCRIP.nc'
        s=f(var)
        while s.rank()>2: ## Gets only lat/lon
            s=s[0]
        f.close()
        ## Creates the obs SCRIP File
        toScripFormat(s,finobs)
        
        qry = 'select accro_orig, version_orig from modelversions,models where models.id=modelversions.model'
        res=SQL.dbquery(qry)
        for m in res:
            ## Prepares strings
            rmp1='rmp_NAR_to_MARS_conserv.nc'
            rmp2='rmp_MARS_to_NAR_conserv.nc'
            tit1='NAR to MARS'
            tit2='MARS to NAR'
            fin='remap_grid_NAR.nc'
            finobs='remap_grid_MARS.nc'
            scrip_str=scrip_in % ( fin , finobs, rmp1, rmp2, tit1, tit2 )
            dbpth=os.path.join(DB.root,'grids',mod_accro)
            rmp1db=os.path.join(dbpth,rmp1)
            test = not os.path.exists(rmp1db)
            test = True
            if test:
                print 'Creating remap files between %s and %s' % (mod_accro,acc)
                ## Reads in grid file
                fnm=os.path.join(dbpth,mod_accro+'_grid.nc')
                f=cdms.open(fnm)
                sm=f('mask')[0]
                sm=MV.masked_equal(sm,0)
                toScripFormat(sm,fin)
                fscr=open('scrip_in','w')
                print >> fscr,scrip_str
                fscr.close()
                ln=os.popen(scrip_exec).readlines()
                print 'Moving remap files'
                ln=os.popen('mv %s %s' % (rmp1,rmp1db)).readlines()
                print ln
                os.popen('mv %s %s' % (rmp2,os.path.join(dbpth,rmp2))).readlines()
                os.remove('scrip_in')
                os.remove(fin)
        os.remove(finobs)
                
    print 'Done'
    ## Bellow lines to remap using remap files  

    print 'Remapping using : ',rmp1
    frmp1=cdms.open(rmp1)
    remapper=regrid.readRegridder(frmp1)
    sst_nar_int2D=remapper(s)
    ##     fout=cdms.open(acc+'_remapped_to_1x1.nc','w')
    ##     fout.write(s1)
    ##     fout.close()



    # Read the SCRIP regridder
    #regridf = regrid.readRegridder(fremap)

    # Mettre dans une boucle sur les differents pas de temps
    
    # Regrid the variable
    #sst_nar_int2D = regridf(sst_nar)

    return sst_nar_int2D
    #-- Fin interpolation horizontale
    #----------------------------------------------------



