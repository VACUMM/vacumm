from vacumm.data.satellite.satellite import Satellite
import ConfigParser

class Sst(Satellite) :
    """ SST satellite """
    def __init__(self,cfg=None):
        import os, cdtime
        Satellite.__init__(self)

	SCRIPT_DIR=os.getcwd()
        self.SCRIPT_DIR = SCRIPT_DIR

        #-- retrocompatibilite
        if cfg is None:
            config = ConfigParser.RawConfigParser()
            config.read(os.path.join(SCRIPT_DIR,'config.cfg'))	
            andeb = config.getint('Time Period', 'andeb')
            anfin = config.getint('Time Period', 'anfin')		
            mdeb = config.getint('Time Period', 'mdeb')
            mfin = config.getint('Time Period', 'mfin')
            jdeb = config.getint('Time Period', 'jdeb')
            jfin = config.getint('Time Period', 'jfin')
            hdeb = config.getint('Time Period', 'hdeb')
            hfin = config.getint('Time Period', 'hfin')
            
            self.WORKDIR = config.get('Env', 'workdir')
        else:
            
            andeb = cfg['Time Period']['andeb']
            anfin = cfg['Time Period']['anfin']
            mdeb = cfg['Time Period']['mdeb']
            mfin = cfg['Time Period']['mfin']
            jdeb = cfg['Time Period']['jdeb']
            jfin = cfg['Time Period']['jfin']
            hdeb = cfg['Time Period']['hdeb']
            hfin = cfg['Time Period']['hfin']

            self.WORKDIR = cfg['Env']['workdir']

        #print "%(andeb)s-%(mdeb)s-%(jdeb)s %(hdeb)s"%vars()
        #print "%(anfin)s-%(mfin)s-%(jfin)s %(hfin)s"%vars()

        # Conversion en "Component Time"
        self.ctdeb=cdtime.comptime(andeb,mdeb,jdeb,hdeb,0,0)
        self.ctfin=cdtime.comptime(anfin,mfin,jfin,hfin,0,0)


        DIR_SST=os.path.join(self.WORKDIR,'SST')    # repertoire de stockage des donnees     
        if os.path.isdir(DIR_SST)==False:
            os.mkdir(DIR_SST)
        self.WORKDIR=os.path.join(self.WORKDIR,'SST')
            
        self.name="Satellite Sea surface Temperature"
        self.shortname="SST"
        self.units=" "

        #print "Workdir in sst.py"
        #print self.WORKDIR
