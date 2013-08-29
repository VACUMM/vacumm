#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, os.path as P, pprint, sys
import vacumm.misc.config as C

def main(argv=None):
    specfile = '%s.ini'%(P.splitext(__file__)[0])
    cfgfile = '%s.cfg'%(P.splitext(__file__)[0])
    
    if not '--cfgfile' in argv:
        argv += ['--cfgfile', cfgfile]
    
    cfgm = C.ConfigManager(specfile)
    
    print '***'
    print '*** You can use the --help flag with this script to interact with it'
    print '*** to test setting config values with command line options and/or'
    print '*** the config file'
    print '***'
    
    print '\n### Using the config manager\n'
    
    cfg = cfgm.defaults()
    print '\n*** Default config from specfile %r:\n\n%s\n'%(specfile, pprint.pformat(cfg.dict()))
    
    # FIXME: load method which currently apply defaults, make it as an option ?
    #cfg = cfgm.load(cfgfile)
    # The ConfigManager currently set defaults from its specifications when using its load method.
    # So in this example we'll directly use configobj to show you only the content of the configuration file
    import configobj
    cfg = configobj.ConfigObj(cfgfile)
    print '\n*** Content of the config file %r:\n\n%s\n'%(cfgfile, pprint.pformat(cfg.dict()))
    
    optcfg = argcfg = None
    
    print '\n### Using the config manager with optparse\n'
    
    try:
    
        cfg = optcfg = cfgm.opt_parse(args=argv)
        print '\n*** Command line config only:\n\n%s\n'%(pprint.pformat(cfg.dict()))
        
        cfg = cfgm.opt_parse(args=argv, cfgfilepatch=True)
        print '\n*** Merged config: file and command:\n\n%s\n'%(pprint.pformat(cfg.dict()))
        
        cfg = cfgm.opt_parse(args=argv, patch=True, cfgfilepatch=True)
        print '\n*** Merged config: default, file and command:\n\n%s\n'%(pprint.pformat(cfg.dict()))
        
        # FIXME: because the load method applies defaults, the real order is command, default and file
        #cfg = cfgm.opt_parse(args=argv, patch=True, cfgfilepatch='after')
        #print '\n*** Merged config: default, command and file:\n\n%s\n'%(pprint.pformat(cfg.dict()))
    
    except SystemExit: pass
    
    print '\n### Using the config manager with argparse\n'
    
    try:
    
        cfg = argcfg = cfgm.arg_parse(args=argv)
        print '\n*** Command line config only:\n\n%s\n'%(pprint.pformat(cfg.dict()))
        
        cfg = cfgm.arg_parse(args=argv, cfgfilepatch=True)
        print '\n*** Merged config: file and command:\n\n%s\n'%(pprint.pformat(cfg.dict()))
        
        cfg = cfgm.arg_parse(args=argv, patch=True, cfgfilepatch=True)
        print '\n*** Merged config: default, file and command:\n\n%s\n'%(pprint.pformat(cfg.dict()))
    
    except SystemExit: pass
    
    assert optcfg == argcfg, 'optparse and argparse configurations differ !'

if __name__ =='__main__':
    main(sys.argv[1:])

