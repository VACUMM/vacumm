#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, os.path, pprint, sys
import vacumm.misc.config as C

def main(argv=None):
	
	if '--specfile' in argv:
		i = argv.index('--specfile')
		specfile = argv[i+1]
		argv = argv[:i]+argv[i+2:]
	else:
		specfile = '%s.ini'%(os.path.splitext(__file__)[0])
	
	if not '--cfgfile' in argv:
		cfgfile = '%s.cfg'%(specfile)
		argv += ['--cfgfile', cfgfile]
	
	cfgm = C.ConfigManager(specfile)
	
	print '***'
	print '*** Usage: %s --specfile specfile.ini --cfgfile cfgfile.cfg'%os.path.basename(sys.argv[0])
	print '***'
	print '*** You can use the --help flag with this script to interact with it'
	print '*** to test setting config values with command line options and/or'
	print '*** the config file'
	print '***'
	
	cfg = cfgm.defaults()
	print '\n*** Default config from specfile %r:\n\n%s\n'%(specfile, pprint.pformat(cfg.dict()))
	
	# TODO: fix load method which currently apply defaults, make it as an option
	#cfg = cfgm.load(cfgfile)
	import configobj
	cfg = configobj.ConfigObj(cfgfile)
	print '\n*** Content of the config file %r:\n\n%s\n'%(cfgfile, pprint.pformat(cfg.dict()))
	
	cfg = cfgm.opt_parse(args=argv)
	print '\n*** Command line config only:\n\n%s\n'%(pprint.pformat(cfg.dict()))
	
	cfg = cfgm.opt_parse(args=argv, cfgfilepatch=True)
	print '\n*** Merged config: file and command:\n\n%s\n'%(pprint.pformat(cfg.dict()))
	
	cfg = cfgm.opt_parse(args=argv, patch=True, cfgfilepatch=True)
	print '\n*** Merged config: default, file and command:\n\n%s\n'%(pprint.pformat(cfg.dict()))
	
	# TODO: because the load method applies defaults, the real order is command, default and file
	cfg = cfgm.opt_parse(args=argv, patch=True, cfgfilepatch='after')
	print '\n*** Merged config: default, command and file:\n\n%s\n'%(pprint.pformat(cfg.dict()))

if __name__ =='__main__':
	main(sys.argv[1:])

