#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright or © or Copr. Actimar/IFREMER (2010-2017)
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

