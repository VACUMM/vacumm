#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import vacumm.misc.log as L

l1 = L.Logger(name='L1', level='info', format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
l2 = L.Logger(name='L2', level='debug', format='[%(asctime)s %(name)s %(levelname)s] %(message)s', date_format='%H:%M:%S')
l1.info('info from logger l1')
l1.debug('debug from logger l1')
l2.info('info from logger l2')
l2.debug('debug from logger l2')
l1.notice('logger l1 setting config from logger l2')
l1.config(l2)
l1.info('info from logger l1')
l1.debug('debug from logger l1')

lf = '%s.%s.log'%(os.path.splitext(sys.argv[0])[0], os.getpid())
l = L.Logger(
    name=os.path.basename(sys.argv[0]),
    level='debug',
    logfile=dict(
        filePath=lf, maxFileSize=1024*50, maxFileBkp=1,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        date_format='%Y-%m-%d %H:%M:%S %Z',
        level='info'
    ),
    console=dict(
        stream=sys.stdout, colorize=False,
        format='%(name)s - %(levelname)s - %(message)s',
        #level='verbose'
    )
)
l.info('logging to %s', lf)
l.debug('debug')
l.verbose('verbose')
l.warning('warning')
with file(lf) as lff:
    l.info('content of %s:\n%s', lf, lff.read(2**20))
if os.path.isfile(lf):
    l.debug('removing log file %s', lf)
    os.remove(lf)

