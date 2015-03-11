#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright or © or Copr. Actimar/IFREMER (2010-2015)
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


__author__ = 'Jonathan Wilkins'
__email__ = 'wilkins@actimar.fr'
__date__ = '2010-11-02'
__doc__ = 'Exceptions utilities'


import sys, traceback, types


class ExceptionDebugger(object):
    '''
    Usefull for debugging.

    Retreive details about exceptions like locals for each frames.

    This provides informations about the program state (variables) for each frame
    of the call stack.
    '''

    def __init__(self, exit=True, nframe=None):
        self.exit, self.nframe = exit, nframe

    def getDetailedExceptionInfo(self, etype=None, evalue=None, tb=None):
        '''
        Print the usual traceback information, followed by a (huge!) listing of all the
        local variables in each frame.

        **TODO**: use StringIO

        '''
        ls = u'-'*80
        info = ''
        try:
            if etype is None or evalue is None or tb is None:
                fmt_tb = traceback.format_exc()
            else:
                fmt_tb = u''.join(traceback.format_exception(etype, evalue, tb))
            if tb is None:
                tb = sys.exc_info()[2]
            if tb is None:
                return u'(Error while formating exception traceback, sys.exc_info returned no traceback object)'
            while 1:
                if not tb.tb_next:
                    break
                tb = tb.tb_next
            stack = []
            f = tb.tb_frame
            while f:
                stack.append(f)
                f = f.f_back
            stack.reverse()
            info = u'%s\nException details, locals by frame, innermost last'%(ls)
            if self.nframe is None:
                self.nframe = len(stack)
            for frame in stack[-self.nframe:]:
                #info = u'%s\n%s\nFrame %s in %s at line %s\n' % (info, ls, frame.f_code.co_name, frame.f_code.co_filename, frame.f_lineno)
                info = u'%s\n%s\n%s\n' % (info, ls, ''.join(traceback.format_stack(frame, 1)))
                for key in sorted(frame.f_locals.keys()):
                    value = frame.f_locals[key]
                    if type(value) == types.ModuleType or type(value) == types.FunctionType:
                        continue
                    if key.startswith('__') and key.endswith('__'):
                        continue
                    info = u'%s\t%20s '%(info, key)
                    try:
                        s = str(value)
                        if len(s) > 800:
                            s = s[:800] + ' (...)'
                        s = s.replace('\n', '\n\t'+' '*22)
                        s = '%s = %s'%(type(value), s)
                        info = u'%s%s\n'%(info, s)
                    except Exception, e:
                        info = u'%s<ERROR WHILE PRINTING VALUE: %s>\n'%(info, e)
        except Exception, e:
            info = u'%s<ERROR WHILE PRINTING EXCEPTION DETAILS: %s>\n'%(info, e)
        info = u'%s\n%s\n%s\n%s\n'%(info, ls, fmt_tb, ls)
        return info

    def exceptHook(self, etype, evalue, tb):
        '''Exception hook used when installed by :meth:`bindExceptHook`'''
        print>>sys.stderr, self.getDetailedExceptionInfo(etype, evalue, tb)
        if self.exit:
            sys.exit(1)

    def bindExceptHook(self):
        '''Install this classe's exception traceback feature to sys.excepthook.
        :meth:`getDetailedExceptionInfo` will then be called when an exception
        is not handled by user's code (and then before exiting).
        '''
        sys.excepthook = self.exceptHook


# Default exception debbuger
_excdbg = ExceptionDebugger()

def bindExceptHook(*a, **k):
    '''See :meth:`ExceptionDebugger.bindExceptHook`'''
    _excdbg.bindExceptHook(*a, **k)


def getDetailedExceptionInfo(*a, **k):
    '''See :meth:`ExceptionDebugger.getDetailedExceptionInfo`'''
    return _excdbg.getDetailedExceptionInfo(*a, **k)



