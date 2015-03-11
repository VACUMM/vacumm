#!/usr/bin/env python
# -*- coding: utf8 -*-
# File: parcalls.py
# Date : 06/2008
# Authors: Guillaume SICOT, Nicolas THOMAS (Actimar)
# Desc: Contains tools to annotate functions with elapsed execution time, number of calls...
# Usefull for performance testing.
# Actimar: http://www.vacumm.fr

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
"""
Authors: Guillaume SICOT, Nicolas THOMAS (Actimar)
Tools to annotate functions with elapsed execution time, number of calls...
Usefull for performance testing.
Actimar: http://www.vacumm.fr

######################
#Usage example
######################

#Create a test function
def timeme():
    time.sleep(1)

#Annotate the function
timeme = timeit(timeme)

#Call the function
for i in range(5):
    timeme()

#Print annotations
print_functimes(timeme,no_call=True)

"""

import sys
import time
import resource


def timeit(f):
    """ Annotate a function with its elapsed execution time, number of calls...
    @param f : function to annotate
    @type f : Function
    @return : Annotated function
    @rtype : Function
    """

    #Create the annotated function
    def timed_f(*args, **kwargs):
        #Get current "real time" and "user processing" time
        t1 = time.time()
        user1 = __user_time__()
        #Run f
        try:
            res = f(*args, **kwargs)
        #Get new "real time" and "user processing" time
        finally:
            t2 = time.time()
            user2 = __user_time__()

        #Compute delta_t for "real time" and "user processing"
        dt = t2-t1
        duser = user2-user1

        #Annotate f :
        timed_f.fct_name=f.__name__
        timed_f.func_maxtime = max(dt, getattr(timed_f, 'func_maxtime', 0))
        timed_f.func_mintime = min(dt, getattr(timed_f, 'func_mintime',sys.maxint))
        timed_f.user_maxtime = max(duser, getattr(timed_f, 'user_maxtime', 0))
        timed_f.user_mintime = min(duser, getattr(timed_f, 'user_mintime',sys.maxint))
        timed_f.func_numcalls = getattr(timed_f, 'func_numcalls',0)+1
        timed_f.func_avgtime = ((dt+((timed_f.func_numcalls-1) * getattr(timed_f, 'func_avgtime',0))) / timed_f.func_numcalls)
        timed_f.user_avgtime = ((duser+((timed_f.func_numcalls-1) * getattr(timed_f, 'user_avgtime',0))) / timed_f.func_numcalls)

        return res

    #Init some annotations for the case where the function is never called
    function=timed_f
    function.func_numcalls = getattr(timed_f, 'func_numcalls',0)
    function.fct_name=f.__name__

    #Return the annotated function
    return function

def print_functimes(f,no_call=False):
    """ Write to stdout timed function annotations
    @param f : Annotated function
    @type f : Function
    @keyparam no_call : if "True" write  'Function not called' if f has not been called (else write nothing),[Default : False]
    """
    #Format header
    if f.func_numcalls>1:
        line_length = 36
    else:
        line_length = 32

    header='- fct : '+f.fct_name+' -'
    sep_line = '-'
    while len(sep_line)<len(header):
        sep_line+='-'
    while len(header)<line_length:
        header='-'+header+'-'
        sep_line=sep_line+'--'

    #Format results
    if f.func_numcalls==0:
        if no_call:
            print sep_line
            print header
            print 'Function not called'
            print sep_line
    elif f.func_numcalls==1:
        print sep_line
        print header
        print 'Running time: %.3f secs' % f.func_avgtime
        print 'Running user time: %.3f secs' % f.user_avgtime
        print sep_line
    else:
        print sep_line
        print header
        print 'Avg running time: %.3f secs' % f.func_avgtime
        print 'Max running time: %.3f secs' % f.func_maxtime
        print 'Min running time: %.3f secs' % f.func_mintime
        print 'Avg running user time: %.3f secs' % f.user_avgtime
        print 'Max running user time: %.3f secs' % f.user_maxtime
        print 'Min running user time: %.3f secs' % f.user_mintime
        print 'Number of calls: %d' % f.func_numcalls
        print sep_line

def __user_time__():
    """Get the current user processing time"""
    return resource.getrusage(resource.RUSAGE_SELF)[0]

# ######################
# #Usage example
# ######################
# #Create a test function
# def timeme():
#     time.sleep(1)

# #Annotate the function
# timeme = timeit(timeme)

# #Call the function
# for i in range(5):
#     timeme()

# #Print annotations
# print_functimes(timeme,no_call=True)

