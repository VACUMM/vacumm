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

# ==================================================================================================

import os, re, StringIO, sys

# ==================================================================================================

__author__ = 'Jonathan Wilkins'
__email__ = 'wilkins@actimar.fr'
__date__ = '2011-03-01'
__doc__ = 'Load / save fortran namelist string / file.'


# ==================================================================================================

# FIXME:
# - [Check if still applicable] Fix badly handled comment when comment characters are in a quoted (string) value
# TODO:
# - Recheck namelist format(s?) and allow specifying differents characters for namelist start/end,
#   comments start and array indexing.
# - Allow N-dimensionnal arrays (1d currently).
# - Allow store of arrays into python lists rather than dicts.
# - Store namelist and variables in python lists (keep order)
# - Eventually store comments ?

# ==================================================================================================

class Namelist(dict):
    '''Handles the fortran namelist format as below:

        .. code-block:: none

            ! header comment

            &namelist1 b=.true. i=0,f=0.1, s="ab" , ss="abcd" /    ! inline comment

            ! standalone comment
            # another "comment"

            &namelist2
            integer=1500, float=1.5e3    ! inline comment
            str1="A : B"
            str2="/path/to/file.txt"
            array(0)=0
            array(1)=1
            /

        Empty lines, spaces and comments (from ":", "!" or "#" characters to the end of line) are allowed.

        Namelists start tag is the "&" character and end tag is "/", both on their own line or not.

        Variables are delimited with space(s) or comma.

        Namelists (sections) are stored and made accessible in a dictionnary fashion as this class
        inherit from dict.
            - Keys are the namelists names, values are the variables of the corresponding namelist.
            - Namelist variables are also stored in dicts.
            - If the variable is an array like variable, values are also stored in a dict.
        Variables types are automatically handled, available types are:
            - boolean: stored as bool python variable
            - integer: as int
            - real/scientific: as float
            - array: as dict

        .. note:
            - When loading namelists, previously loaded namelists are overriden, not cleared (use the parent class clear method if needed).
            - Only 1-dimensionnal array like variables are allowed

        .. note::
            - **Variables cannot be inserted directly into the Namelist**, you must access/setup the namelist level before (see examples)
              that's because there is no global variable in a namelist file !

        :Examples:

            >>> n = Namelist()
            >>> # The following is not correct and would results in a ValueError:
            >>> n['myvar'] = 0
            ... ValueError: Cannot set a global variable in namelist, value must be a dict
            >>> # That's the good way to do things:
            >>> n['mynamelist'] = {}
            >>> n['mynamelist']['myvar'] = 0

            >>> n = Namelist()
            >>> n['namelist1'] = dict(integer=1, real=1.0, string='hello', array=dict(a=0,b=1))
            >>> print n.save_string()
            ... &namelist1
            ...   real = 1.0
            ...   integer = 1
            ...   array(a) = 0
            ...   array(b) = 1
            ...   string = "hello"
            ...   /

            >>> n = Namelist.from_file('mynamelist.txt')


    '''

    _maxreadsize = 2 * (2**20) # n * Mio

    def __init__(self):
        dict.__init__(self)

    def __str__(self):
        return '\n\n'.join('%s:\n%s'%(namelist, '\n'.join(('  %s = %s (%s)'%(k,v,v.__class__.__name__) for k,v in variables.iteritems()))) for namelist, variables in self.iteritems())

    def __setitem__(self, key, value):
        if not isinstance(value, dict):
            raise ValueError('Cannot set a global variable in namelist, value must be a dict')
        dict.__setitem__(self, key, value)

    @staticmethod
    def remove_comments(s, strip=True):
        '''Internal method for comments filtering.
           Return s without fortran comments.
           If strip is True, also eliminate empty lines and leading and trailing spaces
        '''
        comment = r'^\s*[:|!|#].*$'
        #comment = r'''(?<!(?:"|')).*?[:|!|#].*(?<!(?:"|'))'''
        lines = (re.sub(comment, '', l) for l in s.splitlines())
        if strip: lines = (l.strip() for l in lines if l.strip())
        return '\n'.join(lines)

    @classmethod
    def from_file(cls, filepath, *a, **k):
        '''Factory method populating the created namelist with a file.'''
        o = cls(*a, **k)
        o.load_file(filepath)
        return o

    @classmethod
    def from_string(cls, string, *a, **k):
        '''Factory method populating the created namelist with a string.'''
        o = cls(*a, **k)
        o.load_string(string)
        return o

    def load_string(self, string):
        '''Load namelists from a string.'''
        #print '[start of raw namelist]\n%s\n[end of raw namelist]'%(s)
        # Clean comments / spaces / empty lines
        s = self.remove_comments(string)
        #print '[start of filtered namelist]\n%s\n[end of filtered namelist]'%(s)
        # Patterns
        namelist = re.compile(r'''\s*&(?P<namelist>[a-zA-Z][a-zA-Z0-9_]*)\s+(?P<statements>(?:(?:"|').*?(?:"|')|[^/])*)/''', re.MULTILINE|re.DOTALL|re.VERBOSE)
        variable = re.compile(r'\s*(?P<variable>[a-zA-Z][a-zA-Z0-9_]*)(?:\((?P<index>\d+)\))?\s*=\s*(?P<value>[^\s,]+)(?:\s+|,)', re.MULTILINE|re.DOTALL|re.VERBOSE)
        boolean = re.compile(r'\.?(?P<boolean>true|false|t|f)\.?', re.I)
        integer = re.compile(r'(?P<integer>[-+]?\d+)')
        real =  re.compile(r'(?P<real>[-+]?(\d+\.\d*|\.\d+)([eE][-+]?\d+)?)')
        string = re.compile(r'''(?:"|')(?P<string>.*?)(?:"|')''')
        # Values converters. Order is important !
        casts = ((string,str), (boolean,lambda v: v.lower().strip().startswith('t') and True or False), (real,float), (integer,int))
        # Loop on namelists
        for nm in namelist.finditer(s):
            ngd = nm.groupdict()
            #print (' namelist: %s '%(ngd['namelist'])).center(80, '=')
            #print 'start:', vm.start(), 'end:', vm.end()
            #print 'groups:', vm.groups(), 'groupdict:', vgd
            self[ngd['namelist']] = {}
            # Loop on variables
            for vm in variable.finditer(ngd['statements']):
                vgd = vm.groupdict()
                #print (' variable: %s '%(vgd['variable'])).center(80, '-')
                #print 'start:', vm.start(), 'end:', vm.end()
                #print 'groups:', vm.groups(), 'groupdict:', vgd
                index = vgd.get('index', None)
                # Initialize array (dict) if needed
                if index is not None and vgd['variable'] not in self[ngd['namelist']]:
                    self[ngd['namelist']][vgd['variable']] = {}
                # Find and convert value
                for p,c in casts:
                    m = p.search(vgd['value'])
                    if m:
                        if index is not None:
                            self[ngd['namelist']][vgd['variable']][index] = c(m.groups()[0])
                        else:
                            self[ngd['namelist']][vgd['variable']] = c(m.groups()[0])
                        break

    def save_string(self):
        '''Return the (fortran) formated namelists'''
        s = StringIO.StringIO()
        for namelist,variables in self.iteritems():
            s.write('&%s\n'%namelist)
            for vn,vv in variables.iteritems():
                if isinstance(vv, dict):
                    for i,vi in vv.iteritems():
                        if isinstance(vi, (str,unicode)):
                            s.write('  %s(%s) = "%s"\n'%(vn,i,vi))
                        elif isinstance(vi, bool):
                            s.write('  %s(%s) = .%s.\n'%(vn,i,str(vi).lower()))
                        else:
                            s.write('  %s(%s) = %s\n'%(vn,i,vi))
                else:
                    if isinstance(vv, (str,unicode)):
                        s.write('  %s = "%s"\n'%(vn,vv))
                    elif isinstance(vv, bool):
                        s.write('  %s = .%s.\n'%(vn,str(vv).lower()))
                    else:
                        s.write('  %s = %s\n'%(vn,vv))
            s.write('/\n\n')
        r = s.getvalue()
        s.close()
        return r

    def load_file(self, filepath):
        '''Load namelists from a text file'''
        self._filepath = filepath
        fh = file(self._filepath, 'rU')
        fc = fh.read(self._maxreadsize)
        fh.close()
        self.load_string(fc)

    def save_file(self, filename, append=False):
        '''Save the (fortran) formatted namelists to a file'''
        fh = file(filename, append and 'a' or 'w')
        fh.write(self.save_string())
        fh.close()


# ==================================================================================================

if __name__ == '__main__':
    '''Test this module'''

    namelist = '''\
! header comment

    &namelist1 b=.true. i=0,f=0.1, s="ab" , ss="abcd" /    ! inline comment

  # standalone &n a="comment" /

&namelist2
  integer=1500, float=1.5e3    ! inline comment
  str1 = "A : B"
  str2="/path/to/file.txt"
  array(0)=0
  array(1) = 1 /

 &namfine
  head_fine = './IN/head.quib3'/

    '''

    namelist = '''\
! header comment

    &namelist1 b=.true. i=0,f=0.1, s="ab" , ss="abcd" /    ! inline comment

&namelist2
  integer=1500, float=1.5e3    ! inline comment
  str1 = "A:B"
  str2='./path/to/file.txt'
  array(0)=0
  array(1) = 1 /

'''

    if len(sys.argv) == 2 and sys.argv[1] in ('-t', '-test', '--test'):
        print ' namelist in input '.center(80, '=')
        print namelist
        print
        nf = Namelist.from_string(namelist)
    elif len(sys.argv) == 2:
        nf = Namelist.from_file(sys.argv[1])
    else:
        print 'usage: %(prog)s (namelistfile|--test)'%dict(prog=os.path.split(sys.argv[0])[1])
        sys.exit(1)
    print ' namelist as a descriptive string (__str__ method) '.center(80, '=')
    print nf
    print
    print ' namelist as formated string (save_string method) '.center(80, '=')
    print nf.save_string()
    print



