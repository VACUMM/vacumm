# -*- coding: utf8 -*-
"""Docstring refomatting utilities for VACUMM help

.. warning:: This module requires python 2.7+ to work.
"""
# Copyright or © or Copr. Actimar/IFREMER (2012-2015)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is a computer program whose purpose is to [describe
# functionalities and technical features of your software].
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


from textwrap import dedent
import re, inspect

class Docstring2Params(dict):
    """Scan the docstring of an object to load parameters names and descriptions
    and reformat them.


    The docstring must follow this example::

        '''Description

        Some text

        :Section: text

        :Params:

            - **arg**: text
              text
            - **opt**, optional: text

        Some text or sections.

        '''

    In particular:

        - Parameters are described as in this example, on one or more lines.
        - They must be inside a section whose name contains "param".

    :Params:

        - **obj**: Typically a function, method or class.
        - **select**, optional: A single or a list of parameter names to
          restrict scan.

    :Example:

        >>> dsp = Docstring2Params(myfunc)
        >>> dsp.format('firstparam')
        >>> params = dsp.asdict(select=['secondpar', 'thirdpar'], indent=8)
    """

    # Regular expressions
    RE_SECTION_MATCH = re.compile(r'^:([^:]+):.*').match
    RE_INDENT_MATCH = re.compile(r"^(\s*)(\S.*)?$").match
    RE_PARAM = re.compile(r"^(\s+)-\s*\*\*?\s*([^\*]+)\s*\*\*?([^:]*:.*)$")

    def __init__(self, obj, select=None, prefix=None, sections = ['param']):

        dict.__init__(self)

        # Format docstring
        if isinstance(obj, basestring): obj = eval(obj)
        self._target = obj
        self._docstring = doc = obj.__doc__
        if doc is None:
            self.params = {}
            return
        nn = doc.index("\n")
        doc = doc[nn+1:]
        doc = dedent(doc)
        docs = doc.splitlines()
        if select is not None and isinstance(select, basestring):
            select = [select]
        if prefix is None: prefix = ''

        # Loop on lines
        param = 0
        ifc = None
        for iline, line in enumerate(docs):

            # Section header
            m = self.RE_SECTION_MATCH(line)
            if m:
                if True in [secname in m.group(1).lower() for secname in sections]: # Params section header
                    param = 1
                    continue
                elif param: # New section, so the end
                    break

            # Still no param section found
            if not param: continue

            # Param declaration
            mp = self.RE_PARAM.match(line)
            if mp:
                param = mp.group(2)
                if select is None or param in select:
                    if prefix:
                        line = self.RE_PARAM.sub(r'\1- **%s\2**\3'%prefix, line)
                    if not (prefix+param) in self:
                        self[prefix+param] = []
                    self[prefix+param].append(line)
                if ifc is None: ifc = len(mp.group(1)) # Official indent
                continue

            # Valid line for a param description continuation
            if not isinstance(param, basestring): continue
            mi = self.RE_INDENT_MATCH(line)
            ns = len(mi.group(1))
            rem = mi.group(2)
            if rem is None or ns >= ifc:
                if ns == ifc and rem is not None:
                    if rem.startswith('-'): # Description without parameter
                        param = 2
                    else: # No longer describing a parameter
                        param = 1
                elif param != 2 and (select is None or param in select):
                    self[prefix+param].append(line)
                continue

            # End of parameters
            break

        # Remove trail empty line of the last param
        if self and isinstance(param, basestring) and len(self[prefix+param])>1:
            p = self[prefix+param]
            while p[-1].strip()=='': del p[-1]


    def format(self, select=None, indent=4):
        """Reformat one or several parameters

        :Params:

            - **indent**, optional: Indentation string for all but the first line.

              - String: Used as is.
              - Integer:  ``ident`` times the space char.
              - ``False``: ``''``.
              - Class method: 12 spaces.
              - Other object: 8 spaces.

            - **select**, optional: List of parameter names to format.
              If None, all available parameters are used.


        """
        # Selection of parameters
        if select is None:
            select = self.keys()
            select.sort()
        elif isinstance(select, basestring):
            select = [select]
        # Indentation
        if indent is True or indent is None: indent = 4
        elif not indent: indent = ''
        elif isinstance(indent, int): indent = ' '*indent
        elif not isinstance(indent, basestring):
            indent = (12*' ') if inspect.ismethod(indent) else (8*' ')
        out = []
        for param in select:
            if param not in self: continue
            if len(self[param])==1:
                out.append(self[param][0].strip())
            else:
                out.extend(dedent('\n'.join(self[param])).splitlines())
                if not self[param][-1].strip():
                    out.append('')
        return ('\n'+indent).join(out)


    def asfmtdict(self, indent=8, select=None):
        """Get parameters as a dictionary of formatted descriptions"""
        out = {}
        for param in self:
            if select is not None and param not in select: continue
            out[param] = self.format(param, indent)
        return out


class DocFiller(object):
    '''Scan objects and return a dictionary whose key are their (formatted) name,
    and value are a dictionary of their attributes description

    If the object is a method like ``Map.contour``, the associated key is
    ``Map_contour``. Therefore, classe names are supposed to not have any
    underscore inside.

    Parameter description are indented with 12 spaces for methods
    and 8 spaces for other objects.


    :Example:

        >>> df = DocFiller(MyClass1.my_method, MyClass2, myfync1)


        ::

            @df.docfill
            def myfunc2(para, parb, parc, pard)
                """My description

                :Params:

                    {MyClass1_my_method[para]}
                    {MyClass2[parb]}
                    {MyClass2[parc]}
                    {myfync1[pard]}
    '''
    def __init__(self, *objs, **aobjs):
        self.content = {}
        self.verbose = aobjs.pop('verbose', True)
        if objs: self.scan(*objs, **aobjs)

    def scan(self, *objs, **aobjs):
        aliases = dict([(o, a) for a, o in aobjs.items()])
        for obj in objs+tuple(aliases.keys()):
            key =  obj.__name__
            if inspect.ismethod(obj):
                key = obj.im_class.__name__+"_"+key
            prefix = aliases[obj] if obj in aliases else ''
            self.content[key] = Docstring2Params(obj, prefix=prefix)

    def formatted(self, indent):
        """Format all parameter descriptions

        :Params:

            - **indent**, optional: Indentation string for all but the first line.

              - String: Used as is.
              - Integer:  ``ident`` times the space char.
              - ``False``: ``''``.
              - Class method: 12 spaces.
              - Other object: 8 spaces.

        """
        out = {}
        for key, val in self.content.items():
            out[key] = val.asfmtdict(indent=indent)
        return out

    def docfill(self, obj):
        '''Fill docstring of an object with collected parameter descriptions

        It is typically used to format function docstrings using a decorator.

        .. warning::

            It cannot be applied to methods and classes.

        '''
#        print obj.__doc__
#        print self.formatted(obj)
#        obj.__doc__ = obj.__doc__.format(**self.formatted(obj))
        try:
            obj.__doc__ = obj.__doc__.format(**self.formatted(obj))
        except KeyError, e:
            if self.verbose: print 'Missing key for docfill: '+e.message
#            print 'When formatting doc:'
#            print obj.__doc__
        except Exception, e:
            if self.verbose: print 'Docfill error:', e.message
        return obj

from vacumm import docfiller_verbose
docfiller = DocFiller(verbose=docfiller_verbose)
docfill = docfiller.docfill

if __name__=='__main__':

    df = DocFiller(Docstring2Params.format)

    @df.docfill
    def myfunc(indent):
        """Description

        :Params:

            {Docstring2Params_format[indent]}
        """
        pass

    print myfunc.__doc__


