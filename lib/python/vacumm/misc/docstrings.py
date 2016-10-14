# -*- coding: utf8 -*-
"""Docstring reformatting utilities for VACUMM internal usage

.. warning::

    It requires python 2.7 and :mod:`sphinx.ext.napoleon` to work.
"""
# Copyright or © or Copr. Actimar/IFREMER (2010-2016)
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
from textwrap import dedent
import re, inspect

from sphinx.ext.napoleon import GoogleDocstring, NumpyDocstring
from sphinx.ext.napoleon import Config
from sphinx.util.docstrings import prepare_docstring

def get_indent_int(line):
    for i, s in enumerate(line):
        if not s.isspace():
            return i
    return len(line)
def get_indent_str(line):
    return ' '*get_indent_int(line)

class FieldsGetter(NumpyDocstring):
    """A crude derivation of :class:`NumpyDocstring` to intercept field declarations"""
    vacumm_fields = []
    def _parse_parameters_section(self, section):
        fields = self.vacumm_fields = self._consume_fields()
        if self._config.napoleon_use_param:
            lines = []
            for _name, _type, _desc in fields:
                field = ':param %s: ' % _name
                lines.extend(self._format_block(field, _desc))
                if _type:
                    lines.append(':type %s: %s' % (_name, _type))
            return lines + ['']
        else:
            return self._format_fields('Parameters', fields)

class Docstring2Params(dict):
    """Inspect a docstring to retreive and format field declarations"""
    def __init__(self, obj, select=None, prefix=None):

#        fg = FieldsGetter(prepare_docstring(getattr(obj, '__doc__', '')))
        doc = inspect.getdoc(obj) or ''
        fg = FieldsGetter(doc)
        prefix = prefix or ''
        for _name, _type, _desc in fg.vacumm_fields:
            self[prefix+_name] = (_type, _desc)


    def format(self, select=None, indent=1):
        """Reformat one or several fields

        Parameters
        ----------
        indent: optional, string, int, object
            Indentation string for all but the first line.
            A single indentation is four spaces.

            - String: Used as is.
            - Integer:  ``indent`` times the 4 spaces.
            - ``False``: ``''``.
            - Other object: detected from indentation at object declaration, which
              inspected from source code.

        select: optional, strings
            List of parameter names to format.
            If None, all available parameters are used.

        Return
        ------
        string
            Formatted declarations of fields
        """
        # Selection of parameters
        if select is None:
            select = self.keys()
            select.sort()
        elif isinstance(select, basestring):
            select = [select]

        # Indentation
        indent_unit = 4*' '
        if indent is True or indent is None:
            indent = 1
        elif not indent:
            indent = ''
        if not isinstance(indent, (basestring, int)): # object
            indent = get_indent_str(inspect.getsource(indent))
        if isinstance(indent, int):
            indent = indent_unit * indent
        indent = indent + indent_unit # obj indentation + 1 unit

        # Loop on params
        out = []
        for param in select:
            if param not in self: continue
            _type, desc = self[param]
            ss = param
            if _type:
                ss = ss + ': ' + _type
            ss = ss + '\n'
            if desc:
                ss = ss + indent + indent_unit + ('\n'+indent+indent_unit).join(desc)
            out.append(ss)
        return ('\n'+indent).join(out)


    def asfmtdict(self, indent=1, select=None):
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

                Parameters
                ----------

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
        """Format all parameter declarations

        Parameters
        ----------
        indent: optional: string, int, object
            Indentation string for all but the first line.
            A single indentation as four spaces.

            - String: Used as is.
            - Integer:  ``indent`` times the 4 spaces.
            - ``False``: ``''``.
            - Other object: detected from indentation at object declaration, which
              inspected from source code.

        """
        out = {}
        for key, val in self.content.items():
            out[key] = val.asfmtdict(indent=indent)
        return out

    def docfill(self, obj, indent=None):
        '''Fill docstring of an object with collected parameter descriptions

        It is typically used to format function docstrings using a decorator.

        '''
        indent = obj if indent is None else indent
        try:
            obj.__doc__ = obj.__doc__.format(**self.formatted(indent))
        except KeyError, e:
            if self.verbose: print 'Missing key for docfill: '+e.message
        except Exception, e:
            if self.verbose: print 'Docfill error:', e.message
        return obj

    __call__ = docfill

from vacumm.config import VACUMM_CFG
docfiller = DocFiller(verbose=VACUMM_CFG['vacumm.misc.docstrings']['verbose'])
docfill = docfiller.docfill

if __name__=='__main__':

    df = DocFiller(Docstring2Params.format)
    print df.content

    @df
    def myfunc(arg):
        """Description

        Parameters
        ----------

        {Docstring2Params_format[indent]}
        """
        pass

    print myfunc.__doc__


