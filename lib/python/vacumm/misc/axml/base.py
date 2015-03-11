# -*- coding: utf8 -*-
"""XML utilities"""
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


#from elementtree import ElementTree as ET
from xml.etree import ElementTree as ET

import operator,os, copy
import xml.dom.minidom
import re
from string import Template

__all__ = ['special_nodes', 'Base', 'def_setget', 'isloadable']


_msg_badtag = "Cannot load an element with a different tag name '%s' (!= '%s')"


# ################
# Useful functions
# ################

def def_setget(param, description, default='_nodefault_', formatter=None):
    """Helps defining a set_param/get_param couplet of methods
    @usage def_setget(<parameter_name>,<description>,default=<default_value>)
    """
    # Define the big macro for both set and get methods
    macro = '''def set_%(param)s(self, %(param)s):
    """Set %(description)s"""
    self._set_("%(param)s", %(param)s, formatter=formatter)
def get_%(param)s(self%(tdefault_value)s, formatter=formatter):
    """Get %(description)s"""
    return self._get_("%(param)s"%(tdefault_call)s)'''
    # Guess fillers
    tdefault_value = tdefault_call = ''
    if default != '_nodefault_':
        tdefault_value = ", default=" + serialize(default)
        tdefault_call = ", default=default"
    return macro % vars()


def check_unicode(value):
    """If it's a string, be sure that is unicode"""
    if isinstance(value, str):
        value = value.decode('utf8')
    return value

def serialize(value, formatter=None, pure_string=False):
    """Serialize a simple python object"""
    # String must be unicode
    value = check_unicode(value)
    # Using a formatter
    if formatter is not None:
        return formatter(value)
    # Unicode string case
    if isinstance(value, unicode):
        if not pure_string:
            try:
                tmp = '"%s"'%value
                eval(tmp)
                value = tmp
            except:
                value = "'%s'"%value
        return value.strip(' \t\n')
    # Generic case
    return str(value)

def unserialize(text, formatter=None, pure_string=False):
    """Unserialize a simple python object"""
    if text is None: return
    # Using a formatter
    return formatter(text, reverse=True)
    # Guess real value
    text = text.strip(' \t\n')
    if pure_string:
        # Strings are not protected
        try:
            value = eval(text)
        except:
            value = text
    else:
        value = eval(text) #FIXME: eval breaks unicode
    # String case
    value = check_unicode(value)
    return value

def isloadable(toload):
    """Check if toload is can be laoded into an xml node"""
    # Direct loading or file loading
    if toload is None or ET.iselement(toload) or isinstance(toload, file): return True
    if isinstance(toload, (str, unicode)):
        # *ml path file
        if os.path.exists(toload) and toload.endswith('ml'): return True
        # xml like string
        toload = toload.lower().strip(' \t\n')
        if toload.startswith('<') and toload.endswith('>'): return True
    # Not loadable
    return False

# ##########
# Main class
# ##########

class Base(ET._ElementInterface,ET.ElementTree, object):
    mandatory_params = []
    pure_string = False
    def __init__(self, toload=None, tag=None, **attrib):

        # Define name of the node
        if tag is None:
            tag = self.__class__.__name__
        self.tag = tag

        # Scan attributes
        initial_values = {}
        for att, val in attrib.items():
            # Check if real attributes
            set_param = 'set_'+att
            if hasattr(self, set_param) and callable(getattr(self, set_param)):
                # Value is for a method
                initial_values[getattr(self, set_param)] = val
                del attrib[att]
            else:
                # Simple attribute converted to string
                attrib[att] = check_unicode(attrib[att])
                if not isinstance(attrib[att], unicode):
                    attrib[att] = str(attrib[att])

        # First: create an empty element
        ET._ElementInterface.__init__(self, self.tag, attrib)
        self._root = self

        # Second: try to load something into an temporary element
        if toload is None: # Nothing ot load, except default parameters
            for met, val in initial_values.items(): met(val)
            return
        else:# We want unicode
            toload = check_unicode(toload)
        if ET.iselement(toload):
            # Load directly from a element
            print 'xml from element'
            assert self.tag == toload.tag, _msg_badtag % (toload.tag,self.tag)
            tmp_element = copy.copy(toload)

        elif (isinstance(toload, unicode) and os.path.exists(toload)) or isinstance(toload, file):
            # Load from a file
            print 'xml from file'
            tmp_element = ET.parse(toload).getroot()

        elif isinstance(toload, unicode):
            # Load from a string
            print 'xml from string'
            tmp_element = ET.fromstring(toload)

        else:
            raise TypeError, "Cant load nothing else but a element, a file name or descriptor, or a xml string."

        # Third: copy text, children and attributes
        self.text = tmp_element.text
        for att, val in tmp_element.attrib.items():
            val = check_unicode(val)
            if not isinstance(val, unicode):
                tmp_element.attrib[att] = str(val)
        self.attrib.update(tmp_element.attrib)
        for child in tmp_element:
            self.append(child)
        del tmp_element

        # Fourth: load special nodes
        self._load_()

        # Finally load initial values
        for met, val in initial_values.items(): met(val)

    def append(self, element):
        """Append a child to this node"""
        self._children.append(element)
        element._root = self._root
        self._parent = self

    def parent(self):
        """Get the parent of this node"""
        return getattr(_parent, None)

    def _root_(self, root=False):
        """Get the root of this tree"""
        if root is True:
            return self._root
        if root in [None, False]:
            return self
        return root

    # Set and get content
    def _load_(self,basenode=None):
        """Loop through the tree checking if we have special nodes"""
        if basenode is None: basenode = self
        special_node_names = [n.__name__ for n in special_nodes]
        for i,node in enumerate(basenode):
            if not isinstance(node,Base) and node.tag in special_node_names:
                basenode[i] = special_nodes[special_node_names.index(node.tag)](toload=node)
            else:
                self._load_(node)
#       print 'before _strip_', ET.tostring(self)
        self.strip()

    def strip(self, root=None):
        """Strip out extra spaces and newlines from the node tree"""
        resubs = [re.compile('^\s+'), re.compile('\s+$')]
        for snode in self._root_(root).getiterator():
            if snode.text is not None:
#               print 'text avant', snode.text
                for res in resubs:
                    snode.text = res.sub('', snode.text)
#               print 'text apres', snode.text


    def _set_(self, param, value, formatter=None, **attrib):
        """Set the text of a node"""
        # Try to find element and create it if does not exists
        param_node = self.find(param)
        if param_node is None:
            param_node = ET.SubElement(self, param)
        # Formatter or serialization
        func = '%s_%s'%(param, formatter)
        if formatter is not None and hasattr(self, func) and callable(getattr(self, func)):
            text = getattr(self, func)(check_unicode(value))
        else:
            text =  serialize(value, formatter=formatter)
        # We update it
        param_node.text = serialize(value, formatter=formatter, pure_string=self.pure_string)
        param_node.attrib.update(attrib)

    def _get_(self, param, default='raise', root=None, formatter=None):
        """Get a node text as string"""
        try:
            param_node = self._root_(root).find(param)
        except:
            param_node = None
        if param_node is None:
            if default == 'raise':
                raise ValueError('No such tag: '+param)
            else:
                return default
        else:
            # Formatter
            func = '%s_%s'%(param, formatter)
            if formatter is None and hasattr(self, func) and callable(getattr(self, func)):
                return getattr(self, func)(param_node.text, reverse=True)
            # Normal unserialization
            return unserialize(param_node.text, formatter=formatter, pure_string=self.pure_string)


    def set(self, key, value):
        """Set an attribute

        Value is safely converted to string
        """
        if not isinstance(value): value = str(value)
        self.attrib[key] = value

    # Input/ouput
    def __str__(self):
        return self.toxml()

    def prettify(self,indent='\t', newl='\n', encoding="utf-8"):
        """Convert to string with pretty indendation
        @keyparam indent: Indentation char [default is tab: '\\t']
        @keyparam newl: Newline char [default is unix: '\\n']
        @keyparam encoding: encoding defaults to 'utf-8'
        """
        return xml.dom.minidom.parseString(str(self)).toprettyxml(indent=indent,
            newl=newl, encoding=encoding).replace('&quot;', '"')
    toprettyxml = prettify

    def ptf(self, indent='\t', newl='\n'):
        pass

    def toxml(self, encoding="utf-8"):
        """Convert to string
        @keyparam encoding: encoding defaults to 'utf-8'
        """
        return ET.tostring(self, encoding=encoding)

    def write(self, file, encoding="utf-8", pretty=False, **kwargs):
        """Write to file
        @keyparam encoding: encoding [default: 'utf-8']
        """
        if not hasattr(file, "write"):
            file = open(file, "wb")
        if pretty:
            file.write(self.prettify(encoding=encoding, **kwargs))
            file.close()
        else:
            if encoding == "utf-8" or encoding == "us-ascii":
                file.write("<?xml version='1.0' encoding='%s'?>\n" % encoding)
            ET.ElementTree.write(self, file, encoding=encoding)
    save = write

    # Checks

    def isset(self, params=None):
        """Is this node fully set, thus usable
        @keyparam params: if None, check if theses parameters (list of names) are set.
        """
        if params is None: return True
        return bool(check_params(params))
    def checkparams(self, params=None, mode=None):
        """Return the list of params that are not set from a specified list of names.
        @param params: list of parameter names to check.
        @keyparam mode: What to do. If None, just return the needed parameters, else raise an error'
        """
        if self.isempty(): return False
        if params is None: params = self.mandatory_params
        params_needed = []
        if not isinstance(params, (list, tuple)):
            params = [params]
        for param_name in params:
            try:
                self._get_text_(param_name, None)
            except:
                params_needed.append(param_name)
        if mode == 'raise':
            assert params_needed, 'This instance needs the following parameters to be set (probably using .set_<parameter>(<value>): %s' % params_needed
        if params_needed: return params_needed
        return False

    def isempty(self):
        """Is this node empty"""
        return len(self) == 0


try:
    len(special_nodes)
except:
    special_nodes = []
special_nodes.append(Base)

