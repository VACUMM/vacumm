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
from __future__ import absolute_import
if __name__ == '__main__':
    import sys; del sys.path[0]

__author__ = 'Jonathan Wilkins'
__email__ = 'wilkins@actimar.fr'
__doc__ = '''
**This module provides XML serialization/deserialization features.**

It is intented to provide a way to define python classes able to load/save data
from/to XML elements, strings or files. The idea is to make a close but non restrictive
link between XML and Python objects models with the least possible effort.

Below is an example to briefly introduce how to use this module.
This could be simplified because there are shortcuts for default type, values
and other stuff but we will be as explicit as possible to demonstrate how this
module work.

Suppose you have or want the following XML syntax:

.. code-block:: xml

    <Report>
        <MapFigure showgrid="True" output="meteo.png">
            <Title>Temperature and Wind</Title>
            <BoundingBox lonmin="-10" lonmax="10" latmin="40" latmax="60">
            <ColorBar levels="0,5,10,15,20,25,30"/>
            <Data file="/path/to/data/file" variable="temperature"/>
            <Data file="/path/to/data/file" variable="wind" plot="windbarbs"/>
        </Map>
    </Report>

This can be serialized with the following python definitions:

.. code-block:: python

    class BoundingBox(XmlConfig):
        xml_attributes = {
            'lonmin':{type:float, default:-180},
            'lonmax':{type:float, default:180},
            'latmin':{type:float, default:-90},
            'latmax':{type:float, default:90},
        }

    def parseFloatList(floatListString):
        return map(float, floatListString.split(','))

    class ColorBar(XmlConfig):
        xml_attributes = {
            'levels':{type:parseFloatList, default:[0, 10, 20, 30]},
        }

    class Data(XmlConfig):
        xml_attributes = {
            'file':{type:unicode, default:None},
            'variable':{type:unicode, default:None},
            'plot':{type=unicode, default:'fillcolor'},
        }

    class MapFigure(XmlConfig):
        xml_attributes = {
            'showgrid':{type:bool, default:False},
            'temperature':{type:unicode, default:'output.png'},
        }
        xml_childnodes = {
            'boundingBox':{'type':BoundingBox, single:True},
            'colorBar':{'type':BoundingBox, single:True},
            'data':{'type':Data},
        }
        xml_textnodes = {
            'title':{'type':unicode, 'tagName':'Title', single:True},
        }

    class Report(XmlConfig):
        xml_childnodes = {
            'mapFigures':{type:MapFigure},
        }

To load a file 'report.xml' containing the XML code above, you can either do:

.. code-block:: python

    r = Report()
    r.load_xml_file('report.xml')

or

.. code-block:: python

    r = Report.from_xml_file('report.xml')

To save a xml file from the loaded object above you can do:

.. code-block:: python

    r.save_xml_file('report.xml')

.. warning::

    * This module make use of xml.etree.ElementTree and works with python >= 2.6.
    * Some XML functionnalities are not handled, there is no namespace
      support although namespaces can be ignored when loading XML data.
      When specifying an attribute name as 'namespace:name' the namespace is ignored
      when loading but is kept when saving.
    * When loading a list of elements into a dictionnary (e.g. using an attribute value as key),
      the order will be kept if the module collections.OrderedDict is available (python > 2.7)

'''

#__all__ = ['XmlConfig', 'XmlConfigList', 'XmlConfigDict']

import codecs, os, re, sys, traceback, urllib2

try: import xml.etree.cElementTree as ElementTree
except ImportError: import xml.etree.ElementTree as ElementTree


# Debugging stuff
_debug = len(os.environ.get('XML_DEBUG', ''))>0
_debug_classes = os.environ.get('XML_DEBUG_CLASS').split(',') if 'XML_DEBUG_CLASS' in os.environ else None
if _debug:
    import pprint
    def log(o, m, *a):
        if o and _debug_classes and o.__class__.__name__ not in _debug_classes: return
        sys.stderr.write('[XML %s] %s\n'%(
            o.__class__.__name__,
            m%tuple(pprint.pformat(aa, indent=2, depth=None) if isinstance(aa, (list,tuple,dict)) else aa for aa in a)))
else:
    def log(*a, **k): pass
log.__doc__ = '''Debugging log function'''


# Check presence of xml.dom.ext for pretty printing
_xmldomext = None
try:
    from xml.dom.ext.reader import Sax2
    from xml.dom.ext import PrettyPrint
    _xmldomext = True
except: pass


def pretty_indent(elem, level=0, indent = '  '):
    '''
    Indent an ElementTree element avoiding extraneous spaces around node text.
    '''
    i = "\n" + level*indent
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + indent
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            pretty_indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


# Use ordered dicts when possible
Dict = dict
try:
    from collections import OrderedDict
    Dict = OrderedDict
except Exception, e:
    log(sys.modules[__name__], 'collections not available: %s', e)


_true = ['y', 'yes', 'true', 'on', '1']
_false = ['n', 'no', 'false', 'off', '0']
def parse_bool(s):
    '''
    Boolean string parser, accepts:
        - %s as True
        - %s as False
    '''%(', '.join(_true), ', '.join(_false))
    _s, s = s, str(s).lower().strip()
    if s in _true:
        return True
    if s in _false:
        return False
    raise ValueError('Could not parse bool from string %r'%_s)


# A "do nothing" function
_emptyfunc = lambda x: x


def register(obj, attributes=None, elements=None, entries=None):
    '''
    Populate an object with data

    :Params:
        - **obj**: the object to be manipulated
        - **attributes**, optional: dict, used with setattr on obj
        - **elements**, optional: list, used with obj.extend
        - **entries**, optional: dict, used with obj.update

    :Return: obj
    '''
    if attributes:
        for k,v in attributes.iteritems():
            setattr(obj, k, v)
    if elements: obj.extend(elements)
    if entries: obj.update(entries)
    return obj


class _no_default(object):
    '''
    Type for representating a 'no default value' to allow None as default value
    later in this code.
    '''
    pass


def usafe(o):
    '''
    Return a unicode textual representation of the given object
    '''
    if isinstance(o, basestring) and not isinstance(o, unicode):
        return unicode(o, encoding='UTF-8', errors='replace')
    elif not isinstance(o, basestring):
        return unicode('%s'%o)
    return o



def xml_element_name(name):
    '''
    Build an xml element name from the given name, this function returns a
    CamelCase name with any non alphanumeric characters omitted (eg. 'Class_A'
    => 'ClassA', 'class_b' => 'Classb' and class_C => 'ClassC')
    '''
    return re.sub(r'[^a-zA-Z0-1]+', '', ''.join(map(lambda s:s[0].upper()+s[1:], name.split())))


#def xml_attribute_name(name):
#    name = xml_element_name(name)
#    return name[0].lower()+name[1:]


#def create_xml_document():
#    return ElementTree.ElementTree()


def create_xml_element(name, attributes=None, childnodes=None, skipNone=True, **kwargs):
    '''
    Create an xml element with given attributes and children

    :Params:
        - **name**: the element tag name
        - **attributes**, optional: dict of attributes name:value
        - **childnodes**, optional: list/tuple/dict of ElementTree.Element object
        - **skipNone**, optional: skip attributes which value is None

    :Return: a ElementTree.Element object
    '''
    elt = ElementTree.Element(name)
    if attributes is not None:
        for attrName, attrValue in attributes.iteritems():
            if attrValue is None and skipNone: continue
            attrValue = '%s'%(attrValue,)
            elt.set(attrName, attrValue)
    if childnodes is not None:
        if not isinstance(childnodes, (dict, list, tuple)):
            childnodes = (childnodes,)
        if isinstance(childnodes, dict):
            for k,c in childnodes.iteritems():
                elt.append(c)
        else:
            for c in childnodes:
                elt.append(c)
    return elt


def splitns(name):
    '''
    Split a 'namespace:name' string into ('namespace', 'name') or (None, 'name')
    if the name doesnt contains a namespace.
    '''
    n = name.split(':')
    namespace, name = n[:2] if len(n) > 1 else (None, name)
    return namespace, name


def findchildren(elt, name):
    '''
    Find children of an element ignoring namespaces

    :Params:
        - **elt**: element to search in
        - **name**: searched element tag name

    :Return: a generator on ElementTree.Element objects
    '''
    ns, name = splitns(name)
    for e in elt.getchildren():
        if e.tag == name: yield e
        elif e.tag[:1] == '{':
            namespace, tagname = e.tag[1:].rsplit('}', 1)
            #if ns and ns != namespace: continue
            if tagname == name: yield e


def getattribute(elt, name, default=None):
    '''
    Get the attribute value of an element ignoring namespaces

    :Params:
        - **elt**: element to search in
        - **name**: searched attribute name
        - **default**, optional: default value when elt has no such attribute

    :Return: attribute value or default
    '''
    ns, name = splitns(name)
    for an, av in elt.attrib.iteritems():
        if an == name: return av
        elif an[:1] == '{':
            namespace, attname = an[1:].rsplit('}', 1)
            #if ns and ns != namespace: continue
            if attname == name: return av
    return default


def hasattribute(elt, name):
    '''
    Check the presence of an element attribute ignoring namespaces

    :Params:
        - **elt**: element to search in
        - **name**: searched attribute name

    :Return: boolean
    '''
    return getattribute(elt, name, default=None) != None


# Internal name of the flag used to track if an attribute was loaded from xml element
__loaded_from_elt__ = '__loaded_from_elt__'
def load_xml_attributes(elt=None, attributes=None, dst=None, create=False):
    '''
    Deserialize attributes from an element

    :Params:
        - **elt**: element to search in
        - **attributes**: specification for attributes deserialization
            This may be either:
                - **a list of attribute names**
                - **a dict of attributes** {name:spec, ...} with name as a string and spec as:
                    - **a type function**
                    - **a list or tuple** (type function, default value)
                    - **a dict** {'type':function, 'default':value}
        - **dst**: an object that will be populated with the attributes
        - **create**: whether to create attributes in dst with the default value when they are not in elt

    The type must be a callable taking the attribute value as single string argument and returning the converted value
    With no type function, attributes are loaded without modification (as string)

    :Return: dict of deserialized attributes
    '''
    if elt is None: elt = ElementTree.Element('')
    m = Dict()
# load all attributes as raw data ? must be configurable
#    if attributes is None:
#        for a in elt.attrib.keys():
#            m[str(a)] = getattribute(elt, a)a)
#    elif isinstance(attributes, (list,tuple)):
    if isinstance(attributes, (list,tuple)):
        for a in attributes:
            if hasattribute(elt, a): m[str(a)] = getattribute(elt, a)
            elif create: m[str(a)] = None
    elif isinstance(attributes, dict):
        for a,dspec in attributes.iteritems():
            nsa = a
            ns, a = splitns(nsa)
            if isinstance(dspec, (list,tuple)): dtype, ddef = dspec[:]
            elif isinstance(dspec, dict): dtype, ddef = dspec.get('type', None), dspec.get('default', None)
            else: dtype, ddef = dspec, None
            #if dtype is bool: dtype = parse_bool
            log(dst, 'load_xml_attributes processing "%s" with spec %s', a, dspec) #####
            if hasattribute(elt, nsa):
                try:
                    m[str(a)] = (getattr(dtype, 'parse', dtype) or _emptyfunc)(getattribute(elt, nsa))
                    dspec[__loaded_from_elt__] = True
                except Exception, e:
                    log(dst, '*** load_xml_attributes failed loading attribute "%s" of element "%s": %s', a, elt.tag, e) #####
                    if create:
                        log(dst, '*** load_xml_attributes setting default value for "%s": "%s"', a, ddef) #####
                        m[str(a)] = ddef
            elif create and dst is not None and not hasattr(dst, a): m[str(a)] = ddef
    if dst is not None:
        register(dst, m)
    return m


def load_xml_childnodes(elt=None, childnodes=None, dst=None, create=False):
    '''
    Deserialize children from an element

    :Params:
        - **elt**: element to search in
        - **childnodes**: specification for children deserialization
            This must be a dict of **object attribute** {name:spec, ...} with name as a string and spec as:
                - **a type function**
                - **a list or tuple** (type function, default value)
                - **a dict with the entries**:
                    - **type**: a XmlConfig subclass
                    - **default**: default value, expected as instance or list or dict of XMlConfig object(s), depending on the single/key parameters
                    - **create**: whether to create the instance or list when there are no such children in elt
                    - **args**: positionnal arguments for XmlConfig children creation
                    - **kwargs**: named arguments for XmlConfig children creation
                    - **single**: whether we expect a single instance or a list of instance as result in dst (default is False). This may also be an integer specifying the child 0-based index.
                    - **self**: whether we want to register loaded children in dst as dict entries (default is False)
                    - **key**: which children attribute to use as key when loading with self=True
        - **dst**: an object that will be populated with the children
        - **create**: whether to create attributes in dst with the default value when they are not in elt

    The type must be a subclass of XmlConfig which will be used to deserialize the children elements

    :Return: a tuple (attributes, elements, entries) of the loaded data
    '''
    if elt is None: elt = ElementTree.Element(dst.__class__.__name__)
    attrs, elements, entries = Dict(), [], Dict()
    for dname,dspec in childnodes.iteritems():
        dname = str(dname)
        if isinstance(dspec, (list,tuple)):
            (dtype, ddef), dcreate, dargs, dkwargs, dsingle, dself, dkey = (
                dspec, False, [], {}, False, False, None)
        elif isinstance(dspec, dict):
            dtype, ddef, dcreate, dargs, dkwargs, dsingle, dself, dkey = (
                dspec.get('type'), dspec.get('default', _no_default),
                dspec.get('create', not dspec.get('single', False)), dspec.get('args', []), dspec.get('kwargs', {}),
                dspec.get('single', False), dspec.get('self', False), dspec.get('key', None))
        else:
            dtype, ddef, dcreate, dargs, dkwargs, dsingle, dself, dkey = dspec, _no_default, False, [], {}, False, False, None
        if dsingle is True: isingle = 0
        elif dsingle is not None and dsingle is not False: dsingle,isingle = True,int(dsingle)
        else: isingle = 0
        if not isinstance(dargs, (list, tuple)): dargs = [dargs]
        log(dst, 'load_xml_childnodes processing "%s" with spec %s', dname, dspec) #####
        l = [(getattribute(c, dkey) if dkey else None, c) for c in findchildren(elt, dtype.tag_name())]
        if len(l):
            if dsingle:
                attrs[dname] = dtype.from_xml_elt(l[isingle][1], *dargs, **dkwargs)
            else:
                objs = [(k, dtype.from_xml_elt(c, *dargs, **dkwargs)) for k,c in l]
                if dself:
                    if dkey: entries.update(Dict(objs))
                    else: elements.extend([c for k,c in objs])
                else:
                    if dkey: attrs[dname] = Dict(objs)
                    else: attrs[dname] = [c for k,c in objs]
        elif not dself:
            if ddef is not _no_default: attrs[dname] = ddef
            elif create:
                if dcreate:
                    # Create default except when same class (avoid recursion)
                    if dsingle:
                        if issubclass(dtype, type(dst)):
                            print __file__,'*** avoided recursion in creation of attribute "%s" with type "%s" on object with same type'%(dname, dtype)
                        else: attrs[dname] = dtype(*dargs, **dkwargs)
                    else: attrs[dname] = Dict() if dkey else []
                else:
                    if dsingle: attrs[dname] = None
                    else: attrs[dname] = Dict() if dkey else []
    if dst is not None:
        register(dst, attrs, elements, entries)
    log(dst, 'load_xml_childnodes done:\n  attributes: %s\n  elements: %s\n  entries: %s', attrs, elements, entries) #####
    return attrs, elements, entries


def load_xml_textnodes(elt=None, textnodes=None, dst=None, create=False, xml_element_name=xml_element_name):
    '''
    Deserialize text node children from an element

    :Params:
        - **elt**: element to search in
        - **textnodes**: specification for text node children deserialization
            This must be a dict of object attribute {name:spec, ...} with name as a string and spec as:
                - **the text nodes tag name**
                - **a list or tuple** (tag name, default value)
                - **a dict with the entries**:
                    - **tagName**: the text nodes tag name (default is the textnode name formatted by xml_element_name)
                    - **type**: a factory function taking the text node content as argument
                    - **default**: default value, expected as instance or list depending on the single parameter
                    - **args**: positionnal arguments for children creation
                    - **kwargs**: named arguments for children creation
                    - **single**: whether we expect a single instance or a list of instance as result in dst (default is False). This may also be an integer specifying the child 0-based index.
        - **dst**: an object that will be populated with the text node children
        - **create**: whether to create attributes in dst with the default value when they are not in elt

    The type must be a callable taking the text node value as single string argument and returning the converted value.

    With no type function, text nodes are loaded without modification (as string)

    :Return: None
    '''
    if elt is None: elt = ElementTree.Element(dst.__class__.__name__)
    for name,spec in textnodes.iteritems():
        log(dst, 'load_xml_textnodes processing "%s" with spec %s', name, spec) #####
        tagName, factory, default, single = xml_element_name(name), _emptyfunc, _no_default, False
        if isinstance(spec, (list,tuple)):
            tagName, default = spec
        elif isinstance(spec, dict):
            tagName, factory, default, single = (spec.get('tagName', tagName),
                spec.get('type', factory), spec.get('default', default), spec.get('single', single))
        elif spec:
            tagName = spec
        if not factory: factory = _emptyfunc
        if factory is bool: factory = parse_bool
        if single is True: isingle = 0
        elif single is not None and single is not False: single,isingle = True,int(single)
        else: isingle = 0
        childElements = [c for c in findchildren(elt, tagName)]
        if childElements:
            if single:
                setattr(dst, name, factory(childElements[isingle].text if childElements[isingle].text else ''))
            else:
                setattr(dst, name, [factory(c.text if c.text else '') for c in childElements])
        elif default is not _no_default:
            setattr(dst, name, default)
        elif create:
            if single: setattr(dst, name, None)
            else: setattr(dst, name, [])


# Internal name used to pass the xml element from XmlConfig.from_xml_elt to XmlConfig.__init__ through XmlConfigMetaClass.__call__
__load_xml_elt__ = '__load_xml_elt__'
class XmlConfigMetaClass(type):
    '''
    Metaclass used to pass the xml element from XmlConfig.from_xml_elt to
    XmlConfig.__init__. This allow XmlConfig descendants to use the loaded
    data in their __init__ straight after XmlConfig.__init__ have been called

    If your implementation already inherit from a metaclass based class, you may use the code at:
    http://code.activestate.com/recipes/204197-solving-the-metaclass-conflict/
    '''
    def __call__(cls, *args, **kwargs):
        e = kwargs.pop(__load_xml_elt__, None)
        o = cls.__new__(cls, *args, **kwargs)
        setattr(o, __load_xml_elt__, e)
        o.__init__(*args, **kwargs)
        return o


class XmlConfig(object):
    '''
    Base class for XML serializations/deserializations.

    Subclasses can define the following **class attributes** which define the serialization specifications:
        - **xml_tag_name**: The corresponding xml element name (default is the class name formatted by :meth:`xml_element_name`).
        - **xml_attributes**: Attribute specifications, see :func:`load_xml_attributes`
        - **xml_childnodes**: Children elements specifications, see :func:`load_xml_childnodes`
        - **xml_textnodes**: Text children elements specifications, see :func:`load_xml_textnodes`

    **These attributes can also be passed as named arguments to the constructor to extend or overwrite
    those declared in the class definition.**

    **XML attributes, childnodes and textnodes are stored in this instance (object) attributes.**

    '''
    __metaclass__ = XmlConfigMetaClass

    xml_tag_name = ''
    xml_attributes = Dict()
    xml_childnodes = Dict()
    xml_textnodes = Dict()

    def __init__(self, *args, **kwargs):
        '''
        Can receive extra specification in kwargs with the same keys as the class special attributes
        (xml_tag_name, xml_attributes,xml_childnodes and xml_textnodes)

        kwargs can contain the **xml_create_attributes** flag specifying whether to initalize expected
        attributes/childnodes/textnodes with None/empty-list/empty-dict depending on specifications (default is True)

        Other parameters in kwargs are used as initialization values when they are known by (present in) the given in specifications
        (the order of processing is then attributes, childnodes, textnodes), the user is then responsible of the types of these initial values.

        This class also maintain the currently loaded/saved file if any, see set_current_file and get_current_file
        '''
        log(self, '__init__:\n  args: %s\n  kwargs: %s)', args, kwargs) #####
        # Attributes
        self.xml_attributes = self.xml_attributes_specs(self.xml_attributes)
        xml_attributes = kwargs.pop('xml_attributes', None)
        if xml_attributes:
            self.xml_attributes.update(self.xml_attributes_specs(xml_attributes))
        # Child nodes
        self.xml_childnodes = self.xml_childnodes_specs(self.xml_childnodes)
        xml_childnodes = kwargs.pop('xml_childnodes', None)
        if xml_childnodes:
            self.xml_childnodes.update(self.xml_childnodes_specs(xml_childnodes))
        # Text nodes
        self.xml_textnodes = self.xml_textnodes_specs(self.xml_textnodes)
        xml_textnodes = kwargs.pop('xml_textnodes', None)
        if xml_textnodes:
            self.xml_textnodes.update(self.xml_textnodes_specs(xml_textnodes))
        # Force creation of object attributes
        # TODO: create defaults only when they dont already exists or wont be created by load_xml_elt below
        if kwargs.pop('xml_create_attributes', True):
            load_xml_attributes(None, self.xml_attributes, self, True)
            load_xml_childnodes(None, self.xml_childnodes, self, True)
            load_xml_textnodes(None, self.xml_textnodes, self, True, xml_element_name=self.xml_element_name)
        for nsa,s in self.xml_attributes.items():
            ns, a = splitns(nsa)
            if a in kwargs:
                v = kwargs[a]
                t = s.get('type', None)
                if t:
                    try:
                        v = t(v)
                    except:
                        raise TypeError('Unsupported type %r, %s required'%(v.__class__, t))
                setattr(self, a, v)
        for a,s in self.xml_childnodes.items():
            if a in kwargs:
                v = kwargs[a]
                t = s.get('type', None)
                if s.get('single', None):
                    if t and not isinstance(v, t):
                        raise TypeError('Unsupported type %r, %s required'%(v.__class__, t))
                else:
                    if s.get('key', None):
                        if t and not isinstance(v, dict):
                            raise TypeError('Unsupported type %r, %r required'%(v.__class__, dict))
                        for vv in v.values():
                            if t and not isinstance(vv, t):
                                raise TypeError('Unsupported type %r, %r required'%(vv.__class__, t))
                    else:
                        if t and not isinstance(v, (list,tuple)):
                            raise TypeError('Unsupported type %r, %r required'%(v.__class__, (list,tuple)))
                        for i, vv in enumerate(v):
                            if isinstance(vv, dict):
                                v[i] = t(**vv)
                            elif t and not isinstance(vv, t):
                                raise TypeError('Unsupported type %r, %r required'%(vv.__class__, t))
                setattr(self, a, v)
        for a,s in self.xml_textnodes.items():
            if a in kwargs:
                setattr(self, a, kwargs[a])
        #import pdb; pdb.set_trace()
        log(self, '__init__ done:\n  xml_attributes: %s\n  xml_childnodes: %s\n  xml_textnodes: %s', self.xml_attributes, self.xml_childnodes, self.xml_textnodes) #####
        # Check if we are created from factory function from_xml_elt, then load the element here,
        # allowing the implementation to use some loaded data in its __init__ constructor
        e = getattr(self, __load_xml_elt__, None)
        if e is not None:
            self.load_xml_elt(e)
            delattr(self, __load_xml_elt__)

    # TODO: validate modifications using setattr; more actions needed for nested lists/dicts (for children)
    #def __setattr__(self, k, v):

    @classmethod
    def dict(cls, *a, **k):
        '''Wrapper to collections.OrderedDict if available, dict otherwise'''
        return Dict(*a, **k)

    @classmethod
    def xml_attributes_specs(cls, attrs, copy=True): # copy=False has odd side effects depending on usages ... !
        '''Fix attrs as used internally'''
        if isinstance(attrs, dict):
            if copy: attrs = attrs.copy()
        elif isinstance(attrs, (list,tuple)):
            attrs = Dict([(a,{}) if isinstance(a, basestring) else a for a in attrs])
        else:
            raise TypeError('%s: Invalid attributes specifications: %r'%(cls, attrs))
        for a,spec in attrs.iteritems():
            d = spec
            if isinstance(spec, (list, tuple)):
                d = {'type':spec[0], 'default':spec[1]}
            elif not isinstance(spec, dict):
                d = {'type':spec, 'default':None}
            if d.get('type', None) is bool:
                d['type'] = parse_bool
                d['formatter'] = lambda b: str(b).lower()
            attrs[a] = d
        return attrs

    @classmethod
    def xml_childnodes_specs(cls, childnodes, copy=True): # copy=False has odd side effects depending on usages ... !
        '''Fix childnodes as used internally'''
        if isinstance(childnodes, dict):
            if copy: childnodes = childnodes.copy()
        else:
            raise TypeError('%s: Invalid childnodes specifications: %r'%(cls, childnodes))
        for c,spec in childnodes.iteritems():
            d = spec
            if isinstance(spec, (list, tuple)):
                d = {'type':spec[0], 'default':spec[1]}
            elif not isinstance(spec, dict):
                d = {'type':spec}
            childnodes[c] = d
        return childnodes

    @classmethod
    def xml_textnodes_specs(cls, textnodes, copy=True): # copy=False has odd side effects depending on usages ... !
        '''Fix textnodes as used internally'''
        if isinstance(textnodes, dict):
            if copy: textnodes = textnodes.copy()
        else:
            raise TypeError('%s: Invalid textnodes specifications: %r'%(cls, textnodes))
        for t,spec in textnodes.iteritems():
            d = spec
            if isinstance(spec, (list, tuple)):
                d = {'type':spec[0], 'default':spec[1]}
            elif not isinstance(spec, dict):
                d = {'type':spec, 'default':None}
            if d.get('type', None) is bool: d['type'] = parse_bool
            textnodes[t] = d
        return textnodes

    def set_current_file(self, f):
        '''Set the currently used xml file. This method automatically called when loading/saving xml file'''
        self.xml_current_file = f

    def get_current_file(self):
        '''Get the currently used xml file (from last load/save or manual set)'''
        return getattr(self, 'xml_current_file', '')

    @classmethod
    def tag_name(cls):
        '''
        Return the xml element name for this class or instance.
        If cls.xml_tag_name is not set, return class name formatted by :func:`xml_element_name`
        '''
        if cls.xml_tag_name: return cls.xml_tag_name
        return cls.xml_element_name((isinstance(cls, type) and cls or cls.__class__).__name__)

    @classmethod
    def xml_element_name(cls, *a, **k):
        return xml_element_name(*a, **k)

    @classmethod
    def from_xml_doc(cls, d, *a, **k):
        '''Factory method instanciating from xml dom document.'''
        return cls.from_xml_elt(d.getroot() if isinstance(d, ElementTree.ElementTree) else d, *a, **k)

    @classmethod
    def from_xml_file(cls, f, *a, **k):
        '''Factory method instanciating from xml file.'''
        log(cls, 'from_xml_file "%s"', f) #####
        o = cls.from_xml_doc(ElementTree.parse(f), *a, **k)
        o.set_current_file(f)
        return o

    @classmethod
    def from_xml_str(cls, s, *a, **k):
        '''Factory method instanciating from xml string.'''
        return cls.from_xml_doc(ElementTree.fromstring(s), *a, **k)

    @classmethod
    def from_url(cls, url, *a, **k):
        '''Factory method instanciating from xml string returned by url.'''
        return cls.from_xml_str(urllib2.urlopen(url).read(-1), *a, **k)

    @classmethod
    def from_xml_elt(cls, e, *a, **k):
        '''Factory method instanciating from xml  element.'''
        try:
            o = cls(*a, __load_xml_elt__=e, **k)
        except Exception, x:
            xmlmax = 2**10*100
            tb = traceback.format_exc()
            try:
                xml = ElementTree.tostring(e)[:xmlmax]
            except Exception, ee: xml = '[%s]'%ee
            raise x.__class__('Failed to create %s from xml element:\n%s\n\n%s'%(
                cls.__name__, xml, tb))
        return o

    def load_xml_doc(self, d, *a, **k):
        '''Load from an ElementTree object'''
        e = d.getroot() if isinstance(d, ElementTree.ElementTree) else d
        self.load_xml_elt(e, *a, **k)

    def load_xml_file(self, f, *a, **k):
        '''Load from a file'''
        log(self, 'load_xml_file "%s"', f) #####
        self.load_xml_doc(ElementTree.parse(f), *a, **k)
        self.set_current_file(f)

    def load_xml_str(self, s, *a, **k):
        '''Load from a string'''
        self.load_xml_doc(ElementTree.fromstring(s), *a, **k)

    def load_xml_elt(self, e, *a, **k):
        '''Load configurations from an xml element.'''
        if not hasattr(self, __load_xml_elt__):
            self.pre_load_xml(e, *a, **k)
        load_xml_attributes(e, self.xml_attributes, self)
        load_xml_childnodes(e, self.xml_childnodes, self)
        load_xml_textnodes(e, self.xml_textnodes, self, xml_element_name=self.xml_element_name)
        log(self, 'load_xml_elt done:\n  elt: %s\n  args: %s\n  kwargs: %s\n  vars:%s', e, a, k, vars(self)) #####
        if not hasattr(self, __load_xml_elt__):
            self.post_load_xml(e, *a, **k)

    def pre_load_xml(self, *a, **k):
        '''Called before load_xml_elt is executed'''
        pass

    def post_load_xml(self, *a, **k):
        '''Called after load_xml_elt is executed'''
        pass

    def to_xml_elt(self, **kwargs):
        '''Dump object to an xml element.'''
        dumpnonloaded = kwargs.get('dumpnonloaded', True) # can be overriden in attr spec
        dumpdefault = kwargs.get('dumpdefault', True) # can be overriden in attr spec
        # Attributes
        attributes = Dict()
        for a, spec in self.xml_attributes.iteritems():
            nsa = a
            ns, a = splitns(nsa)
            if hasattr(self, a):
                v = getattr(self, a)
                if v is None: continue
                if not spec.get(__loaded_from_elt__, spec.get('dumpnonloaded', dumpnonloaded)): continue
                d = spec.get('default', None)
                if v == d and not spec.get('dumpdefault', dumpdefault): continue
                f = spec.get('formatter', None)
                if f: v = f(v)
                #attributes['{%s}%s'%(ns,a) if ns else a] = usafe(v)
                attributes['%s:%s'%(ns,a) if ns else a] = usafe(v)
        childnodes = []
        # Text nodes
        for name,spec in self.xml_textnodes.iteritems():
            v = getattr(self, name, None)
            if v is None: continue
            if not spec: spec = {}
            tagName = spec.get('tagName', self.xml_element_name(name))
            f = spec.get('formatter', None)
            if spec.get('single', False):
                if f: v = f(v)
                e = ElementTree.Element(tagName)
                e.text = usafe(v)
                childnodes.append(e)
            else:
                if not isinstance(v, (list, tuple)):
                    raise TypeError('List of textnodes values expected for %s.%s (spec: %r), got %r'%(self.tag_name(), name, spec, type(v)))
                for t in v:
                    if f: t = f(t)
                    e = ElementTree.Element(tagName)
                    e.text = usafe(t)
                    childnodes.append(e)
        # Child nodes
        for child,spec in self.xml_childnodes.iteritems():
            tmpchildnodes = []
            log(self, 'to_xml_elt\n  %s', vars(self)) #####
            if spec.get('self', False):
                if spec.get('key', None):
                    for k,o in self.iteritems():
                        e = o.to_xml_elt(**kwargs)
                        if not e.attrib.has_key(spec['key']):
                            e.set(spec['key'], usafe(k))
                        tmpchildnodes.append(e)
                else:
                    tmpchildnodes.extend((o.to_xml_elt(**kwargs) for o in self))
            elif hasattr(self, child):
                if spec.get('single', False):
                    if getattr(self, child) is not None:
                        #print self, child, getattr(self, child)
                        tmpchildnodes.append(getattr(self, child).to_xml_elt(**kwargs))
                else:
                    if spec.get('key', None):
                        for k,o in getattr(self, child).iteritems():
                            e = o.to_xml_elt(**kwargs)
                            if not e.attrib.has_key(spec['key']):
                                e.set(spec['key'], k)
                            tmpchildnodes.append(e)
                    else:
                        #print self, child, spec, getattr(self, child), kwargs
                        tmpchildnodes.extend((o.to_xml_elt(**kwargs) for o in getattr(self, child)))
            # Exclude empty nodes (no attributes and no child)
            childnodes.extend(c for c in tmpchildnodes if len(c.attrib) or len(c.getchildren()))
        return create_xml_element(self.tag_name(), attributes=attributes, childnodes=childnodes, **kwargs)

    def to_xml_doc(self, **kwargs):
        '''Dump object to an ElementTree'''
        return ElementTree.ElementTree(self.to_xml_elt(**kwargs))

    def to_xml_stream(self, stream, pretty=True, indent='  ', encoding='UTF-8', **kwargs):
        '''Dump object to a file object like stream.'''
        close = False
        if isinstance(stream, basestring):
            close = True
            if _xmldomext: stream = file(stream, 'w')
            else: stream = codecs.open(stream, mode='w', encoding=encoding, errors='replace')
        try:
            e = self.to_xml_elt(**kwargs)
            if pretty:
                if _xmldomext:
                    PrettyPrint(Sax2.Reader().fromString(ElementTree.tostring(e)),
                        stream=stream, encoding=encoding, indent=indent, preserveElements=None)
                else:
#                    minidom.parseString(
#                        ElementTree.tostring(e)).writexml(
#                            stream, addindent=indent, newl='\n')
                    pretty_indent(e)
                    stream.write(ElementTree.tostring(e))
            else:
                d = ElementTree.ElementTree(e)
                #d.write(stream, xml_declaration=True, method="xml")
                d.write(stream, encoding=encoding, xml_declaration=True, method="xml")
        finally:
            if close: stream.close()
        return e

    def to_xml_str(self, **kwargs):
        '''Dump object to a xml string.'''
        from StringIO import StringIO
        s = StringIO()
        self.to_xml_stream(s, **kwargs)
        return s.getvalue()

    def to_xml_file(self, f, **kwargs):
        '''Dump object to a xml file.'''
        log(self, 'to_xml_file "%s"', f) #####
        e = self.to_xml_stream(f, **kwargs)
        self.set_current_file(f)
        return e

    def to_obj(self):#def to_obj(self, obj=None):
        def to_obj(obj):
            if isinstance(obj, XmlConfig):
                d = dict()
                for a,s in obj.xml_attributes.iteritems():
                    if hasattr(obj, a):
                        v = getattr(obj, a)
                        f = s.get('formatter', None)
                        if f: v = f(v)
                        d[a] = v
                for a,s in obj.xml_childnodes.iteritems():
                    if s.get('self', None):
                        if s.get('key', None):
                            dd = {}
                            for aa,vv in obj.iteritems():
                                dd[aa] = vv.to_obj()
                            d['__items__'] = dd
                        else:
                            l = []
                            for vv in obj:
                                l.append(vv.to_obj())
                            d['__items__'] = l
                    elif hasattr(obj, a):
                            if s.get('single', None):
                                d[a] = getattr(obj, a).to_obj()
                            else:
                                if s.get('key', None):
                                    d[a] = dict((kk, vv.to_obj()) for kk,vv in getattr(obj, a).iteritems())
                                else:
                                    d[a] = [vv.to_obj() for vv in getattr(obj, a)]
                for a,s in obj.xml_textnodes.iteritems():
                    d[a] = getattr(obj, a)
                    if s:
                        f = s.get('formatter', None)
                        if f: d[a] = f(d[a])
                return d
            if isinstance(obj, (list,tuple)):
                return [to_obj(o) for o in obj]
            if isinstance(obj, dict):
                return dict(((k,to_obj(v)) for k,v in obj.iteritems()))
            return obj
        return to_obj(self)

    def set_obj(self, obj):
        for k,s in self.xml_attributes.items():
            if k in obj:
                v = obj[k]
                if v is None: pass
                elif s.get('type', None): v = s['type'](v)
                setattr(self, k, v)
        # TODO: childnodes ?
        for k,s in self.xml_textnodes.items():
            if k in obj:
                v = obj[k]
                if v is None: pass
                elif s.get('type', None): v = s['type'](v)
                setattr(self, k, v)
        return obj


#    # test method to validate to_obj
#    def to_obj_to_json(self):
#        import pprint
#        obj = self.to_obj()
#        print pprint.pformat(obj, indent=2)
#        import json
#        return json.JSONEncoder(sort_keys=True, indent=4).encode(obj)

# Need update for textnodes (at least..:)
#    def to_json(self):
#        import json
#        class JSONEncoder(json.JSONEncoder):
#            def default(self, obj):
#                if isinstance(obj, XmlConfig):
#                    d = dict(((a, getattr(obj, a)) for a in obj.xml_attributes.iterkeys() if hasattr(obj, a)))
#                    for a,v in obj.xml_childnodes.iteritems():
#                        if v.get('self', None):
#                            if v.get('key', None):
#                                dd = {}
#                                for aa,vv in obj.iteritems():
#                                    dd[aa] = self.default(vv)
#                                d['__items__'] = dd
#                            else:
#                                l = []
#                                for vv in obj:
#                                    l.append(self.default(vv))
#                                d['__items__'] = l
#                        elif hasattr(obj, a):
#                            d[a] = self.default(getattr(obj, a))
#                    return d
#                #return json.JSONEncoder.default(self, obj)
#                return obj
#        return JSONEncoder(sort_keys=True, indent=4).encode(self)


class XmlConfigList(list, XmlConfig):
    '''Convenient class inheriting list and XmlConfig'''

    def __init__(self, *args, **kwargs):
        list.__init__(self, *args)
        XmlConfig.__init__(self, **kwargs)


class XmlConfigDict(Dict, XmlConfig):
    '''Convenient class inheriting Dict and XmlConfig to avoid metaclasses conflict when Dict is collections.OrderedDict'''
    if hasattr(Dict, '__metaclass__'):
        class XmlConfigDictMeta(Dict.__metaclass__, XmlConfig.__metaclass__): pass
        __metaclass__ = XmlConfigDictMeta

    def __init__(self, *args, **kwargs):
        Dict.__init__(self, *args)
        XmlConfig.__init__(self, **kwargs)


def _test():
    print ' TEST XML '.center(80, '=')
    class SKlass(XmlConfig): pass
    class Klass(XmlConfig):
        def __init__(self):
            XmlConfig.__init__(self,
                xml_attributes=self.dict([
                    ('name',(None,None)),
                    ('default',(None,'default')),
                    ('boolean',(bool,False)),
                ]),
                xml_childnodes=self.dict([
                    ('childnodes',{'type':Klass, 'create':True}),
                    ('single_childnode',{'type':SKlass, 'single':True, 'create':True}),
                ]),
                xml_textnodes=self.dict([
                    ('textnode1', {'single':True}),
                    ('textnodes', None),
                ]),
            )
    x = '''
<Klass>
    <Klass name="a"/>
    <Klass name="b" boolean="true"/>
    <SKlass/>
    <Textnode1>textnode1</Textnode1>
    <Textnodes>text1</Textnodes>
    <Textnodes>text2</Textnodes>
</Klass>
    '''
    print ' xml '.center(80, '-')
    print x
    o = Klass.from_xml_str(x)
    print (' object from_xml_str %s '%(id(o))).center(80, '-')
    print '\n'.join('%-20s : %-20s : %s'%(k,type(v),repr(v)) for k,v in o.__dict__.iteritems() if not k.startswith('xml_'))
    print (' object to_xml_str %s '%(id(o))).center(80, '-')
    print o.to_xml_str()
    print (' object to_obj %s '%(id(o))).center(80, '-')
    import pprint
    print pprint.pformat(o.to_obj())
    print (' object from_xml_str(to_xml_str).to_xml_str %s '%(id(o))).center(80, '-')
    print Klass.from_xml_str(o.to_xml_str()).to_xml_str()

if __name__ == '__main__':
    _test()


