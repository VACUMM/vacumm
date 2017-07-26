#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from PyQt4 import QtCore

from vacumm.misc.bases import Class, Object as _Object, classmaker
from vacumm.misc.xml import XmlConfig, XmlConfigDict
from vacumm.misc.log import logger # for use as default logger (e.g. in utils)


config_dir = os.path.join(os.path.expanduser('~'), '.cfgui')
preferences_configuration_file = os.path.join(config_dir, 'preferences.xml')
sessions_configuration_file = os.path.join(config_dir, 'sessions.xml')


class Object(_Object):
	def __init__(self, *args, **kwargs):
		_Object.__init__(self, *args, **kwargs)
		self.logger.set_format('[%(asctime)s %(name)s.%(funcName)s %(levelname)-8s] %(message)s')
		#self.logger.set_format('[%(asctime)s %(name)s %(levelname)-8s] %(message)s')


class XmlConfigObject(Object, XmlConfig):
	
	__metaclass__ = classmaker()
	
	def __init__(self, *args, **kwargs):
		Object.__init__(self, **kwargs)
		XmlConfig.__init__(self, *args, **kwargs)
	
	@classmethod
	def xml_element_name(cls, name):
		return name


class XmlConfigDictObject(Object, XmlConfigDict):
	
	__metaclass__ = classmaker()
	
	def __init__(self, *args, **kwargs):
		Object.__init__(self, **kwargs)
		XmlConfigDict.__init__(self, *args, **kwargs)


class QtObjectMeta(QtCore.pyqtWrapperType, Class):
	
	def __init__(cls, name, bases, dct):
		QtCore.pyqtWrapperType.__init__(cls, name, bases, dct)
		Class.__init__(cls, name, bases, dct)


class QtObject(Object):
	
	__metaclass__ = QtObjectMeta
	

