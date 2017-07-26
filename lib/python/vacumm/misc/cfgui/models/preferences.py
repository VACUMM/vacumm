#!/usr/bin/env python
# -*- coding: utf-8 -*-

from vacumm.misc.cfgui import XmlConfigObject


class Preferences(XmlConfigObject):
	'''
	Simple structure contaning the user's preferences.
	These preferences are divided in several fields:
		
		* The language
		* The default directory used to open/save data
		* The default specification file used to edit configurations
	'''
	
	xml_textnodes = {
		'save_defaults':{'single':True, 'default':True, 'type':bool},
		'save_comments':{'single':True, 'default':True, 'type':bool},
		'backup_count':{'single':True, 'default':0, 'type':int},
	}
	
	def __init__(self, *args, **kwargs):
		XmlConfigObject.__init__(self, *args, **kwargs)
	
