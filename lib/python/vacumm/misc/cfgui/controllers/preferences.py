#!/usr/bin/env python
# -*- coding: utf-8 -*-
print __name__
import os, sys
from PyQt4 import QtCore, QtGui

from vacumm.misc.file import mkfdirs, rollover

from vacumm.misc.cfgui import Object, preferences_configuration_file
from vacumm.misc.cfgui.models.preferences import Preferences
from vacumm.misc.cfgui.views.preferences import PreferencesDialog


class PreferencesController(Object):
	
	def __init__(self, application, **kwargs):
		Object.__init__(self, **kwargs)
		self.application = application
		self.load_preferences()
		self.preferences_dialog = PreferencesDialog(self)
	
	def show_preferences_dialog(self):
		# Center window
		frame_geometry = self.preferences_dialog.frameGeometry()
		screen = QtGui.QApplication.desktop().screenNumber(QtGui.QApplication.desktop().cursor().pos())
		center_point = QtGui.QApplication.desktop().screenGeometry(screen).center()
		frame_geometry.moveCenter(center_point)
		self.preferences_dialog.move(frame_geometry.topLeft())
		self.preferences_dialog.show()
	
	def set_preferences(self, preferences):
		self.preferences = preferences
	
	def load_preferences(self):
		self.preferences = None
		path = preferences_configuration_file
		if os.path.isfile(path):
			self.notice('Loading preferences from %r', path)
			try:
				self.preferences = Preferences.from_xml_file(path)
			except Exception, e:
				self.exception('Load preferences failed: %s', e)
		if self.preferences is None:
			self.notice('Using default preferences configuration')
			self.preferences = Preferences()
		self.verbose('Preferences:\n%s', self.preferences.to_xml_str())
	
	def save_preferences(self):
		path = preferences_configuration_file
		self.notice('Saving preferences to %r', path)
		self.verbose('Preferences:\n%s', self.preferences.to_xml_str())
		try:
			rollover(path, self.preferences.backup_count)
		except Exception, e:
			self.exception('Rollover of file %r failed: %s', path, e)
		mkfdirs(path)
		self.preferences.to_xml_file(path)

