#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
from PyQt4 import QtCore, QtGui

from vacumm.misc.config import remove_defaults, _walker_remove_all_comments_
from vacumm.misc.file import mkfdirs, rollover

from vacumm.misc.cfgui import Object
from vacumm.misc.cfgui.utils.ui import error_dialog
from vacumm.misc.cfgui.models.sessions import Session
from vacumm.misc.cfgui.views.main import MainWindow


import logging
class StatusBarLogHandler(logging.Handler):
	
	def __init__(self, window, message_msecs=15000):
		self.window = window
		self.message_msecs = message_msecs
		logging.Handler.__init__(self)
		self.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
		self.setLevel(logging.NOTICE)
	
	def emit(self, record):
		self.window.statusBar().showMessage(self.format(record), self.message_msecs)


class MainController(Object):
	
	def __init__(self, application, **kwargs):
		Object.__init__(self, **kwargs)
		self.application = application
		self.main_window = MainWindow(self)
		
		self.logger.addHandler(StatusBarLogHandler(self.main_window))
	
	@property
	def session(self):
		return self.application.session
	
	def show_main_window(self):
		# Center window
		#frame_geometry = self.main_window.frameGeometry()
		#screen = QtGui.QApplication.desktop().screenNumber(QtGui.QApplication.desktop().cursor().pos())
		#center_point = QtGui.QApplication.desktop().screenGeometry(screen).center()
		#frame_geometry.moveCenter(center_point)
		#self.main_window.move(frame_geometry.topLeft())
		
		rect = self.application.qapplication.desktop().screenGeometry()
		rect.setCoords(rect.width()*0.1, rect.height()*0.1, rect.width()*0.9, rect.height()*0.9)
		self.main_window.setGeometry(rect)
		
		self.main_window.show()
	
	def quit(self):
		self.application.quit()
	
	def new_configuration(self):
		self.notice('Loading new configuration')
		self.session.configuration_file = ''
		self.load_default_configuration()
		self.main_window.update_session_status()
		
	def load_default_configuration(self, fromspec=False):
		self.session.session_config = self.session.session_config_manager.defaults()
		self.notice('Loaded default:\n%s\n', self.pformat(self.session.session_config.dict()))
		self.main_window.set_specification(self.session.session_config_manager)
	
	def load_configuration(self, filepath=None, ascurrentfile=True):
		if filepath:
			self.notice('Loading configuration from %r', filepath)
		else:
			filepath = self.session.configuration_file
			self.notice('Reloading configuration from %r', filepath)
		if ascurrentfile:
			self.session.configuration_file = filepath
		self.session.session_config = self.session.session_config_manager.load(filepath)
		self.info('Configuration:\n%s', self.pformat(self.session.session_config.dict()))
		self.main_window.set_configuration(self.session.session_config_manager, self.session.session_config)
		if filepath:
			self.notice('Configuration loaded from %s'%filepath)
		self.main_window.update_session_status()
	
	def get_view_configuration(self):
		return self.main_window.get_configuration(self.session.session_config_manager)
		
	def save_configuration(self, filepath):
		
		# TODO:
		# - exclude all existing __many__ options when copying from self.session.session_config
		# - exclude only sections that do not match name from main_window.get_configuration
		
		# Build default config
		config = self.session.session_config_manager.defaults()
		# Update with config (at least to get loaded expert values)
		config.update(self.session.session_config.dict())
		# Update with UI values
		config.update(self.get_view_configuration())
		# Remove default/comments
		if not self.application.preferences_controller.preferences.save_comments:
			self.logger.verbose('Removing comments')
			config.walk(_walker_remove_all_comments_, call_on_sections=True)
		if not self.application.preferences_controller.preferences.save_defaults:
			self.logger.verbose('Removing defaults')
			remove_defaults(config)
		
		self.notice('Saving configuration to %r', filepath)
		try: rollover(filepath, self.application.preferences_controller.preferences.backup_count)
		except Exception, e: self.exception('Rollover of file %r failed: %s', filepath, e)
		with file(filepath, 'w') as f:
			config.write(f)
		
		self.info('Configuration:\n%s', self.pformat(config.dict()))
		self.session.configuration_file = filepath
		self.session.session_config = config
		self.main_window.set_configuration(self.session.session_config_manager, self.session.session_config)
		self.main_window.update_session_status()

