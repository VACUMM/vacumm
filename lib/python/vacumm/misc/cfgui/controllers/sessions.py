#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
from PyQt4 import QtCore, QtGui

from vacumm.misc.file import mkfdirs, rollover

from vacumm.misc.cfgui import Object, sessions_configuration_file
from vacumm.misc.cfgui.models.sessions import Sessions, Session
from vacumm.misc.cfgui.views.sessions import SessionsDialog


class SessionsController(Object):
	
	def __init__(self, application, **kwargs):
		Object.__init__(self, **kwargs)
		self.application = application
		self.load_sessions()
		self.sessions_dialog = SessionsDialog(self)
		
	def show_sessions_dialog(self):
		# Center window
		frame_geometry = self.sessions_dialog.frameGeometry()
		screen = QtGui.QApplication.desktop().screenNumber(QtGui.QApplication.desktop().cursor().pos())
		center_point = QtGui.QApplication.desktop().screenGeometry(screen).center()
		frame_geometry.moveCenter(center_point)
		self.sessions_dialog.move(frame_geometry.topLeft())
		self.sessions_dialog.show()
	
	def get_current_session(self):
		return self.application.main_controller.session
	
	def set_session(self, session, replace=None):
		'''
		name: if present, remove existing session before setting the new one (rename)
		'''
		if replace in self.sessions:
			self.delete_session(replace)
		self.sessions[session.name] = session
	
	def get_session(self, name):
		session = self.sessions[name]
		return session
	
	def delete_session(self, name):
		session = self.sessions.pop(name)
	
	def load_session(self, session):
		self.application.main_controller.load_session(session)
	
	def load_sessions(self):
		self.sessions = None
		path = sessions_configuration_file
		if os.path.isfile(path):
			self.notice('Loading sessions from %r', path)
			try:
				self.sessions = Sessions.from_xml_file(path)
			except Exception, e:
				self.exception('Load sessions failed: %s', e)
		if self.sessions is None:
			self.notice('Using default sessions configuration')
			self.sessions = Sessions()
		self.verbose('Sessions:\n%s', self.sessions.to_xml_str())
	
	def save_sessions(self):
		path = sessions_configuration_file
		self.notice('Saving sessions to %r', path)
		self.verbose('Sessions:\n%s', self.sessions.to_xml_str())
		try:
			rollover(path, self.application.preferences_controller.preferences.backup_count)
		except Exception, e:
			self.exception('Rollover of file %r failed: %s', path, e)
		mkfdirs(path)
		self.sessions.to_xml_file(path)

