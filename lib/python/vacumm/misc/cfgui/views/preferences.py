#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyQt4 import QtCore, QtGui

from vacumm.misc.cfgui import QtObject
from vacumm.misc.cfgui.models.preferences import Preferences
from vacumm.misc.cfgui.resources.ui.preferences import Ui_PreferencesDialog

class PreferencesDialog(QtObject, Ui_PreferencesDialog, QtGui.QDialog):
	
	def __init__(self, controller, **kwargs):
		self.controller = controller
		QtObject.__init__(self, **kwargs)
		Ui_PreferencesDialog.__init__(self)
		QtGui.QDialog.__init__(self)
		self.setupUi(self)
	
	def show(self):
		self.set_preferences(self.controller.preferences)
		QtGui.QDialog.show(self)
		self.raise_()
		self.activateWindow()
		# XXX
		# defining tristate in designer has no effect
		# setTristate must be called after show
		self.checkbox_save_comments.setTristate(False)
		self.checkbox_save_defaults.setTristate(False)
		# XXX
	
	def set_preferences(self, preferences):
		self.combo_language.clear()
		self.combo_language.addItems(preferences.languages)
		if preferences.language in preferences.languages:
			language_index = preferences.languages.index(preferences.language)
			self.combo_language.setCurrentIndex(language_index)
		self.checkbox_save_defaults.setCheckState(preferences.save_defaults)
		self.checkbox_save_comments.setCheckState(preferences.save_comments)
		self.spinbox_backup_count.setValue(preferences.backup_count)
	
	def get_preferences(self):
		return Preferences(
			language=unicode(self.combo_language.currentText()),
			save_defaults=self.checkbox_save_defaults.isChecked(),
			save_comments=self.checkbox_save_comments.isChecked(),
			backup_count=self.spinbox_backup_count.value()
		)
	
	def accept(self):
		self.controller.set_preferences(self.get_preferences())
		self.controller.save_preferences()
		QtGui.QDialog.accept(self)
	
	def reject(self):
		QtGui.QDialog.reject(self)
