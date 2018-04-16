#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright or © or Copr. Actimar/IFREMER (2010-2017)
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
        #self.combo_language.clear()
        #self.combo_language.addItems(preferences.languages)
        #if preferences.language in preferences.languages:
            #language_index = preferences.languages.index(preferences.language)
            #self.combo_language.setCurrentIndex(language_index)
        self.checkbox_save_defaults.setCheckState(preferences.save_defaults)
        self.checkbox_save_comments.setCheckState(preferences.save_comments)
        self.spinbox_backup_count.setValue(preferences.backup_count)


    def get_preferences(self):
        return Preferences(
            #language=unicode(self.combo_language.currentText()),
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

