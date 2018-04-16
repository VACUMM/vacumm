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

import os, sys
from PyQt4 import QtCore, QtGui

from vacumm.misc.file import mkfdirs, rollover

from vacumm.misc.cfgui import Object
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
        path = self.application.get_preferences_file()

        self.debug('Checking preferences configuration file %r', path)
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
        path = self.application.get_preferences_file()
        self.notice('Saving preferences to %r', path)
        self.verbose('Preferences:\n%s', self.preferences.to_xml_str())

        try:
            rollover(path, self.preferences.backup_count)
        except Exception, e:
            self.exception('Rollover of file %r failed: %s', path, e)

        mkfdirs(path)
        self.preferences.to_xml_file(path)

