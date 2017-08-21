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
        path = self.application.get_sessions_file()

        self.debug('Checking sessions configuration file %r', path)
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
        path = self.application.get_sessions_file()
        self.notice('Saving sessions to %r', path)
        self.verbose('Sessions:\n%s', self.sessions.to_xml_str())

        try:
            rollover(path, self.application.preferences_controller.preferences.backup_count)
        except Exception, e:
            self.exception('Rollover of file %r failed: %s', path, e)

        mkfdirs(path)
        self.sessions.to_xml_file(path)

