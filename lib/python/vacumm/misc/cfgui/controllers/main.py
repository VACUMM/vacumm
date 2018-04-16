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

import configobj

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


    def load_configuration(self, filepath=None):
        if filepath:
            self.notice('Loading configuration from %r', filepath)
        else:
            filepath = self.session.configuration_file
            self.notice('Reloading configuration from %r', filepath)

        self.session.configuration_file = filepath
        self.session.initialize(reloadcfg=True)
        self.info('Configuration:\n%s', self.pformat(self.session.session_config.dict()))

        self.main_window.set_configuration(
            self.session.session_config_manager,
            self.session.session_config,
            self.session.session_config_raw
        )

        if filepath:
            self.notice('Configuration loaded from %s'%filepath)

        self.main_window.update_session_status()


    def get_view_configuration(self):
        return self.main_window.get_configuration(self.session.session_config_manager)


    def save_configuration(self, filepath):

        # Create empty config
        config = configobj.ConfigObj(indent_type='    ')

        # Update with default config if specified
        if self.application.preferences_controller.preferences.save_defaults:
            config.update(self.session.session_config_manager.defaults())

        # Update with UI values
        config.update(self.get_view_configuration())

        # Remove option having the same value as the default if specified
        if not self.application.preferences_controller.preferences.save_defaults:
            self.logger.verbose('Removing defaults')
            remove_defaults(config)

        # Remove comments
        if not self.application.preferences_controller.preferences.save_comments:
            self.logger.verbose('Removing comments')
            config.walk(_walker_remove_all_comments_, call_on_sections=True)

        self.notice('Saving configuration to %r', filepath)
        try:
            rollover(filepath, self.application.preferences_controller.preferences.backup_count)
        except Exception, e:
            self.exception('Rollover of file %r failed: %s', filepath, e)
        with file(filepath, 'w') as f:
            config.write(f)

        self.info('Configuration:\n%s', self.pformat(config.dict()))
        self.session.configuration_file = filepath
        self.session.session_config = config
        self.main_window.set_configuration(
            self.session.session_config_manager,
            self.session.session_config,
            self.session.session_config_raw
        )
        self.main_window.update_session_status()

