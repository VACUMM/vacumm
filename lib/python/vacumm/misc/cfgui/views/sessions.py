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
from vacumm.misc.cfgui.models.sessions import Session
from vacumm.misc.cfgui.resources.ui.sessions import Ui_SessionsDialog
from vacumm.misc.cfgui.utils.ui import question_dialog


class SessionsDialog(QtObject, Ui_SessionsDialog, QtGui.QDialog):


    def __init__(self, controller, **kwargs):
        self.controller = controller
        QtObject.__init__(self, **kwargs)
        Ui_SessionsDialog.__init__(self)
        QtGui.QDialog.__init__(self)
        self.setupUi(self)

        self.combo_name.setAutoCompletion(False)
        self.combo_name.setCompleter(None)

        self.current_session_name = None


    def show(self):
        self.set_session(self.controller.get_current_session())
        QtGui.QDialog.show(self)
        self.raise_()
        self.activateWindow()


    def set_session(self, session):
        self.debug('set_session:\n%s', session.to_xml_str())
        self.current_session_name = session.name
        self.combo_name.setEditText(session.name)
        self.line_specification.setText(session.specification_file or '')
        self.line_configuration.setText(session.configuration_file or '')
        self.update_combo_name()
    
    def get_session(self):
        session = self.controller.application.create_session(
            name=unicode(self.combo_name.currentText()),
            specification_file=unicode(self.line_specification.text()),
            configuration_file=unicode(self.line_configuration.text())
        )
        self.debug('get_session:\n%s', session.to_xml_str())
        return session


    def update_combo_name(self):
        session = self.get_session()
        items = sorted(self.controller.sessions.keys())
        if session.name not in items:
            items.append(unicode(session.name))
        name = unicode(self.combo_name.currentText())
        self.debug('update_combo_name: %(items)r %(name)r', locals())

        self.combo_name.clear()
        self.combo_name.addItems(items)
        #self.combo_name.setEditText(session.name)
        if name in items:
            self.combo_name.setCurrentIndex(items.index(session.name))


    def on_button_new(self):
        self.set_session(self.controller.application.create_session(name='New session'))


    def on_button_save(self):
        session = self.get_session()
        self.controller.set_session(session, replace=self.current_session_name)
        self.controller.save_sessions()
        self.set_session(session)


    def on_button_open(self):
        self.hide()
        self.controller.application.load_session(self.get_session())


    def on_button_delete(self):
        name = unicode(self.combo_name.currentText())
        if name and question_dialog('Delete session "%s"'%name):
            self.controller.delete_session(name)
            self.controller.save_sessions()
            self.update_combo_name()
            self.set_session(self.controller.application.create_session())


    def on_button_close(self):
        #QtGui.QDialog.reject(self)
        self.hide()


    def on_combo_name(self):
        name = unicode(self.combo_name.currentText())
        self.set_session(self.controller.get_session(name))


    def on_button_specification(self):
        path = QtGui.QFileDialog.getOpenFileName(filter=Session.specification_file_filter)
        if not path:
            return
        self.line_specification.setText(path)


    def on_button_configuration(self):
        path = QtGui.QFileDialog.getOpenFileName(filter=Session.configuration_file_filter)
        if not path:
            return
        self.line_configuration.setText(path)

