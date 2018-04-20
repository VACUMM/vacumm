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

import argparse, os, signal, sys, traceback
from PyQt4 import QtCore, QtGui

from vacumm.misc.bases import Object, classinstancemethod
from vacumm.misc.log import Logger

from vacumm.misc.cfgui.utils.ui import error_dialog
from vacumm.misc.cfgui.controllers.main import MainController
from vacumm.misc.cfgui.controllers.preferences import PreferencesController
from vacumm.misc.cfgui.controllers.sessions import SessionsController
from vacumm.misc.cfgui.models.sessions import Session

class Application(Object):


    def __init__(self, **kwargs):
        Object.__init__(self, **kwargs)
        self.configuration_directory = os.path.join(os.path.expanduser('~'), '.cfgui')
        self.plugins = {}


    def run(self, args):

        if args.configuration_directory:
            self.configuration_directory = args.configuration_directory

        QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)

        # Create the application (view)
        self.qapplication = QtGui.QApplication([])

        # Exit, exception and signal handlers
        sys.excepthook = self.excepthook
        self.setup_signals()

        # Hold the session currently in use
        self.session = None

        # Create controllers
        self.preferences_controller = PreferencesController(self)
        self.sessions_controller = SessionsController(self)
        self.main_controller = MainController(self)

        # Setup plugins
        for name, plugin in self.plugins.items():
            self.verbose('Enabling plugin %s', name)
            plugin.enable()

        # Start the gui
        self.main_controller.show_main_window()

        try:
            if args.session:
                self.load_session(args.session)
            elif args.spcfile:
                self.load_session(Session(specification_file=args.spcfile, configuration_file=args.cfgfile))
            else:
                self.load_session(self.create_session())
        except Exception, e:
            msg = 'Error loading session: %s'%e
            self.exception(msg)
            error_dialog(msg, detail=traceback.format_exc(e))

        res = self.qapplication.exec_()

        return res


    def get_configuration_directory(self):
        return self.configuration_directory


    def get_preferences_file(self):
        return os.path.join(self.get_configuration_directory(), 'preferences.xml')


    def get_sessions_file(self):
        return os.path.join(self.get_configuration_directory(), 'sessions.xml')


    def create_session(self, *args, **kwargs):
        self.info('Create session with args: %s, kwargs: %s', args, kwargs)
        session = Session(*args, **kwargs)
        self.info('Created session:\n%s', session.to_xml_str())

        for name, plugin in self.plugins.items():
            if plugin.enabled:
                plugin.session_created(session)
            else:
                self.debug('Ignore disabled plugin %s', name)

        return session


    def load_session(self, session):
        self.notice('Loading session: %r', session)
        if isinstance(session, basestring):
            session = self.sessions_controller.get_session(session)
        self.session = session

        if self.session.specification_file and not os.path.exists(self.session.specification_file):
            error_dialog('Specification file does not exists: %r'%(self.session.specification_file,))

        if self.session.configuration_file and not os.path.exists(self.session.configuration_file):
            error_dialog('Configuration file does not exists: %r'%(self.session.configuration_file,))

        self.session.initialize()

        self.main_controller.main_window.set_specification(self.session.session_config_manager)

        if self.session.session_config:
            self.main_controller.main_window.set_configuration(
                self.session.session_config_manager,
                self.session.session_config,
                self.session.session_config_raw
            )

        self.main_controller.main_window.update_session_status()


    def register_plugin(self, plugin_class):
        self.plugins[plugin_class.__name__] = plugin_class(self)


    def quit(self):
        self.qapplication.quit()


    def excepthook(self, etype, evalue, tb):
        '''
        Display a detailed error message.
        '''
        #self.error('Unhandled exception %s', ''.join(traceback.format_exception(etype, evalue, tb)))
        msg = 'Unhandled exception %s: %s'%(evalue.__class__.__name__, evalue)
        exc = ''.join(traceback.format_exception(etype, evalue, tb))
        self.error('%s\n%s', msg, exc)
        error_dialog(msg, detail=exc)


    def setup_signals(self):
        # Let python interpreter a chance to catch signals, as exposed in QApplication documentation
        if not hasattr(self, '_signal_timer'):
            self._signal_timer = QtCore.QTimer()
            self._signal_timer.timeout.connect(lambda:None)
            self._signal_timer.start(1500)

        signal.signal(signal.SIGINT, self.signal_handler)
        signal.signal(signal.SIGTERM, self.signal_handler)
        #signal.signal(signal.SIGALRM, self.signal_handler)
        signal.signal(signal.SIGHUP, self.signal_handler)
        signal.signal(signal.SIGUSR1, self.signal_handler)
        signal.signal(signal.SIGUSR2, self.signal_handler)


    def signal_handler(self, signum, frame):
        try:
            signame = ', '.join(map(lambda (k,v): k, filter(lambda (k,v): k.startswith('SIG') and v==signum, signal.__dict__.items())))
            self.info('Signal handled signum: %s (%s)', signum, signame)
            if signum == signal.SIGINT: self.quit()
            elif signum == signal.SIGTERM: self.quit()
            elif signum == signal.SIGHUP: pass
            elif signum == signal.SIGUSR1: pass
            elif signum == signal.SIGUSR2: pass
        except Exception, e:
            self.exception(e)


    @classmethod
    def get_argparser(cls):
        parser = argparse.ArgumentParser(description='Configuration editor')
        cls.populate_argparser(parser)
        return parser


    @classmethod
    def populate_argparser(cls, parser):
        parser.add_argument('-s', '--spcfile', help='Specification file to use at startup')
        parser.add_argument('-c', '--cfgfile', help='Configuration file to use at startup')
        parser.add_argument('-n', '--session', help='Name of the session to use at startup')
        parser.add_argument('--configuration_directory', help='Aplication configuration directory')
        Logger.add_argparser_options(parser=parser)


    @classinstancemethod
    def main(obj, cls, args=()):
        parser = cls.get_argparser()
        args = parser.parse_args(args)
        Logger.apply_class_argparser_options(args)

        if obj is None:
            app = cls()
        else:
            app = obj

        return 1 if app.run(args) else 0


class Plugin(Object):


    def __init__(self, application, **kwargs):
        Object.__init__(self, **kwargs)
        self.application = application
        self.enabled = None


    def enable(self):
        self.info('Enabling %s', self.__class__.__name__)
        self.enabled = True


    def disable(self):
        self.info('Disabling %s', self.__class__.__name__)
        self.enabled = False


    def session_created(self, session):
        pass

