#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
		self.plugins = {}
	
	def run(self, args):
		
		QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)
		
		# Create the application (view)
		self.qapplication = QtGui.QApplication([])
		
		# Exit, exception and signal handlers
		sys.excepthook = self.excepthook
		self.setup_signals()
		
		# Create controllers
		self.preferences_controller = PreferencesController(self)
		self.sessions_controller = SessionsController(self)
		self.main_controller = MainController(self)
		
		# Start the gui
		self.main_controller.show_main_window()
		
		try:
			if args.session:
				self.main_controller.load_session(args.session)
			elif args.spcfile:
				self.main_controller.load_session(Session(specification_file=args.spcfile, configuration_file=args.cfgfile))
		except Exception, e:
			msg = 'Error loading session: %s'%e
			self.exception(msg)
			error_dialog(msg, detail=traceback.format_exc(e))
		
		for name, plugin in self.plugins.items():
			self.verbose('Enabling plugin %r (%r)', name, plugin)
			plugin.enable()
		
		res = self.qapplication.exec_()
		
		return res
	
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
		parser.add_argument('-s', '--spcfile', help='Specification file to use at startup')
		parser.add_argument('-c', '--cfgfile', help='Configuration file to use at startup')
		parser.add_argument('-S', '--session', help='Name of the session to use at startup')
		return parser
		
	@classinstancemethod
	def main(obj, cls, args=()):
		parser = cls.get_argparser()
		Logger.add_argparser_options(parser=parser)
		args = parser.parse_args(args)
		Logger.apply_class_argparser_options(args)
		
		if obj is None:
			app = cls()
		else:
			app = obj
		
		return 1 if app.run(args) else 0
	
	
	def register_plugin(self, plugin_class):
		self.plugins[plugin_class.__class__.__name__] = plugin_class(self)


class Plugin(Object):
	
	def __init__(self, application, **kwargs):
		Object.__init__(self, **kwargs)
		self.application = application
	
	def enable(self):
		raise NotImplementedError('')
	
	def disable(self):
		raise NotImplementedError('')
