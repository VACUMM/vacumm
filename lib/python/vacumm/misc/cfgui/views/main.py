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

import collections, os, pprint, re, sys

from PyQt4 import QtCore, QtGui

import configobj

from vacumm.misc.config import ConfigManager, _shelp_, pathname

from vacumm.misc.cfgui import QtObject
from vacumm.misc.cfgui.models.sessions import Session
from vacumm.misc.cfgui.resources.ui.main import Ui_MainWindow
from vacumm.misc.cfgui.utils.ui import confirm_dialog, create_widget, get_widget, set_widget, create_treeview_widget


class ConfigWidgetDict(dict):


    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        self.parent = None # parent ConfigWidgetDict


    def __setitem__(self, key, value):
        value.parent = self
        return dict.__setitem__(self, key, value)


class MainWindow(QtObject, Ui_MainWindow, QtGui.QMainWindow):


    def __init__(self, controller, **kwargs):
        self.controller = controller
        QtObject.__init__(self, **kwargs)
        Ui_MainWindow.__init__(self)
        QtGui.QMainWindow.__init__(self)
        self.setupUi(self)

        self.default_title = self.windowTitle()
        
        self.tabs.tabBar().tabButton(self.tabs.indexOf(self.tab_config), QtGui.QTabBar.RightSide).hide()

        self.model_config = QtGui.QStandardItemModel()
        self.treeview_config.setModel(self.model_config)

        self.treeview_config.setSelectionBehavior(QtGui.QAbstractItemView.SelectItems)

        self.sectionFont = QtGui.QFont()
        self.sectionFont.setWeight(QtGui.QFont.Bold)

        self.sectionIcon = QtGui.QIcon()
        self.sectionIcon.addPixmap(self.style().standardPixmap(QtGui.QStyle.SP_DirClosedIcon),
                QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.sectionIcon.addPixmap(self.style().standardPixmap(QtGui.QStyle.SP_DirOpenIcon),
                QtGui.QIcon.Normal, QtGui.QIcon.On)

        self.optionFont = QtGui.QFont()
        self.optionFont.setWeight(QtGui.QFont.Normal)

        self.optionIcon = QtGui.QIcon()
        self.optionIcon.addPixmap(self.style().standardPixmap(QtGui.QStyle.SP_FileIcon))

        class IgnoreWheelEventFilter(QtCore.QObject):
            def eventFilter(self, source, event):
                if event.type() == QtCore.QEvent.Wheel:
                    return True
                return False
        self.ignoreWheelEventFilter = IgnoreWheelEventFilter()

        self.config_name_validator = QtGui.QRegExpValidator(QtCore.QRegExp('[a-zA-Z_][a-zA-Z_0-9]*'))

        QtGui.QShortcut(QtGui.QKeySequence('ctrl+i'), self, self.log_info)


    def set_application_logo(self, path):
        logo = QtGui.QLabel()
        logo.setPixmap(QtGui.QPixmap(path))
        logo.setAlignment(QtCore.Qt.AlignTop)
        layout = QtGui.QVBoxLayout(self.frame_left)
        layout.addWidget(logo)


    def log_info(self):
        self.logger.info(
            '\n\nconfig widgets:\n\n%s\n\nconfiguration:\n\n%s\n\n',
            pprint.pformat(self.config_widget),
            pprint.pformat(self.get_configuration(self.controller.session.session_config_manager).dict())
        )


    def closeEvent(self, event):
        self.on_menu_file_quit()
        event.ignore()


    def on_menu_file_new(self):
        self.controller.new_configuration()


    def on_menu_file_open(self):
        path = QtGui.QFileDialog.getOpenFileName(filter=Session.configuration_file_filter)
        if not path:
            return
        self.controller.load_configuration(unicode(path))


    def on_menu_file_save(self):
        path = self.controller.session.configuration_file
        if not path:
            path = QtGui.QFileDialog.getSaveFileName(filter=Session.configuration_file_filter)
            if not path:
                return
        self.controller.save_configuration(unicode(path))


    def on_menu_file_save_as(self):
        path = QtGui.QFileDialog.getSaveFileName(filter=Session.configuration_file_filter)
        if not path:
            return
        self.controller.save_configuration(unicode(path))


    def on_menu_file_preferences(self):
        self.controller.application.preferences_controller.show_preferences_dialog()


    def on_menu_file_sessions(self):
        self.controller.application.sessions_controller.show_sessions_dialog()


    def on_menu_file_quit(self):
        if confirm_dialog('Quit ?'):
            self.controller.quit()


    def on_menu_edit_reload(self):
        if confirm_dialog('Reload configuration ?'):
            self.controller.load_configuration()


    def on_menu_edit_load_default(self):
        if confirm_dialog('Load default configuration ?'):
            self.controller.load_default_configuration()


    def on_menu_help_user_guide(self):
        print 'on_menu_help_user_guide'


    def on_menu_help_about(self):
        print 'on_menu_help_about'


    def on_menu_edit_collapse_all(self):
        self.treeview_config.collapseAll()


    def on_menu_edit_expand_all(self):
        self.treeview_config.expandAll()


    def update_session_status(self):
        title = self.default_title
        if self.controller.session.name:
            title += ' - ' + self.controller.session.name
        if self.controller.session.configuration_file:
            title += ' - ' + os.path.abspath(self.controller.session.configuration_file)
        self.setWindowTitle(title)


    def get_tree_label(self, label):
        return label
        #return re.sub(r'[\s_-]+', ' ', label).capitalize()


    def generate_unique_name(self, name, existing, allow=None):
        i = 1
        unique_name = name
        while unique_name in existing or unique_name == allow:
            unique_name = '%s_%d'%(name, i)
            i += 1
        return unique_name


    def autosize_treeview_config(self):
        #self.treeview_config.setColumnWidth(0, 400)
        for column in range(self.model_config.columnCount(QtCore.QModelIndex())):
            self.treeview_config.resizeColumnToContents(column)


    def set_specification(self, cfgman):
        self.verbose('set specification')

        self.model_config.clear()
        self.model_config.setHorizontalHeaderLabels(('Key', 'Value'))

        self.treeview_config.header().setMovable(False)
        #self.treeview_config.header().setResizeMode(0, QtGui.QHeaderView.ResizeToContents)
        #self.treeview_config.header().setResizeMode(1, QtGui.QHeaderView.ResizeToContents)

        cfgspec = cfgman._configspec
        cfgdef = cfgman.defaults()

        section_name = ''
        if self.controller.session.name:
            section_name += self.controller.session.name
        if cfgman._configspecfile:
            section_name += ' - ' + os.path.basename(cfgman._configspecfile)

        self.config_widget = ConfigWidgetDict()

        self.add_section(self.model_config, section_name, cfgman, cfgspec, cfgdef, self.config_widget)

        self.treeview_config.expandAll()
        self.autosize_treeview_config()


    def add_section(self, item_parent, section_name, cfgman, cfgspec, cfgdef, config_widget, can_rename=False, can_remove=False):

        self.verbose('section %s', pathname(cfgspec, section_name))

        item_key = QtGui.QStandardItem()
        item_key.setEditable(False)
        item_key.setIcon(self.sectionIcon)
        #item_key.setFont(self.sectionFont)

        item_value = QtGui.QStandardItem()
        item_value.setEditable(False)

        item_parent.appendRow((item_key, item_value))

        config_widget.item_parent = item_parent
        config_widget.item_key = item_key
        config_widget.item_value = item_value

        frame = QtGui.QFrame()
        layout = QtGui.QHBoxLayout(frame)

        bound_dict_line_to_name = {}

        if can_rename:

            line = QtGui.QLineEdit()
            line.setText(section_name)
            line.setFont(self.sectionFont)
            line.setValidator(self.config_name_validator)

            layout.addWidget(line)

            bound_dict_line_to_name[line] = section_name
            def editingFinished():

                new_name = unicode(line.text())
                previous_name = bound_dict_line_to_name[line]

                # XXX
                current_names = config_widget.parent.keys()
                unique_name = self.generate_unique_name(new_name, config_widget.parent.keys(), previous_name)
                # XXX

                self.debug('before rename section: previous name %r, new name %r, current names %r', previous_name, new_name, current_names)

                if new_name == previous_name:
                    self.debug('section name did not change')
                    return
                if new_name in current_names:
                    self.debug('section name %r already used', new_name)

                self.verbose('rename section %r to %r', previous_name, unique_name)
                line.setText(unique_name)

                config_widget.parent[unique_name] = config_widget.parent.pop(previous_name)

                bound_dict_line_to_name[line] = unique_name
                self.debug('after rename section: config_widget %s bound_dict_line_to_name %s', config_widget.keys(), bound_dict_line_to_name.values())

            line.editingFinished.connect(editingFinished)

            def resizeLine():
                text, font_metrics = line.text(), line.fontMetrics()
                line.setMinimumWidth(max(font_metrics.width(' '*60), font_metrics.width(text)))
                self.autosize_treeview_config()
            resizeLine()
            line.textEdited.connect(resizeLine)

        else:

            label = QtGui.QLabel(self.get_tree_label(section_name))
            label.setFont(self.sectionFont)

            layout.addWidget(label)

        if can_remove:

            removeButton = QtGui.QPushButton('Remove', None)
            removeButton.setIcon(QtGui.QIcon.fromTheme('list-remove'))
            removeButton.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum))

            layout.addWidget(removeButton)

            def removeClicked():
                self.debug('before remove section: config_widget %s bound_dict_line_to_name %s', config_widget.parent.keys(), bound_dict_line_to_name.values())
                item_key.parent().removeRow(item_key.row())
                name = unicode(line.text())
                bound_dict_line_to_name.pop(line)
                config_widget.parent.pop(name)
                self.verbose('remove many section %s', name)
                self.debug('after remove section: config_widget %s bound_dict_line_to_name %s', config_widget.parent.keys(), bound_dict_line_to_name.values())
            removeButton.clicked.connect(removeClicked)

        self.treeview_config.setIndexWidget(item_key.index(), frame)

        # Build options
        for option_name in cfgspec.scalars:
            try:
                self.verbose('option %s', pathname(cfgspec, option_name))

                # Retrieve this option specification (type, default, doc and other properties, ...)
                optspec = cfgman.getspec(cfgspec, option_name)
                self.debug('specification:\n%s', optspec)

                if option_name in ('__many__', '___many___'):

                    self.add_many_option(item_key, option_name, optspec, config_widget)

                else:

                    self.add_option(item_key, option_name, optspec, config_widget)

            except:
                self.exception('Error creating option %s', pathname(cfgspec, option_name))

        # Build subsections
        for section_name in cfgspec.sections:

            if section_name in ('__many__', '___many___'):

                self.add_many_section(item_key, cfgman, cfgspec[section_name], cfgdef.get(section_name, configobj.ConfigObj()), config_widget)

            else:

                config_widget[section_name] = ConfigWidgetDict()

                self.add_section(item_key, section_name, cfgman, cfgspec[section_name], cfgdef.get(section_name, configobj.ConfigObj()), config_widget[section_name])


    def add_many_section(self, item_parent, cfgman, cfgspec, cfgdef, config_widget):

        item_key = QtGui.QStandardItem()
        item_key.setEditable(False)

        item_value = QtGui.QStandardItem()
        item_value.setEditable(False)

        item_parent.appendRow((item_key, item_value))

        frame = QtGui.QFrame()
        layout = QtGui.QHBoxLayout(frame)

        addButton = QtGui.QPushButton('Add section', None)
        addButton.setIcon(QtGui.QIcon.fromTheme('list-add'))
        addButton.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum))
        layout.addWidget(addButton)

        self.treeview_config.setIndexWidget(item_key.index(), frame)

        def on_add_section():

            manysection_name = self.generate_unique_name(unicode('new_section'), config_widget.keys())
            self.debug('add many section %s', manysection_name)

            config_widget[manysection_name] = ConfigWidgetDict()

            self.add_section(item_parent, manysection_name, cfgman, cfgspec, cfgdef, config_widget[manysection_name], can_rename=True, can_remove=True)

        addButton.clicked.connect(on_add_section)


    def add_option(self, item_parent, option_name, optspec, config_widget, can_rename=False, can_remove=False):

        if not optspec.funcname:
            self.warning('option %s: unsupported option specification:\n%s', option_name, pprint.pformat(optspec))

        help = ''
        try:
            help = _shelp_(optspec, option_name, format='%(shelp)s\nDefault: %(default)r', undoc='')
        except Exception, e:
            self.debug('Failed to get help for option %r: %s %s', option_name, e.__class__.__name__, e)

        item_key = QtGui.QStandardItem()
        item_key.setEditable(False)
        item_key.setIcon(self.optionIcon)
        item_key.setFont(self.optionFont)
        item_key.setToolTip(help)

        item_value = QtGui.QStandardItem()
        item_value.setEditable(False)
        item_value.setToolTip(help)

        item_parent.appendRow((item_key, item_value))

        frame = QtGui.QFrame()
        layout = QtGui.QHBoxLayout(frame)

        def bound_enable(state):
            self.debug('option %s enabled: %s', option_name, state)
            if config_widget[option_name].frame:
                config_widget[option_name].frame.setEnabled(state)
            config_widget[option_name].setEnabled(state)
            config_widget[option_name].checkbox.setCheckState(state)
            config_widget[option_name].checkbox.setTristate(False)

        checkbox = None
        if not can_rename and not can_remove:
            checkbox = QtGui.QCheckBox()
            layout.addWidget(checkbox)
            checkbox.setCheckState(True)
            checkbox.setTristate(False)
            checkbox.setToolTip('Enable or disable this option. When disabled, the option will not appear in the output configuration file unless "Save defaults" is enabled in preferences.')
            def stateChanged():
                bound_enable(not not checkbox.checkState())
            checkbox.stateChanged.connect(stateChanged)

        bound_dict_line_to_name = {}

        if can_rename:

            line = QtGui.QLineEdit()
            line.setText(option_name)
            line.setValidator(self.config_name_validator)

            bound_dict_line_to_name[line] = option_name
            def editingFinished():
                new_name = unicode(line.text())
                previous_name = bound_dict_line_to_name[line]

                current_names = config_widget.keys()
                unique_name = self.generate_unique_name(new_name, current_names, previous_name)

                self.debug('before rename option: previous name %r, new name %r, current names %r', previous_name, new_name, current_names)

                if new_name == previous_name:
                    self.debug('option name did not change')
                    return
                if new_name in current_names:
                    self.debug('option name %r already used', new_name)

                self.verbose('rename option %r to %r', previous_name, unique_name)
                line.setText(unique_name)

                config_widget[unique_name] = config_widget.pop(previous_name)

                bound_dict_line_to_name[line] = unique_name
                self.debug('after rename option: current names %s bound_dict_line_to_name %s', config_widget.keys(), bound_dict_line_to_name.values())
            line.editingFinished.connect(editingFinished)

            def resizeLine():
                text, font_metrics = line.text(), line.fontMetrics()
                line.setMinimumWidth(max(font_metrics.width(' '*60), font_metrics.width(text)))
                self.autosize_treeview_config()
            resizeLine()
            line.textEdited.connect(resizeLine)

            layout.addWidget(line)

        else:

            label = QtGui.QLabel(self.get_tree_label(option_name))

            layout.addWidget(label)

        if can_remove:

            removeButton = QtGui.QPushButton('Remove', None)
            removeButton.setIcon(QtGui.QIcon.fromTheme('list-remove'))
            removeButton.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum))

            def removeClicked():
                self.debug('before rename option: config_widget %s bound_dict_line_to_name %s', config_widget.keys(), bound_dict_line_to_name.values())
                item_key.parent().removeRow(item_key.row())
                name = unicode(line.text())
                config_widget.pop(name)
                bound_dict_line_to_name.pop(line)
                self.verbose('remove many option %s', name)
                self.debug('after rename option: config_widget %s bound_dict_line_to_name %s', config_widget.keys(), bound_dict_line_to_name.values())
            removeButton.clicked.connect(removeClicked)

            layout.addWidget(removeButton)

        self.treeview_config.setIndexWidget(item_key.index(), frame)

        widget = create_treeview_widget(optspec, self.treeview_config, item_value)

        widget.installEventFilter(self.ignoreWheelEventFilter)

        widget.item_parent = item_parent
        widget.item_key = item_key
        widget.item_value = item_value
        widget.checkbox = checkbox
        widget.enable = bound_enable

        config_widget[option_name] = widget

        self.autosize_treeview_config()


    def add_many_option(self, item_parent, option_name, optspec, config_widget):

        help = ''
        try:
            help = _shelp_(optspec, option_name, format='%(shelp)s\nDefault: %(default)r', undoc='')
        except Exception, e:
            self.debug('Failed to get help for option %r: %s %s', option_name, e.__class__.__name__, e)

        item_key = QtGui.QStandardItem()
        item_key.setEditable(False)
        item_key.setToolTip(help)

        item_value = QtGui.QStandardItem()
        item_value.setEditable(False)
        item_value.setToolTip(help)

        item_parent.appendRow((item_key, item_value))

        frame = QtGui.QFrame()
        layout = QtGui.QHBoxLayout(frame)

        addButton = QtGui.QPushButton('Add option', None)
        addButton.setIcon(QtGui.QIcon.fromTheme('list-add'))
        addButton.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum))
        layout.addWidget(addButton)

        self.treeview_config.setIndexWidget(item_key.index(), frame)

        def on_add_option():

            manyoption_name = self.generate_unique_name(unicode('new_option'), config_widget.keys())
            self.debug('add many option %s', manyoption_name)

            self.add_option(item_parent, manyoption_name, optspec, config_widget, can_rename=True, can_remove=True)

        addButton.clicked.connect(on_add_option)


    def set_configuration(self, cfgman, config, config_raw):

        self.verbose('set configuration')

        self.debug('config:\n%s', pprint.pformat(config.dict()))

        self.debug('begin ui widgets:\n%s', pprint.pformat(self.config_widget))

        def set_section(cfgspec, config, config_raw, config_widget):

            section_name = cfgspec.name
            self.verbose('section %s', pathname(cfgspec, section_name))

            # Remove any many option/section
            self.logger.debug('clean many options/sections')
            for name, dict_or_widget in config_widget.items():
                if name not in cfgspec.scalars and name not in cfgspec.sections:
                    self.logger.debug('remove many %s %r', 'section' if isinstance(dict_or_widget, dict) else 'option', name)
                    dict_or_widget.item_key.parent().removeRow(dict_or_widget.item_key.row())
                    config_widget.pop(name)

            # Process regular options
            for option_name in cfgspec.scalars:
                try:

                    if option_name in ('__many__', '___many___'):
                        continue

                    self.verbose('set option %s', pathname(cfgspec, option_name))

                    if option_name not in config_widget:
                        self.error('option %s ignored, not found in view', pathname(cfgspec, option_name))
                        continue

                    enabled = config_raw is not None and option_name in config_raw
                    
                    config_widget[option_name].enable(enabled)

                    if option_name not in config:
                        self.verbose('option %s ignored, not found in configuration', pathname(cfgspec, option_name))
                        continue

                    optspec = cfgman.getspec(cfgspec, option_name)
                    optionvalue = config.get(option_name, optspec.default)
                    optiontype = optionvalue.__class__.__name__

                    self.debug('option %s = %r (%s)', pathname(cfgspec, option_name), optionvalue, optiontype)
                    set_widget(optspec, config_widget[option_name], optionvalue)

                except:
                    self.exception('Error setting option %s', pathname(cfgspec, option_name))

            # Process regular sections
            for section_name in cfgspec.sections:

                if section_name in ('__many__', '___many___'):
                    continue

                if section_name not in config_widget:
                    self.error('section %s ignored, not found in view', pathname(cfgspec, section_name))
                    continue

                if section_name not in config:
                    self.verbose('section %s ignored, not found in configuration', pathname(cfgspec, section_name))
                    continue

                set_section(
                    cfgspec[section_name],
                    config[section_name],
                    None if config_raw is None else config_raw.get(section_name, None),
                    config_widget[section_name]
                )

            # Process many options
            for option_name in config.scalars:

                # If present in cfgspec.scalars, it means this is not a many option
                if option_name in cfgspec.scalars:
                    continue
                #if option_name in ('__many__', '___many___'):
                    #continue

                self.verbose('add many option %s', pathname(config, option_name))

                for manyname in ('__many__', '___many___'):
                    if not isinstance(cfgspec.get(manyname, None), (dict, type(None))):
                        optspec = cfgman.getspec(cfgspec, manyname)
                        break
                else:
                    self.error('many option %s ignored, specification not found', pathname(config, option_name))
                    continue

                optionvalue = config.get(option_name, optspec.default)
                optiontype = optionvalue.__class__.__name__

                self.debug('option %s = %r (%s)', pathname(config, option_name), optionvalue, optiontype)
                self.add_option(config_widget.item_key, option_name, optspec, config_widget, can_rename=True, can_remove=True)
                set_widget(optspec, config_widget[option_name], config[option_name])

            # Process many sections
            for section_name in config.sections:

                # If present in cfgspec.sections, it means this is not a many section
                if section_name in cfgspec.sections:
                    continue
                # When using ConfigObj(unrepr=False), which is the default, an option with a dict value will
                # causes this option to appear in config.sections, so we check it is not present in cfgspec.scalars
                if section_name in cfgspec.scalars:
                    continue

                self.verbose('add many section %s', pathname(config, section_name))

                for manyname in ('__many__', '___many___'):
                    if isinstance(cfgspec.get(manyname, None), (configobj.ConfigObj, configobj.Section)):
                        secspec = cfgspec[manyname]
                        break
                else:
                    self.error('many section %s ignored, specification not found', pathname(config, section_name))
                    continue

                config_widget[section_name] = ConfigWidgetDict()

                self.add_section(
                    config_widget.item_key, section_name, cfgman, secspec, 
                    config.get(section_name, configobj.ConfigObj()), config_widget[section_name],
                    can_rename=True, can_remove=True)

                set_section(
                    secspec,
                    config[section_name],
                    None if config_raw is None else config_raw.get(section_name, None),
                    config_widget[section_name]
                )

        cfgspec = cfgman._configspec

        set_section(cfgspec, config, config_raw, self.config_widget)

        self.debug('end ui widgets:\n%s', pprint.pformat(self.config_widget))


    def get_configuration(self, cfgman):

        self.verbose('get configuration')

        self.debug('ui widgets:\n%s', pprint.pformat(self.config_widget))

        cfgspec = cfgman._configspec

        main_config = configobj.ConfigObj()

        def get_section(cfgspec, config, config_widget):

            section_name = cfgspec.name
            self.verbose('section %s', pathname(cfgspec))

            for name, dict_or_widget in config_widget.items():

                if not isinstance(dict_or_widget, dict):

                    option_name = name
                    self.verbose('option %s', pathname(cfgspec, option_name))

                    if dict_or_widget.checkbox and not dict_or_widget.checkbox.checkState():
                        self.verbose('option %s ignored, state is disabled', pathname(cfgspec, option_name))
                        continue

                    if option_name not in cfgspec:
                        for manyname in ('__many__', '___many___'):
                            if not isinstance(cfgspec.get(manyname, None), (dict, type(None))):
                                optspec = cfgman.getspec(cfgspec, manyname)
                                break
                        else:
                            self.error('many option %s ignored, specification not found', pathname(cfgspec, option_name))
                            continue
                    else:
                        optspec = cfgman.getspec(cfgspec, option_name)

                    optionvalue = get_widget(optspec, dict_or_widget)
                    config[option_name] = optionvalue

                    optiontype = optionvalue.__class__.__name__
                    self.debug('option %s = %r (%s)', pathname(cfgspec, option_name), optionvalue, optiontype)

                else:
                    section_name = name

                    if section_name not in cfgspec:
                        for manyname in ('__many__', '___many___'):
                            if isinstance(cfgspec.get(manyname, None), (configobj.ConfigObj, configobj.Section)):
                                secspec = cfgspec[manyname]
                                break
                        else:
                            self.error('many section %s ignored, specification not found', pathname(cfgspec, section_name))
                            continue
                    else:
                        secspec = cfgspec[section_name]

                    config[section_name] = configobj.Section(config, config.depth+1, main_config, indict={}, name=section_name)

                    get_section(secspec, config[section_name], dict_or_widget)

                    if not len(config[section_name]):
                        self.verbose('remove empty section %s', pathname(cfgspec, section_name))
                        config.pop(section_name)

        get_section(cfgspec, main_config, self.config_widget)

        self.debug('config:\n%s', pprint.pformat(main_config.dict()))

        return main_config

