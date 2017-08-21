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

import datetime, os, pprint, sys, traceback
from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QMessageBox

from vacumm.misc.cfgui import logger


def dialog(text=None, title=None, info=None, detail=None, icon=None, buttons=None, default_button=None, escape_button=None):
    '''
        text: string
        icon: QMessageBox.Information, ...
        buttons: QMessageBox.Ok, QMessageBox.Cancel, ...
        title: string
        info: string
        detail: string
        return: QMessage's button code
    '''
    msg = QMessageBox()
    if text: msg.setText('%s'%text)
    if title: msg.setWindowTitle('%s'%title)
    if info: msg.setInformativeText('%s'%info)
    if detail: msg.setDetailedText('%s'%detail)
    if icon: msg.setIcon(icon)
    if buttons: msg.setStandardButtons(buttons)
    if default_button: msg.setDefaultButton(default_button)
    if escape_button: msg.setDefaultButton(escape_button)
    res = msg.exec_()
    return res


def info_dialog(*args, **kwargs):
    kwargs.update(dict(icon=QMessageBox.Information))
    return dialog(*args, **kwargs)

def question_dialog(*args, **kwargs):
    kwargs.update(dict(icon=QMessageBox.Question))
    return dialog(*args, **kwargs)

def confirm_dialog(*args, **kwargs):
    kwargs.update(dict(icon=QMessageBox.Question, buttons=QMessageBox.Yes|QMessageBox.No))
    return dialog(*args, **kwargs) == QMessageBox.Yes

def warning_dialog(*args, **kwargs):
    kwargs.update(dict(icon=QMessageBox.Warning))
    return dialog(*args, **kwargs)

def error_dialog(*args, **kwargs):
    kwargs.update(dict(icon=QMessageBox.Critical))
    return dialog(*args, **kwargs)

def exception_dialog(*args, **kwargs):
    kwargs.update(dict(icon=QMessageBox.Critical, detail=traceback.format_exc()))
    return dialog(*args, **kwargs)


def create_widget(spec, parent=None, getframe=False):
    '''
    Create a widget from the given ConfigObj specification.
    
    :Params:
        - **spec**: the widget specification:
            - **funcname**: describe the data type, one of:
                - **checkbox | flag | bool | boolean**: A CheckBox
                - **integer | int | long**: A SpinCtrl
                - **number | real | float**: A FloatSpin
                - **text | string | str | unicode**: a TextCtrl
                - **radio**: Radio buttons
                - **checkboxes | checklistbox | flags | bools | booleans**: A set of CheckListBox
                - **selector | select | combobox | combo | list | tuple'**: A ComboBox
            - **value|default**: the widget initial value or selection index (for radio and selector only) or indexes (for checkboxes only)
            - **args**: a list of choices (for radio, checkboxes and selector only)
            - **kwargs**: dict:
                - **min**: minimum value (for integer and float only)
                - **max**: maximum value (for integer and float only)
                - **increment**: increment value (float only)
                - **digits**: number of digits value (float only)
        - **parent**: widget's parent (qt parent)
        
        - **getframe**:
            special case such as file or directory types requires the creation of a 
            widget(text), a button and to group them into a frame in which case using getframe=True
            will return (frame, widget)
    
    :Return: the created qt widget
    '''
    
    widget_type = 'string' if spec.iterable or spec.funcname in ('dict',) else spec.funcname
    
    frame = None
    
    if widget_type in ('bool', 'boolean', 'checkbox'):
    
        widget = QtGui.QCheckBox(parent)
        widget.setTristate(False)
    
    # TODO: find acceptable setMinimum/setMaximum for QSpinBox/QDoubleSpinBox, qt default is 0-99
    # if out of acceptable range, qt may not behave as expected (sys.maxsize/sys.float_info.max are too large)
    elif widget_type in ('int', 'integer'):
    
        widget = QtGui.QSpinBox(parent)
        vmin = spec.kwargs.get('min', -1000000000)
        if vmin is not None:
            widget.setMinimum(int(vmin))
        vmax = spec.kwargs.get('max', 1000000000)
        if vmax is not None:
            widget.setMaximum(int(vmax))
        step = spec.kwargs.get('step', 1)
        if step is not None:
            widget.setSingleStep(int(step))
    
    elif widget_type in ('float', 'real'):
    
        widget = QtGui.QDoubleSpinBox(parent)
        vmin = spec.kwargs.get('min', -1000000000)
        if vmin is not None:
            widget.setMinimum(float(vmin))
        vmax = spec.kwargs.get('max', 1000000000)
        if vmax is not None:
            widget.setMaximum(float(vmax))
        step = spec.kwargs.get('step', 1)
        if step is not None:
            widget.setSingleStep(float(step))
        dec = spec.kwargs.get('decimals', 2)
        if dec is not None:
            widget.setDecimals(int(dec))
    
    elif widget_type in ('str', 'unicode', 'string', 'text'):
    
        widget = QtGui.QLineEdit(parent)
    
    elif widget_type in ('select', 'selector', 'combobox', 'combo', 'option', 'list', 'tuple'):
    
        widget = QtGui.QComboBox(parent)
        for item in spec.args:
            widget.addItem(item)
    
    elif widget_type in ('datetime',):#, 'date'):
    
        widget = QtGui.QDateTimeEdit()
        fmt = spec.kwargs.get('fmt', '%Y-%m-%d %H:%M:%S')
        widget.setCalendarPopup(True)
        widget.setDisplayFormat('yyyy-MM-dd hh:mm:ss')
    
    elif widget_type in ('path', 'file', 'directory'):
        select_file = widget_type in ('path', 'file')
        
        frame = QtGui.QFrame(parent)
        line = QtGui.QLineEdit(frame)
        button = QtGui.QPushButton("Browse")
        button.setIcon(QtGui.QIcon.fromTheme('document-open' if select_file else 'folder-open'))
        def clicked():
            if select_file:
                path = QtGui.QFileDialog.getOpenFileName(frame, directory=os.path.dirname(unicode(line.text())))
            else:
                path = QtGui.QFileDialog.getExistingDirectory(frame, directory=line.text())
            if path:
                line.setText(path)
        button.clicked.connect(clicked)
        layout = QtGui.QHBoxLayout(frame)
        layout.addWidget(line)
        layout.addWidget(button)
        
        widget = line
    
    else:
        logger.warning('Widget type is not supported %r, using it as raw text (spec: %s)', widget_type, pprint.pformat(spec))
        widget = QtGui.QLineEdit(parent)
        #raise TypeError("Failed to create a widget, unexpected parameter type: %s"%widget_type)
    
    set_widget(spec, widget, None)
    
    widget = (frame, widget) if getframe else widget
    
    return widget


def create_treeview_widget(spec, treeview, item):
    frame, widget = create_widget(spec, parent=treeview, getframe=True)
    widget.frame = frame
    if frame:
        treeview.setIndexWidget(item.index(), frame)
    else:
        treeview.setIndexWidget(item.index(), widget)
    return widget


def set_widget(spec, widget, value):
    
    if value is None:
        value = spec.get('value', spec.get('default', None))
    
    if isinstance(widget, QtGui.QCheckBox):
        widget.setCheckState(False if value is None else parse_bool(value))
        # XXX: force no tri-state, it seems the checkbox must have been added to a layout (or parent or ??) to take into account setTristate
        widget.setTristate(False)
    
    elif isinstance(widget, (QtGui.QSpinBox, QtGui.QDoubleSpinBox)):
        widget.setValue(0 if value in (None, 'None') else float(value))
    
    elif isinstance(widget, QtGui.QComboBox):
        index = widget.findText('' if value is None else value)
        widget.setCurrentIndex(index)
    
    elif isinstance(widget, QtGui.QLineEdit):
        # If this option type is a list, values will be shown and retrieved from UI as a string
        if spec.iterable:
            value = '' if value in (None, 'None') else ','.join(map(unicode, value))
        if spec.funcname == 'dict':
            value = repr(value)
        widget.setText('' if value is None else value)
    
    elif isinstance(widget, QtGui.QDateTimeEdit):
        if value in (None, 'None'):
            value = datetime.datetime.fromtimestamp(0)
        if not isinstance(value, datetime.datetime):
            fmt = spec.kwargs.get('fmt', '%Y-%m-%d %H:%M:%S')
            value = datetime.datetime.strptime(value, fmt)
        widget.setDateTime(value)
    
    else:
        raise TypeError("Failed to set widget value, unexpected type: %s"%type(widget))


def get_widget(spec, widget):
    
    if isinstance(widget, QtGui.QCheckBox):
        value = widget.isChecked() # != checkState
    
    elif isinstance(widget, (QtGui.QSpinBox, QtGui.QDoubleSpinBox)):
        value = widget.value()
    
    elif isinstance(widget, QtGui.QComboBox):
        value = unicode(widget.currentText())
    
    elif isinstance(widget, QtGui.QLineEdit):
        value = unicode(widget.text())
        if spec.funcname == 'dict':
            value = eval(value, {}, {})
    
    elif isinstance(widget, QtGui.QDateTimeEdit):
        fmt = spec.kwargs.get('fmt', '%Y-%m-%d %H:%M:%S')
        value = widget.dateTime().toPyDateTime().strftime(fmt)
    
    else:
        raise TypeError("Failed to get widget value, unexpected type: %s"%type(widget))
    
    # If this option type is a list, values will be shown and retrieved from UI as a string
    if spec.iterable:
        value = filter(None, map(lambda s: s.strip(), value.split(',')))
    
    return value


_true = ['y', 'yes', 'true', 'on', '1']
_false = ['n', 'no', 'false', 'off', '0']
def parse_bool(s):
    '''Boolean string parser'''
    _s, s = s, str(s).lower().strip()
    if s in _true:
        return True
    if s in _false:
        return False
    raise ValueError('Could not parse bool from string "'+_s+'"')

