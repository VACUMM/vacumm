# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'sessions.ui'
#
# Created: Wed Jul 26 10:05:58 2017
#      by: PyQt4 UI code generator 4.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_SessionsDialog(object):
    def setupUi(self, SessionsDialog):
        SessionsDialog.setObjectName(_fromUtf8("SessionsDialog"))
        SessionsDialog.resize(608, 203)
        self.verticalLayout = QtGui.QVBoxLayout(SessionsDialog)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label_name = QtGui.QLabel(SessionsDialog)
        self.label_name.setObjectName(_fromUtf8("label_name"))
        self.gridLayout.addWidget(self.label_name, 0, 0, 1, 1)
        self.label_specification = QtGui.QLabel(SessionsDialog)
        self.label_specification.setObjectName(_fromUtf8("label_specification"))
        self.gridLayout.addWidget(self.label_specification, 1, 0, 1, 1)
        self.line_specification = QtGui.QLineEdit(SessionsDialog)
        self.line_specification.setObjectName(_fromUtf8("line_specification"))
        self.gridLayout.addWidget(self.line_specification, 1, 1, 1, 1)
        self.line_configuration = QtGui.QLineEdit(SessionsDialog)
        self.line_configuration.setObjectName(_fromUtf8("line_configuration"))
        self.gridLayout.addWidget(self.line_configuration, 2, 1, 1, 1)
        self.label_configuration = QtGui.QLabel(SessionsDialog)
        self.label_configuration.setObjectName(_fromUtf8("label_configuration"))
        self.gridLayout.addWidget(self.label_configuration, 2, 0, 1, 1)
        self.combo_name = QtGui.QComboBox(SessionsDialog)
        self.combo_name.setEditable(True)
        self.combo_name.setObjectName(_fromUtf8("combo_name"))
        self.gridLayout.addWidget(self.combo_name, 0, 1, 1, 1)
        self.button_specification = QtGui.QPushButton(SessionsDialog)
        icon = QtGui.QIcon.fromTheme(_fromUtf8("document-open"))
        self.button_specification.setIcon(icon)
        self.button_specification.setObjectName(_fromUtf8("button_specification"))
        self.gridLayout.addWidget(self.button_specification, 1, 2, 1, 1)
        self.button_configuration = QtGui.QPushButton(SessionsDialog)
        icon = QtGui.QIcon.fromTheme(_fromUtf8("document-open"))
        self.button_configuration.setIcon(icon)
        self.button_configuration.setObjectName(_fromUtf8("button_configuration"))
        self.gridLayout.addWidget(self.button_configuration, 2, 2, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.button_new = QtGui.QPushButton(SessionsDialog)
        icon = QtGui.QIcon.fromTheme(_fromUtf8("document-new"))
        self.button_new.setIcon(icon)
        self.button_new.setObjectName(_fromUtf8("button_new"))
        self.horizontalLayout_2.addWidget(self.button_new)
        self.button_save = QtGui.QPushButton(SessionsDialog)
        icon = QtGui.QIcon.fromTheme(_fromUtf8("document-save"))
        self.button_save.setIcon(icon)
        self.button_save.setObjectName(_fromUtf8("button_save"))
        self.horizontalLayout_2.addWidget(self.button_save)
        self.button_open = QtGui.QPushButton(SessionsDialog)
        icon = QtGui.QIcon.fromTheme(_fromUtf8("document-open"))
        self.button_open.setIcon(icon)
        self.button_open.setObjectName(_fromUtf8("button_open"))
        self.horizontalLayout_2.addWidget(self.button_open)
        self.button_delete = QtGui.QPushButton(SessionsDialog)
        icon = QtGui.QIcon.fromTheme(_fromUtf8("edit-delete"))
        self.button_delete.setIcon(icon)
        self.button_delete.setObjectName(_fromUtf8("button_delete"))
        self.horizontalLayout_2.addWidget(self.button_delete)
        self.button_close = QtGui.QPushButton(SessionsDialog)
        icon = QtGui.QIcon.fromTheme(_fromUtf8("window-close"))
        self.button_close.setIcon(icon)
        self.button_close.setObjectName(_fromUtf8("button_close"))
        self.horizontalLayout_2.addWidget(self.button_close)
        self.verticalLayout.addLayout(self.horizontalLayout_2)

        self.retranslateUi(SessionsDialog)
        QtCore.QObject.connect(self.button_specification, QtCore.SIGNAL(_fromUtf8("clicked()")), SessionsDialog.on_button_specification)
        QtCore.QObject.connect(self.button_configuration, QtCore.SIGNAL(_fromUtf8("clicked()")), SessionsDialog.on_button_configuration)
        QtCore.QObject.connect(self.combo_name, QtCore.SIGNAL(_fromUtf8("activated(QString)")), SessionsDialog.on_combo_name)
        QtCore.QObject.connect(self.button_new, QtCore.SIGNAL(_fromUtf8("clicked()")), SessionsDialog.on_button_new)
        QtCore.QObject.connect(self.button_save, QtCore.SIGNAL(_fromUtf8("clicked()")), SessionsDialog.on_button_save)
        QtCore.QObject.connect(self.button_open, QtCore.SIGNAL(_fromUtf8("clicked()")), SessionsDialog.on_button_open)
        QtCore.QObject.connect(self.button_delete, QtCore.SIGNAL(_fromUtf8("clicked()")), SessionsDialog.on_button_delete)
        QtCore.QObject.connect(self.button_close, QtCore.SIGNAL(_fromUtf8("clicked()")), SessionsDialog.on_button_close)
        QtCore.QMetaObject.connectSlotsByName(SessionsDialog)

    def retranslateUi(self, SessionsDialog):
        SessionsDialog.setWindowTitle(_translate("SessionsDialog", "Sessions", None))
        self.label_name.setText(_translate("SessionsDialog", "Session name", None))
        self.label_specification.setText(_translate("SessionsDialog", "Specification file", None))
        self.label_configuration.setText(_translate("SessionsDialog", "Configuration file", None))
        self.button_specification.setText(_translate("SessionsDialog", "Browse", None))
        self.button_configuration.setText(_translate("SessionsDialog", "Browse", None))
        self.button_new.setText(_translate("SessionsDialog", "New", None))
        self.button_new.setShortcut(_translate("SessionsDialog", "Ctrl+N", None))
        self.button_save.setText(_translate("SessionsDialog", "Save", None))
        self.button_save.setShortcut(_translate("SessionsDialog", "Ctrl+S", None))
        self.button_open.setText(_translate("SessionsDialog", "Open", None))
        self.button_open.setShortcut(_translate("SessionsDialog", "Ctrl+O", None))
        self.button_delete.setText(_translate("SessionsDialog", "Delete", None))
        self.button_close.setText(_translate("SessionsDialog", "Close", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    SessionsDialog = QtGui.QDialog()
    ui = Ui_SessionsDialog()
    ui.setupUi(SessionsDialog)
    SessionsDialog.show()
    sys.exit(app.exec_())

