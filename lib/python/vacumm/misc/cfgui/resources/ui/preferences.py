# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'preferences.ui'
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

class Ui_PreferencesDialog(object):
    def setupUi(self, PreferencesDialog):
        PreferencesDialog.setObjectName(_fromUtf8("PreferencesDialog"))
        PreferencesDialog.resize(399, 225)
        self.verticalLayout = QtGui.QVBoxLayout(PreferencesDialog)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label_save_defaults = QtGui.QLabel(PreferencesDialog)
        self.label_save_defaults.setObjectName(_fromUtf8("label_save_defaults"))
        self.gridLayout.addWidget(self.label_save_defaults, 0, 0, 1, 1)
        self.checkbox_save_defaults = QtGui.QCheckBox(PreferencesDialog)
        self.checkbox_save_defaults.setText(_fromUtf8(""))
        self.checkbox_save_defaults.setTristate(False)
        self.checkbox_save_defaults.setObjectName(_fromUtf8("checkbox_save_defaults"))
        self.gridLayout.addWidget(self.checkbox_save_defaults, 0, 1, 1, 1)
        self.label_backup_count = QtGui.QLabel(PreferencesDialog)
        self.label_backup_count.setObjectName(_fromUtf8("label_backup_count"))
        self.gridLayout.addWidget(self.label_backup_count, 2, 0, 1, 1)
        self.label_save_comments = QtGui.QLabel(PreferencesDialog)
        self.label_save_comments.setObjectName(_fromUtf8("label_save_comments"))
        self.gridLayout.addWidget(self.label_save_comments, 1, 0, 1, 1)
        self.checkbox_save_comments = QtGui.QCheckBox(PreferencesDialog)
        self.checkbox_save_comments.setText(_fromUtf8(""))
        self.checkbox_save_comments.setTristate(False)
        self.checkbox_save_comments.setObjectName(_fromUtf8("checkbox_save_comments"))
        self.gridLayout.addWidget(self.checkbox_save_comments, 1, 1, 1, 1)
        self.spinbox_backup_count = QtGui.QSpinBox(PreferencesDialog)
        self.spinbox_backup_count.setObjectName(_fromUtf8("spinbox_backup_count"))
        self.gridLayout.addWidget(self.spinbox_backup_count, 2, 1, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.buttonbox = QtGui.QDialogButtonBox(PreferencesDialog)
        self.buttonbox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonbox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonbox.setObjectName(_fromUtf8("buttonbox"))
        self.verticalLayout.addWidget(self.buttonbox)

        self.retranslateUi(PreferencesDialog)
        QtCore.QObject.connect(self.buttonbox, QtCore.SIGNAL(_fromUtf8("accepted()")), PreferencesDialog.accept)
        QtCore.QObject.connect(self.buttonbox, QtCore.SIGNAL(_fromUtf8("rejected()")), PreferencesDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(PreferencesDialog)

    def retranslateUi(self, PreferencesDialog):
        PreferencesDialog.setWindowTitle(_translate("PreferencesDialog", "Configuration", None))
        self.label_save_defaults.setText(_translate("PreferencesDialog", "Save defaults", None))
        self.label_backup_count.setText(_translate("PreferencesDialog", "Backup count", None))
        self.label_save_comments.setText(_translate("PreferencesDialog", "Save comments", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    PreferencesDialog = QtGui.QDialog()
    ui = Ui_PreferencesDialog()
    ui.setupUi(PreferencesDialog)
    PreferencesDialog.show()
    sys.exit(app.exec_())

