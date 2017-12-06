# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'preferences.ui'
#
# Created: Mon Dec  4 17:28:37 2017
#      by: PyQt4 UI code generator 4.6.2
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

class Ui_PreferencesDialog(object):
    def setupUi(self, PreferencesDialog):
        PreferencesDialog.setObjectName("PreferencesDialog")
        PreferencesDialog.resize(399, 225)
        self.verticalLayout = QtGui.QVBoxLayout(PreferencesDialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.label_save_defaults = QtGui.QLabel(PreferencesDialog)
        self.label_save_defaults.setObjectName("label_save_defaults")
        self.gridLayout.addWidget(self.label_save_defaults, 0, 0, 1, 1)
        self.checkbox_save_defaults = QtGui.QCheckBox(PreferencesDialog)
        self.checkbox_save_defaults.setTristate(False)
        self.checkbox_save_defaults.setObjectName("checkbox_save_defaults")
        self.gridLayout.addWidget(self.checkbox_save_defaults, 0, 1, 1, 1)
        self.label_backup_count = QtGui.QLabel(PreferencesDialog)
        self.label_backup_count.setObjectName("label_backup_count")
        self.gridLayout.addWidget(self.label_backup_count, 2, 0, 1, 1)
        self.label_save_comments = QtGui.QLabel(PreferencesDialog)
        self.label_save_comments.setObjectName("label_save_comments")
        self.gridLayout.addWidget(self.label_save_comments, 1, 0, 1, 1)
        self.checkbox_save_comments = QtGui.QCheckBox(PreferencesDialog)
        self.checkbox_save_comments.setTristate(False)
        self.checkbox_save_comments.setObjectName("checkbox_save_comments")
        self.gridLayout.addWidget(self.checkbox_save_comments, 1, 1, 1, 1)
        self.spinbox_backup_count = QtGui.QSpinBox(PreferencesDialog)
        self.spinbox_backup_count.setObjectName("spinbox_backup_count")
        self.gridLayout.addWidget(self.spinbox_backup_count, 2, 1, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.buttonbox = QtGui.QDialogButtonBox(PreferencesDialog)
        self.buttonbox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonbox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonbox.setObjectName("buttonbox")
        self.verticalLayout.addWidget(self.buttonbox)

        self.retranslateUi(PreferencesDialog)
        QtCore.QObject.connect(self.buttonbox, QtCore.SIGNAL("accepted()"), PreferencesDialog.accept)
        QtCore.QObject.connect(self.buttonbox, QtCore.SIGNAL("rejected()"), PreferencesDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(PreferencesDialog)

    def retranslateUi(self, PreferencesDialog):
        PreferencesDialog.setWindowTitle(QtGui.QApplication.translate("PreferencesDialog", "Configuration", None, QtGui.QApplication.UnicodeUTF8))
        self.label_save_defaults.setText(QtGui.QApplication.translate("PreferencesDialog", "Save defaults", None, QtGui.QApplication.UnicodeUTF8))
        self.label_backup_count.setText(QtGui.QApplication.translate("PreferencesDialog", "Backup count", None, QtGui.QApplication.UnicodeUTF8))
        self.label_save_comments.setText(QtGui.QApplication.translate("PreferencesDialog", "Save comments", None, QtGui.QApplication.UnicodeUTF8))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    PreferencesDialog = QtGui.QDialog()
    ui = Ui_PreferencesDialog()
    ui.setupUi(PreferencesDialog)
    PreferencesDialog.show()
    sys.exit(app.exec_())

