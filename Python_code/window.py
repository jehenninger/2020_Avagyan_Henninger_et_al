# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/Users/jon/PycharmProjects/zbow/bin/window.ui'
#
# Created by: PyQt5 UI code generator 5.11.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1271, 856)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.centralwidget.sizePolicy().hasHeightForWidth())
        self.centralwidget.setSizePolicy(sizePolicy)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_15 = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_15.sizePolicy().hasHeightForWidth())
        self.label_15.setSizePolicy(sizePolicy)
        self.label_15.setObjectName("label_15")
        self.verticalLayout.addWidget(self.label_15)
        self.fileLabel = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.fileLabel.sizePolicy().hasHeightForWidth())
        self.fileLabel.setSizePolicy(sizePolicy)
        self.fileLabel.setText("")
        self.fileLabel.setObjectName("fileLabel")
        self.verticalLayout.addWidget(self.fileLabel)
        self.viewOutliersButton = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.viewOutliersButton.sizePolicy().hasHeightForWidth())
        self.viewOutliersButton.setSizePolicy(sizePolicy)
        self.viewOutliersButton.setObjectName("viewOutliersButton")
        self.verticalLayout.addWidget(self.viewOutliersButton)
        self.removeOutliers = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.removeOutliers.sizePolicy().hasHeightForWidth())
        self.removeOutliers.setSizePolicy(sizePolicy)
        self.removeOutliers.setObjectName("removeOutliers")
        self.verticalLayout.addWidget(self.removeOutliers)
        self.parameterTable = QtWidgets.QTableWidget(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.parameterTable.sizePolicy().hasHeightForWidth())
        self.parameterTable.setSizePolicy(sizePolicy)
        self.parameterTable.setShowGrid(False)
        self.parameterTable.setObjectName("parameterTable")
        self.parameterTable.setColumnCount(0)
        self.parameterTable.setRowCount(0)
        self.parameterTable.horizontalHeader().setStretchLastSection(True)
        self.parameterTable.verticalHeader().setVisible(False)
        self.parameterTable.verticalHeader().setStretchLastSection(False)
        self.verticalLayout.addWidget(self.parameterTable)
        self.updateParams = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.updateParams.sizePolicy().hasHeightForWidth())
        self.updateParams.setSizePolicy(sizePolicy)
        self.updateParams.setObjectName("updateParams")
        self.verticalLayout.addWidget(self.updateParams)
        self.horizontalLayout.addWidget(self.groupBox)
        self.groupBox_4 = QtWidgets.QGroupBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_4.sizePolicy().hasHeightForWidth())
        self.groupBox_4.setSizePolicy(sizePolicy)
        self.groupBox_4.setTitle("")
        self.groupBox_4.setObjectName("groupBox_4")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_4)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.groupBox_2 = QtWidgets.QGroupBox(self.groupBox_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_2.sizePolicy().hasHeightForWidth())
        self.groupBox_2.setSizePolicy(sizePolicy)
        self.groupBox_2.setAutoFillBackground(True)
        self.groupBox_2.setTitle("")
        self.groupBox_2.setObjectName("groupBox_2")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.groupBox_2)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.clusterOptionsBox = QtWidgets.QGroupBox(self.groupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.clusterOptionsBox.sizePolicy().hasHeightForWidth())
        self.clusterOptionsBox.setSizePolicy(sizePolicy)
        self.clusterOptionsBox.setObjectName("clusterOptionsBox")
        self.formLayout = QtWidgets.QFormLayout(self.clusterOptionsBox)
        self.formLayout.setObjectName("formLayout")
        self.label_2 = QtWidgets.QLabel(self.clusterOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        self.label_2.setScaledContents(True)
        self.label_2.setObjectName("label_2")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_2)
        self.clusterSampleSize = QtWidgets.QLineEdit(self.clusterOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.clusterSampleSize.sizePolicy().hasHeightForWidth())
        self.clusterSampleSize.setSizePolicy(sizePolicy)
        self.clusterSampleSize.setObjectName("clusterSampleSize")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.clusterSampleSize)
        self.label_3 = QtWidgets.QLabel(self.clusterOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy)
        self.label_3.setScaledContents(True)
        self.label_3.setObjectName("label_3")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_3)
        self.clusterOnData = QtWidgets.QComboBox(self.clusterOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.clusterOnData.sizePolicy().hasHeightForWidth())
        self.clusterOnData.setSizePolicy(sizePolicy)
        self.clusterOnData.setObjectName("clusterOnData")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.clusterOnData)
        self.label_4 = QtWidgets.QLabel(self.clusterOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_4.sizePolicy().hasHeightForWidth())
        self.label_4.setSizePolicy(sizePolicy)
        self.label_4.setScaledContents(True)
        self.label_4.setObjectName("label_4")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.label_4)
        self.clusterMinClusterSize = QtWidgets.QLineEdit(self.clusterOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.clusterMinClusterSize.sizePolicy().hasHeightForWidth())
        self.clusterMinClusterSize.setSizePolicy(sizePolicy)
        self.clusterMinClusterSize.setObjectName("clusterMinClusterSize")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.clusterMinClusterSize)
        self.label_5 = QtWidgets.QLabel(self.clusterOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_5.sizePolicy().hasHeightForWidth())
        self.label_5.setSizePolicy(sizePolicy)
        self.label_5.setScaledContents(True)
        self.label_5.setObjectName("label_5")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.label_5)
        self.clusterMinSamples = QtWidgets.QLineEdit(self.clusterOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.clusterMinSamples.sizePolicy().hasHeightForWidth())
        self.clusterMinSamples.setSizePolicy(sizePolicy)
        self.clusterMinSamples.setObjectName("clusterMinSamples")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.clusterMinSamples)
        self.clusterPushButton = QtWidgets.QPushButton(self.clusterOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.clusterPushButton.sizePolicy().hasHeightForWidth())
        self.clusterPushButton.setSizePolicy(sizePolicy)
        self.clusterPushButton.setObjectName("clusterPushButton")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.clusterPushButton)
        self.evaluateClusteringCheckBox = QtWidgets.QCheckBox(self.clusterOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.evaluateClusteringCheckBox.sizePolicy().hasHeightForWidth())
        self.evaluateClusteringCheckBox.setSizePolicy(sizePolicy)
        self.evaluateClusteringCheckBox.setObjectName("evaluateClusteringCheckBox")
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.evaluateClusteringCheckBox)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.formLayout.setItem(7, QtWidgets.QFormLayout.FieldRole, spacerItem)
        self.highlightClusterPushButton = QtWidgets.QPushButton(self.clusterOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.highlightClusterPushButton.sizePolicy().hasHeightForWidth())
        self.highlightClusterPushButton.setSizePolicy(sizePolicy)
        self.highlightClusterPushButton.setObjectName("highlightClusterPushButton")
        self.formLayout.setWidget(9, QtWidgets.QFormLayout.FieldRole, self.highlightClusterPushButton)
        self.joinClusterPushButton = QtWidgets.QPushButton(self.clusterOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.joinClusterPushButton.sizePolicy().hasHeightForWidth())
        self.joinClusterPushButton.setSizePolicy(sizePolicy)
        self.joinClusterPushButton.setObjectName("joinClusterPushButton")
        self.formLayout.setWidget(10, QtWidgets.QFormLayout.FieldRole, self.joinClusterPushButton)
        self.splitCluster = QtWidgets.QPushButton(self.clusterOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.splitCluster.sizePolicy().hasHeightForWidth())
        self.splitCluster.setSizePolicy(sizePolicy)
        self.splitCluster.setObjectName("splitCluster")
        self.formLayout.setWidget(11, QtWidgets.QFormLayout.FieldRole, self.splitCluster)
        self.horizontalLayout_2.addWidget(self.clusterOptionsBox)
        self.plotOptionsBox = QtWidgets.QGroupBox(self.groupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plotOptionsBox.sizePolicy().hasHeightForWidth())
        self.plotOptionsBox.setSizePolicy(sizePolicy)
        self.plotOptionsBox.setObjectName("plotOptionsBox")
        self.formLayout_2 = QtWidgets.QFormLayout(self.plotOptionsBox)
        self.formLayout_2.setObjectName("formLayout_2")
        self.label_6 = QtWidgets.QLabel(self.plotOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_6.sizePolicy().hasHeightForWidth())
        self.label_6.setSizePolicy(sizePolicy)
        self.label_6.setScaledContents(True)
        self.label_6.setObjectName("label_6")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_6)
        self.label_7 = QtWidgets.QLabel(self.plotOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_7.sizePolicy().hasHeightForWidth())
        self.label_7.setSizePolicy(sizePolicy)
        self.label_7.setScaledContents(True)
        self.label_7.setObjectName("label_7")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_7)
        self.scatterColorOption = QtWidgets.QComboBox(self.plotOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.scatterColorOption.sizePolicy().hasHeightForWidth())
        self.scatterColorOption.setSizePolicy(sizePolicy)
        self.scatterColorOption.setObjectName("scatterColorOption")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.scatterColorOption)
        self.label_8 = QtWidgets.QLabel(self.plotOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_8.sizePolicy().hasHeightForWidth())
        self.label_8.setSizePolicy(sizePolicy)
        self.label_8.setScaledContents(True)
        self.label_8.setObjectName("label_8")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_8)
        self.scatterScaleOption = QtWidgets.QComboBox(self.plotOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.scatterScaleOption.sizePolicy().hasHeightForWidth())
        self.scatterScaleOption.setSizePolicy(sizePolicy)
        self.scatterScaleOption.setObjectName("scatterScaleOption")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.scatterScaleOption)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.formLayout_2.setItem(3, QtWidgets.QFormLayout.LabelRole, spacerItem1)
        self.label_9 = QtWidgets.QLabel(self.plotOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_9.sizePolicy().hasHeightForWidth())
        self.label_9.setSizePolicy(sizePolicy)
        self.label_9.setScaledContents(True)
        self.label_9.setObjectName("label_9")
        self.formLayout_2.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.label_9)
        self.label_10 = QtWidgets.QLabel(self.plotOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_10.sizePolicy().hasHeightForWidth())
        self.label_10.setSizePolicy(sizePolicy)
        self.label_10.setScaledContents(True)
        self.label_10.setObjectName("label_10")
        self.formLayout_2.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.label_10)
        self.ternColorOption = QtWidgets.QComboBox(self.plotOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ternColorOption.sizePolicy().hasHeightForWidth())
        self.ternColorOption.setSizePolicy(sizePolicy)
        self.ternColorOption.setObjectName("ternColorOption")
        self.formLayout_2.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.ternColorOption)
        self.label_11 = QtWidgets.QLabel(self.plotOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_11.sizePolicy().hasHeightForWidth())
        self.label_11.setSizePolicy(sizePolicy)
        self.label_11.setScaledContents(True)
        self.label_11.setObjectName("label_11")
        self.formLayout_2.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.label_11)
        self.ternScaleOption = QtWidgets.QComboBox(self.plotOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ternScaleOption.sizePolicy().hasHeightForWidth())
        self.ternScaleOption.setSizePolicy(sizePolicy)
        self.ternScaleOption.setObjectName("ternScaleOption")
        self.formLayout_2.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.ternScaleOption)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.formLayout_2.setItem(7, QtWidgets.QFormLayout.LabelRole, spacerItem2)
        self.label = QtWidgets.QLabel(self.plotOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setScaledContents(True)
        self.label.setObjectName("label")
        self.formLayout_2.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.label)
        self.giniCoeff = QtWidgets.QTextBrowser(self.plotOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.giniCoeff.sizePolicy().hasHeightForWidth())
        self.giniCoeff.setSizePolicy(sizePolicy)
        self.giniCoeff.setMaximumSize(QtCore.QSize(100, 30))
        self.giniCoeff.setObjectName("giniCoeff")
        self.formLayout_2.setWidget(8, QtWidgets.QFormLayout.FieldRole, self.giniCoeff)
        self.label_12 = QtWidgets.QLabel(self.plotOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_12.sizePolicy().hasHeightForWidth())
        self.label_12.setSizePolicy(sizePolicy)
        self.label_12.setScaledContents(True)
        self.label_12.setObjectName("label_12")
        self.formLayout_2.setWidget(9, QtWidgets.QFormLayout.LabelRole, self.label_12)
        self.shannonEntropy = QtWidgets.QTextBrowser(self.plotOptionsBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.shannonEntropy.sizePolicy().hasHeightForWidth())
        self.shannonEntropy.setSizePolicy(sizePolicy)
        self.shannonEntropy.setMaximumSize(QtCore.QSize(100, 30))
        self.shannonEntropy.setObjectName("shannonEntropy")
        self.formLayout_2.setWidget(9, QtWidgets.QFormLayout.FieldRole, self.shannonEntropy)
        self.horizontalLayout_2.addWidget(self.plotOptionsBox)
        self.verticalLayout_2.addWidget(self.groupBox_2)
        self.clusterInformationBox = QtWidgets.QGroupBox(self.groupBox_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.clusterInformationBox.sizePolicy().hasHeightForWidth())
        self.clusterInformationBox.setSizePolicy(sizePolicy)
        self.clusterInformationBox.setObjectName("clusterInformationBox")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.clusterInformationBox)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.clusterInformationTable = QtWidgets.QTableWidget(self.clusterInformationBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.clusterInformationTable.sizePolicy().hasHeightForWidth())
        self.clusterInformationTable.setSizePolicy(sizePolicy)
        self.clusterInformationTable.setEditTriggers(QtWidgets.QAbstractItemView.DoubleClicked)
        self.clusterInformationTable.setShowGrid(False)
        self.clusterInformationTable.setObjectName("clusterInformationTable")
        self.clusterInformationTable.setColumnCount(0)
        self.clusterInformationTable.setRowCount(0)
        self.clusterInformationTable.horizontalHeader().setStretchLastSection(True)
        self.clusterInformationTable.verticalHeader().setVisible(False)
        self.clusterInformationTable.verticalHeader().setStretchLastSection(False)
        self.verticalLayout_3.addWidget(self.clusterInformationTable)
        self.verticalLayout_2.addWidget(self.clusterInformationBox)
        self.progressBar = QtWidgets.QProgressBar(self.groupBox_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.progressBar.sizePolicy().hasHeightForWidth())
        self.progressBar.setSizePolicy(sizePolicy)
        self.progressBar.setProperty("value", 24)
        self.progressBar.setObjectName("progressBar")
        self.verticalLayout_2.addWidget(self.progressBar)
        self.horizontalLayout.addWidget(self.groupBox_4)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1271, 22))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.menubar.sizePolicy().hasHeightForWidth())
        self.menubar.setSizePolicy(sizePolicy)
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuData_processing = QtWidgets.QMenu(self.menubar)
        self.menuData_processing.setObjectName("menuData_processing")
        self.menuTools = QtWidgets.QMenu(self.menubar)
        self.menuTools.setObjectName("menuTools")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionLoad = QtWidgets.QAction(MainWindow)
        self.actionLoad.setObjectName("actionLoad")
        self.actionSave = QtWidgets.QAction(MainWindow)
        self.actionSave.setObjectName("actionSave")
        self.actionClear = QtWidgets.QAction(MainWindow)
        self.actionClear.setObjectName("actionClear")
        self.actionRestore = QtWidgets.QAction(MainWindow)
        self.actionRestore.setObjectName("actionRestore")
        self.actionQuit = QtWidgets.QAction(MainWindow)
        self.actionQuit.setObjectName("actionQuit")
        self.actionLoadClusteringSolution = QtWidgets.QAction(MainWindow)
        self.actionLoadClusteringSolution.setObjectName("actionLoadClusteringSolution")
        self.actionMake_cluster_plots = QtWidgets.QAction(MainWindow)
        self.actionMake_cluster_plots.setObjectName("actionMake_cluster_plots")
        self.menuFile.addAction(self.actionLoad)
        self.menuFile.addAction(self.actionLoadClusteringSolution)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionSave)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionClear)
        self.menuFile.addSeparator()
        self.menuData_processing.addAction(self.actionRestore)
        self.menuTools.addAction(self.actionMake_cluster_plots)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuData_processing.menuAction())
        self.menubar.addAction(self.menuTools.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "zbow analysis"))
        self.groupBox.setTitle(_translate("MainWindow", "Data info"))
        self.label_15.setText(_translate("MainWindow", "Loaded file"))
        self.viewOutliersButton.setText(_translate("MainWindow", "view outliers"))
        self.removeOutliers.setText(_translate("MainWindow", "remove outliers"))
        self.updateParams.setText(_translate("MainWindow", "update parameters"))
        self.clusterOptionsBox.setTitle(_translate("MainWindow", "cluster options"))
        self.label_2.setText(_translate("MainWindow", "sample size"))
        self.label_3.setText(_translate("MainWindow", "cluster data"))
        self.label_4.setText(_translate("MainWindow", "HDBSCAN min cluster size"))
        self.label_5.setText(_translate("MainWindow", "HDBSCAN min samples"))
        self.clusterPushButton.setText(_translate("MainWindow", "cluster"))
        self.evaluateClusteringCheckBox.setText(_translate("MainWindow", "Evaluate clustering (slower)"))
        self.highlightClusterPushButton.setText(_translate("MainWindow", "highlight cluster"))
        self.joinClusterPushButton.setText(_translate("MainWindow", "join clusters"))
        self.splitCluster.setText(_translate("MainWindow", "split cluster"))
        self.plotOptionsBox.setTitle(_translate("MainWindow", "plot options"))
        self.label_6.setText(_translate("MainWindow", "3D scatter options"))
        self.label_7.setText(_translate("MainWindow", "color"))
        self.label_8.setText(_translate("MainWindow", "scale"))
        self.label_9.setText(_translate("MainWindow", "ternary plot options"))
        self.label_10.setText(_translate("MainWindow", "color"))
        self.label_11.setText(_translate("MainWindow", "scale"))
        self.label.setText(_translate("MainWindow", "Gini coefficient"))
        self.label_12.setText(_translate("MainWindow", "Shannon Entropy"))
        self.clusterInformationBox.setTitle(_translate("MainWindow", "cluster information"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuData_processing.setTitle(_translate("MainWindow", "Data processing"))
        self.menuTools.setTitle(_translate("MainWindow", "Tools"))
        self.actionLoad.setText(_translate("MainWindow", "Load data"))
        self.actionLoad.setShortcut(_translate("MainWindow", "Ctrl+L"))
        self.actionSave.setText(_translate("MainWindow", "Save data/images"))
        self.actionSave.setShortcut(_translate("MainWindow", "Ctrl+S"))
        self.actionClear.setText(_translate("MainWindow", "Clear session"))
        self.actionClear.setShortcut(_translate("MainWindow", "Ctrl+R"))
        self.actionRestore.setText(_translate("MainWindow", "Restore original data"))
        self.actionQuit.setText(_translate("MainWindow", "Quit"))
        self.actionQuit.setShortcut(_translate("MainWindow", "Ctrl+C"))
        self.actionLoadClusteringSolution.setText(_translate("MainWindow", "Load clustering solution"))
        self.actionMake_cluster_plots.setText(_translate("MainWindow", "Make cluster plots"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

