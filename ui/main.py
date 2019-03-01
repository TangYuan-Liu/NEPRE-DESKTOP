# coding=utf-8

"""
Nepre GUI module
Establised by Siyuan Liu of LiuLab of Beijing Computational Science Research Center.
"""

from PyQt4.QtGui import *
from PyQt4.QtCore import *
from PyQt4.QtWebKit import *
import sys


class MainWindow(QWidget):

    def __init__(self):
        super(MainWindow,self).__init__()
        self.initUI()
    
    def initUI(self):
        self.tabWidget = QWidgets.QTabWidget()

class ContentWidget(QDialog):
    def __init__(self,parent=None):
        super(ContentWidget,self).__init__(parent)
        #self.setStyleSheet("background:royalblue")

class IndexWidget(QDialog):
    def __init__(self,parent=None):
        super(IndexWidget,self).__init__(parent)
        self.setStyleSheet("background:gray")

class nepreUI(QWidget):
    def __init__(self,parent=None):
        super(nepreUI,self).__init__()
        self.initUI()

    def initUI(self):
        # Nepre-F component
        self.txt1 = QLineEdit()
        self.txt1.setStyleSheet("background:linen")
        self.txt4 = QLineEdit()
        self.txt4.setStyleSheet("background:linen")
        self.button1 = QPushButton("Select PDB File")
        self.button1.setStyleSheet("background:linen")
        self.button4 = QPushButton("Go")
        self.button4.setStyleSheet("background:lightsteelblue")
        self.label1 = QLabel("Select Cutoff")
        self.button3 = QPushButton("Select Energy Matrix")
        self.button3.setStyleSheet("background:linen")
        self.checkbox1 = QCheckBox("Use Default Param")
        self.combobox0 = QComboBox()
        self.combobox0.setStyleSheet("background:linen")
        self.combobox0.addItem("6Å")
        self.combobox0.addItem("7Å")
        self.combobox0.addItem("8Å")
        self.combobox0.addItem("9Å")
        self.combobox0.addItem("10Å")
        
        # Nepre-R component
        self.txt2 = QLineEdit()
        self.txt2.setStyleSheet("background:wheat")
        self.txt5 = QLineEdit()
        self.txt5.setStyleSheet("background:wheat")
        self.txt6 = QLineEdit()
        self.txt6.setStyleSheet("background:wheat")
        self.checkbox2 = QCheckBox("Use Default Param")
        self.button2 = QPushButton("Select PDB File")
        self.button2.setStyleSheet("background:wheat")
        self.button5 = QPushButton("Select Energy Matrix")
        self.button5.setStyleSheet("background:wheat")
        self.button6 = QPushButton("Select Radius File")
        self.button6.setStyleSheet("background:wheat")
        self.button7 = QPushButton("Go")
        self.button7.setStyleSheet("background:lightsteelblue")
        # outlog
        self.txt3 = QTextEdit()

        ### Generate energy matrix component
        # Nepre-F energy matrix
        self.button8 = QPushButton("Select Dataset")
        self.button8.setStyleSheet("background:linen")
        self.txt7 = QLineEdit()
        self.txt7.setStyleSheet("background:linen")
        self.button9 = QPushButton("Save Path")
        self.button9.setStyleSheet("background:linen")
        self.txt8 = QLineEdit()
        self.txt8.setStyleSheet("background:linen")
        self.button13 = QPushButton("Go")
        self.button13.setStyleSheet("background:lightsteelblue")
        self.combobox1 = QComboBox()
        self.combobox1.setStyleSheet("background:linen")
        self.combobox1.addItem("6Å")
        self.combobox1.addItem("7Å")
        self.combobox1.addItem("8Å")
        self.combobox1.addItem("9Å")
        self.combobox1.addItem("10Å")

        # Nepre-R energy matrix
        self.button10 = QPushButton("Select Dataset")
        self.button10.setStyleSheet("background:wheat")
        self.txt9 = QLineEdit()
        self.txt9.setStyleSheet("background:wheat")
        self.button11 = QPushButton("Save Path")
        self.button11.setStyleSheet("background:wheat")
        self.txt10 = QLineEdit()
        self.txt10.setStyleSheet("background:wheat")
        self.button12 = QPushButton("Load Radius Parma")
        self.button12.setStyleSheet("background:wheat")
        self.txt11 = QLineEdit()
        self.txt11.setStyleSheet("background:wheat")
        self.button14 = QPushButton("Go")
        self.button14.setStyleSheet("background:lightsteelblue")

        ### define layout
        # left-right, first level
        layout = QHBoxLayout()
        # left part
        layout01 = QVBoxLayout()
        # Nepre part
        layout001 = QVBoxLayout()
        # Nepre-F part
        layout0001 = QHBoxLayout()
        layout0003 = QVBoxLayout()
        layout0004 = QHBoxLayout()
        layout0005 = QHBoxLayout()
        # Nepre-R part
        layout0002 = QHBoxLayout()
        layout0006 = QHBoxLayout()
        layout0007 = QHBoxLayout()
        layout0008 = QHBoxLayout()
        layout0009 = QVBoxLayout()
        # Show pics on assessment part
        layout02 = QVBoxLayout()

        # energy matrix for Nepre-F
        self.groupbox5 = QGroupBox("Nepre-F")
        layout0010 = QVBoxLayout()
        layout0011 = QHBoxLayout()
        layout0012 = QHBoxLayout()
        layout0013 = QHBoxLayout()

        # energy matrix for Nepre-R
        self.groupbox6 = QGroupBox("Nepre-R")
        layout0014 = QVBoxLayout()
        layout0015 = QHBoxLayout()
        layout0016 = QHBoxLayout()
        layout0017 = QHBoxLayout()
        layout0018 = QHBoxLayout()

        ### assemble layout
        # Nepre-F part
        self.groupbox0 = QGroupBox("Nepre-F")
        self.groupbox0.setStyleSheet("background:lavender")
        layout0001.addWidget(self.button1)
        layout0001.addWidget(self.txt1)
        layout0004.addWidget(self.button3)
        layout0004.addWidget(self.txt4)
        layout0004.addWidget(self.checkbox1)

        layout0005.addWidget(self.combobox0)
        layout0005.addSpacing(200)
        layout0005.addWidget(self.button4)
        layout0003.addLayout(layout0001)
        layout0003.addLayout(layout0004)
        layout0003.addLayout(layout0005)
        self.groupbox0.setLayout(layout0003)

        # Nepre-R part
        self.groupbox1 = QGroupBox("Nepre-R")
        layout0002.addWidget(self.button2)
        layout0002.addWidget(self.txt2)
        layout0006.addWidget(self.button5)
        layout0006.addWidget(self.txt5)
        layout0006.addWidget(self.checkbox2)
        layout0007.addWidget(self.button6)
        layout0007.addWidget(self.txt6)
        layout0008.addSpacing(395)
        layout0008.addWidget(self.button7)
        layout0009.addLayout(layout0002)
        layout0009.addLayout(layout0006)
        layout0009.addLayout(layout0007)
        layout0009.addLayout(layout0008)

        self.groupbox1.setLayout(layout0009)
        self.groupbox1.setStyleSheet("background:lavenderblush")

        # Generate energy matrix
        self.groupbox3 = QGroupBox("Generate Energy Matrix")
        layout002 = QVBoxLayout()
        ## Nepre-F
        layout0011.addWidget(self.button8)
        layout0011.addWidget(self.txt7)
        layout0012.addWidget(self.button9)
        layout0012.addWidget(self.txt8)
        layout0013.addWidget(self.combobox1)
        layout0013.addWidget(self.button13)
        layout0010.addLayout(layout0011)
        layout0010.addLayout(layout0012)
        layout0010.addLayout(layout0013)
        self.groupbox5.setLayout(layout0010)
        self.groupbox5.setStyleSheet("background:lavender")

        ## Nepre-R
        layout0015.addWidget(self.button10)
        layout0015.addWidget(self.txt9)
        layout0016.addWidget(self.button11)
        layout0016.addWidget(self.txt10)
        layout0017.addWidget(self.button12)
        layout0017.addWidget(self.txt11)
        layout0018.addSpacing(200)
        layout0018.addWidget(self.button14)
        layout0014.addLayout(layout0015)
        layout0014.addLayout(layout0016)
        layout0014.addLayout(layout0017)
        layout0014.addLayout(layout0018)
        self.groupbox6.setLayout(layout0014)
        self.groupbox6.setStyleSheet("background:lavenderblush")

        # Show pics on assessment part
        self.pixmap = QPixmap("../pics/1.png")
        self.lbl1 = QLabel()
        self.lbl1.setPixmap(self.pixmap)

        # Console part
        self.groupbox4 = QGroupBox("Console")
        layout020 = QHBoxLayout()

        # Assemble layout
        layout001.addWidget(self.groupbox0)
        layout001.addWidget(self.groupbox1)
        layout01.addLayout(layout001)
        layout002.addWidget(self.groupbox5)
        layout002.addWidget(self.groupbox6)
        self.groupbox3.setLayout(layout002)
        layout01.addWidget(self.groupbox3)
        layout020.addWidget(self.txt3)
        self.groupbox4.setLayout(layout020)
        layout02.addWidget(self.lbl1)
        layout02.addWidget(self.groupbox4)
        layout.addLayout(layout01)
        layout.addLayout(layout02)
        self.setLayout(layout)
        
class DataProcess(QWidget):
    def __init__(self,parent=None):
        super(DataProcess,self).__init__()
        self.InitUI()
    
    def InitUI(self):
        # define component
        self.button15 = QPushButton()
        self.button16 = QPushButton()
        self.button17 = QPushButton()

class InstucAndInfo(QWidget):
    def __init__(self,parent=None):
        super(InstucAndInfo,self).__init__()
        self.InitUI()
    
    def InitUI(self):
        # define web viewer component
        self.browser = QWebView()
        self.browser.load(QUrl("./method.html"))
        # define layout
        layout = QHBoxLayout()
        layout.addWidget(self.browser)
        self.setLayout(layout)
    

class TabWidget(QTabWidget):
    def __init__(self,parent=None):
        super(TabWidget,self).__init__(parent)
        self.resize(1000,900)
        self.nepre = nepreUI()
        self.dataprocess = DataProcess()
        self.info = InstucAndInfo()
        self.addTab(self.nepre,"Structure Assessment")
        self.addTab(self.dataprocess, "Data Analyze")
        self.addTab(self.info, "Instruction && Information")
        title = "NEPRE--Scoring Function Based on Neighbourhood Preference Statistics"
        self.setWindowTitle(title)
        



if __name__ == "__main__":
    app = QApplication(sys.argv)
    t = TabWidget()
    t.show()
    app.exec_()

