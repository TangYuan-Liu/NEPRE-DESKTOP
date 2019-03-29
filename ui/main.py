# coding=utf-8

"""
Nepre GUI module
Establised by Siyuan Liu of LiuLab of Beijing Computational Science Research Center.
"""

from PyQt4.QtGui import *
from PyQt4.QtCore import *
from PyQt4.QtWebKit import *
import sys
import os
sys.path.append("..")
import backend.package_func as tools

###################
# print redirection
class OutLog:
    def __init__(self, edit, out=None, color=None):
        """
        edit = QTextEdit
        out = alternate stream(can be the original sys.stdout)
        color = alternate color (color stderr a different color)
        """
        self.edit = edit
        self.out = None
        self.color = color

    def write(self, m):
        if self.color:
            tc = self.edit.textColor()
            self.edit.setTextColor(self.color)
        
        self.edit.moveCursor(QTextCursor.End)
        self.edit.insertPlainText(m)

        if self.color:
            self.edit.setTextColor(tc)
        
        if self.out:
            self.out.write(m)
#######################

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

    def OpenFile(self,line_edit):
        FileDialog = QFileDialog(self)
        filepath = FileDialog.getOpenFileName(self,"Open File")
        FileDialog.setWindowTitle("Open File")
        line_edit.setText(filepath)

    def pearson(self):
        eng_path = self.txt13.text().replace("\\",'/')
        rmsd_path = self.txt14.text().replace("\\",'/')
        eng = []
        rmsd = []
        f = open(eng_path)
        for line in f.readlines():
            temp = line.strip().split()
            eng.append(float(temp[0]))
        f.close()
        f = open(rmsd_path)
        for line in f.readlines():
            temp = line.strip().split()
            rmsd.append(float(temp[0]))
        f.close()
        tools.plot_scatter(eng,rmsd)
        self.ShowPics('pearson')


    def ShowPics(self,name):
        FileDialog = QDialog(self)
        pic_txt = QLabel(self)
        la = QHBoxLayout()
        pic_txt.setPixmap(QPixmap("../cache/pics/" + name + '.png'))
        la.addWidget(pic_txt)
        FileDialog.setLayout(la)
        FileDialog.setGeometry(300,300,640,480)
        FileDialog.setWindowTitle("Amino Acid Statistics")
        FileDialog.exec_()

    def PrimaryExtract(self):
        pdb_path = self.txt21.text().replace("\\",'/')
        f = open(pdb_path)
        primary_seq, amino_dict = tools.primary_extract(f)
        tools.plot_bar(amino_dict)
        print("Primary Structure Sequence:")
        print primary_seq
        print("Amino Acid Statistics:")
        print amino_dict
        self.ShowPics('statistic')

    def RMSD(self):
        ori_path = self.txt19.text().replace("\\",'/')
        tar_path = self.txt20.text().replace("\\",'/')

        print("Calculating the RMSD:")
        print(ori_path)
        print(tar_path)

        f = open("../cache/pymol/rmsd.pml",'w')
        f.write("load " + ori_path + ', ori' + '\n')
        f.write("load " + tar_path + ', tar' + '\n')
        f.write("align ori, tar")
        f.close()
        os.system("PyMOLWin.exe ../cache/pymol/rmsd.pml")

    def Show_3D(self):
        pdb_path = self.txt15.text().replace("\\",'/')
        print("Showing 3D Visual:")
        print(pdb_path)

        f = open("../cache/pymol/3d.pml",'w')
        f.write("load " + pdb_path)
        f.close()
        os.system("PyMOLWin.exe ../cache/pymol/3d.pml")

    def match_segment(self):
        ori_path = self.txt17.text().replace("\\",'/')
        tar_path = self.txt18.text().replace("\\",'/')
        print("Find match segment:")
        print(ori_path)
        print(tar_path)

        f = open("../cache/pymol/match.pml",'w')
        f.write("load " + ori_path + ", na" + '\n')
        f.write("load " + tar_path + ", dcy" + '\n')
        f.write("super na, dcy, object=aln" + '\n')
        f.write("python" + '\n')
        f.write("pairs = cmd.get_raw_alignment('aln')" + '\n')
        f.write("print pairs[0], pairs[-1]" + '\n')
        f.write("print pairs[0][1][1], pairs[-1][1][1], pairs[0][1][1]-pairs[-1][1][1]" + '\n')
        line = '"sel1","na and index %d-%d"%(pairs[0][0][1],pairs[-1][0][1])'
        f.write("cmd.select(" + line + ')' + '\n')
        f.write('cmd.save("../cache/files/match_segment.pdb","sel1")' + '\n')
        f.write("python end")
        f.close()
        pwd = os.getcwd()
        father_path=os.path.abspath(os.path.dirname(pwd)+os.path.sep+".").replace("\\", "/")
        results_path = father_path + '/cache/files/match_segment.pdb'
        print("See results at:" + results_path)
        os.system("PyMOLWin.exe -c ../cache/pymol/match.pml")

    

    def nepref_cal(self):
        if(self.checkbox1.isChecked()):
            pwd = os.getcwd()
            father_path=os.path.abspath(os.path.dirname(pwd)+os.path.sep+".").replace("\\", "/")
            matrix_path = father_path + '/backend/Cutoff/6.npy'
            pdb_path = self.txt1.text().replace("\\",'/')
            print("*****Job Information*****")
            print("Algorithm: NEPRE-F")
            print("Matrix: " + matrix_path)
            print("Cutoff: 6Å")
            print("*************************")
            QApplication.processEvents()
            tools.nepre_f(pdb_path,matrix_path,6)

        else:
            matrix_path = self.txt4.text().replace("\\",'/')
            pdb_path = self.txt1.text().replace("\\",'/')
            cutoff = int(self.combobox0.currentText()[:-2])
            print("*****Job Information*****")
            print("Algorithm: NEPRE-F")
            print("Matrix: " + matrix_path)
            print("Cutoff: " + str(cutoff) + 'Å')
            print("*************************")
            QApplication.processEvents()
            tools.nepre_f(pdb_path,matrix_path,cutoff)
            print("Job Finish!")
            
    def neprer_cal(self):
        if(self.checkbox2.isChecked()):
            pwd = os.getcwd()
            father_path=os.path.abspath(os.path.dirname(pwd)+os.path.sep+".").replace("\\", "/")
            matrix_path = father_path + '/backend/Radius/radius.npy'
            radius_path = father_path + '/backend/Radius/mean_radius.txt'
            pdb_path = self.txt2.text().replace("\\",'/')
            print("*****Job Information*****")
            print("Algorithm: NEPRE-R")
            print("Matrix: " + matrix_path)
            print("*************************")
            QApplication.processEvents()
            tools.nepre_r(pdb_path,matrix_path,radius_path)
            print("Job Finish!")
        else:
            matrix_path = self.txt5.text().replace("\\",'/')
            radius_path = self.txt6.text().replace("\\",'/')
            pdb_path = self.txt2.text().replace("\\",'/')
            print("*****Job Information*****")
            print("Algorithm: NEPRE-R")
            print("Matrix: " + matrix_path)
            print("RadiusFile:" + radius_path)
            print("*************************")
            QApplication.processEvents()
            tools.nepre_r(pdb_path,matrix_path,radius_path)
            print("Job Finish!")


    def initUI(self):
        
        ################## NEPRE-PART ######################
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

        # tools part
        layout03 = QVBoxLayout()

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
        
        

        ####################################################

        ##################### TOOLS PART ###################
        self.button15 = QPushButton("Select Energy")
        self.button15.setStyleSheet("background:linen")
        self.button16 = QPushButton("Select RMSD")
        self.button16.setStyleSheet("background:linen")
        self.button17 = QPushButton("Select PDB")
        self.button17.setStyleSheet("background:wheat")
        self.button18 = QPushButton("Select PDB")
        self.button18.setStyleSheet("background:wheat")
        self.button18.setStyleSheet("background:linen")
        self.button19 = QPushButton("Select Origin")
        self.button19.setStyleSheet("background:wheat")
        self.button20 = QPushButton("Select Target")
        self.button20.setStyleSheet("background:wheat")
        self.button21 = QPushButton("Select Origin")
        self.button21.setStyleSheet("background:linen")
        self.button22 = QPushButton("Select Target")
        self.button22.setStyleSheet("background:linen")
        self.button23 = QPushButton("Select PDB")
        self.button24 = QPushButton("Go")
        self.button24.setStyleSheet("background:lightsteelblue")
        self.button25 = QPushButton("Go")
        self.button25.setStyleSheet("background:lightsteelblue")
        self.button26 = QPushButton("Go")
        self.button26.setStyleSheet("background:lightsteelblue")
        self.button27 = QPushButton("Go")
        self.button27.setStyleSheet("background:lightsteelblue")
        self.button28 = QPushButton("Go")
        self.button28.setStyleSheet("background:lightsteelblue")
        self.button29 = QPushButton("Go")
        self.button29.setStyleSheet("background:lightsteelblue")
        # console
        self.txt12 = QTextEdit()
        # Line Edit
        self.txt13 = QLineEdit()
        self.txt13.setStyleSheet("background:linen")
        self.txt14 = QLineEdit()
        self.txt14.setStyleSheet("background:linen")
        self.txt15 = QLineEdit()
        self.txt15.setStyleSheet("background:wheat")
        self.txt16 = QLineEdit()
        self.txt16.setStyleSheet("background:linen")
        self.txt17 = QLineEdit()
        self.txt17.setStyleSheet("background:wheat")
        self.txt18 = QLineEdit()
        self.txt18.setStyleSheet("background:wheat")
        self.txt19 = QLineEdit()
        self.txt19.setStyleSheet("background:linen")
        self.txt20 = QLineEdit()
        self.txt20.setStyleSheet("background:linen")
        self.txt21 = QLineEdit()
        self.txt21.setStyleSheet("background:wheat")

        # group box
        self.groupbox7 = QGroupBox("Basic Plot")
        self.groupbox8 = QGroupBox("3D Visual")
        self.groupbox9 = QGroupBox("Pair Distribution")
        self.groupbox10 = QGroupBox("Homology Analysis")
        self.groupbox11 = QGroupBox("RMSD")
        self.groupbox12 = QGroupBox("Primary Extracted")
        self.groupbox13 = QGroupBox("Console")

        
        # First class layout
        layout_dp = QVBoxLayout()
        # Second class layout
        layout01_dp = QHBoxLayout()
        layout02_dp = QHBoxLayout()
        # Third class layout
        layout001_dp = QVBoxLayout()
        layout002_dp = QVBoxLayout()
        # Forth class layout
        layout0001_dp = QVBoxLayout()
        layout0002_dp = QVBoxLayout()
        layout0003_dp = QVBoxLayout()
        layout0004_dp = QVBoxLayout()
        layout0005_dp = QVBoxLayout()
        layout0006_dp = QVBoxLayout()
        layout0007_dp = QVBoxLayout()
        # Fifth class layout
        layout00001_dp = QHBoxLayout()
        layout00002_dp = QHBoxLayout()
        layout00003_dp = QHBoxLayout()
        layout00004_dp = QHBoxLayout()
        layout00005_dp = QHBoxLayout()
        layout00006_dp = QHBoxLayout()
        layout00007_dp = QHBoxLayout()
        layout00008_dp = QHBoxLayout()
        layout00009_dp = QHBoxLayout()
        layout000010_dp = QHBoxLayout()
        layout000011_dp = QHBoxLayout()
        layout000012_dp = QHBoxLayout()
        layout000013_dp = QHBoxLayout()
        layout000014_dp = QHBoxLayout()
        layout000015_dp = QHBoxLayout()


        # layout add component
        ######### Basic Plot ##########
        layout00001_dp.addWidget(self.button15)
        layout00001_dp.addWidget(self.txt13)
        layout00002_dp.addWidget(self.button16)
        layout00002_dp.addWidget(self.txt14)
        layout000010_dp.addSpacing(350)
        layout000010_dp.addWidget(self.button24)
        layout0001_dp.addLayout(layout00001_dp)
        layout0001_dp.addLayout(layout00002_dp)
        layout0001_dp.addLayout(layout000010_dp)
        self.groupbox7.setLayout(layout0001_dp)
        self.groupbox7.setStyleSheet("background:lavender")
        #layout001_dp.addWidget(self.groupbox7)

        ######### 3D Visual ############
        layout00003_dp.addWidget(self.button17)
        layout00003_dp.addWidget(self.txt15)
        layout000011_dp.addSpacing(350)
        layout000011_dp.addWidget(self.button25)
        layout0003_dp.addLayout(layout00003_dp)
        layout0003_dp.addLayout(layout000011_dp)
        self.groupbox8.setLayout(layout0003_dp)
        self.groupbox8.setStyleSheet("background:lavenderblush")
        #layout001_dp.addWidget(self.groupbox8)

        ######### Pair Distribution ########
        layout00004_dp.addWidget(self.button18)
        layout00004_dp.addWidget(self.txt16)
        layout000012_dp.addSpacing(350)
        layout000012_dp.addWidget(self.button26)
        layout0004_dp.addLayout(layout00004_dp)
        layout0004_dp.addLayout(layout000012_dp)
        self.groupbox9.setLayout(layout0004_dp)
        self.groupbox9.setStyleSheet("background:lavender")
        #layout001_dp.addWidget(self.groupbox9)

        ########## Homology analysis ###########
        layout00005_dp.addWidget(self.button19)
        layout00005_dp.addWidget(self.txt17)
        layout00006_dp.addWidget(self.button20)
        layout00006_dp.addWidget(self.txt18)
        layout000013_dp.addSpacing(350)
        layout000013_dp.addWidget(self.button27)
        layout0005_dp.addLayout(layout00005_dp)
        layout0005_dp.addLayout(layout00006_dp)
        layout0005_dp.addLayout(layout000013_dp)
        self.groupbox10.setLayout(layout0005_dp)
        self.groupbox10.setStyleSheet("background:lavenderblush")
        #layout002_dp.addWidget(self.groupbox10)

        ######### RMSD #########
        layout00007_dp.addWidget(self.button21)
        layout00007_dp.addWidget(self.txt19)
        layout00008_dp.addWidget(self.button22)
        layout00008_dp.addWidget(self.txt20)
        layout000014_dp.addSpacing(350)
        layout000014_dp.addWidget(self.button28)
        layout0006_dp.addLayout(layout00007_dp)
        layout0006_dp.addLayout(layout00008_dp)
        layout0006_dp.addLayout(layout000014_dp)
        self.groupbox11.setLayout(layout0006_dp)
        self.groupbox11.setStyleSheet("background:lavender")
        #layout002_dp.addWidget(self.groupbox11)

        ######### Primary Structure Extract #######
        layout00009_dp.addWidget(self.button23)
        layout00009_dp.addWidget(self.txt21)
        layout000015_dp.addSpacing(350)
        layout000015_dp.addWidget(self.button29)
        layout0007_dp.addLayout(layout00009_dp)
        layout0007_dp.addLayout(layout000015_dp)
        self.groupbox12.setLayout(layout0007_dp)
        self.groupbox12.setStyleSheet("background:lavenderblush")
        #layout002_dp.addWidget(self.groupbox12)



        #################### TOOLS PART FINISH #############
        layout03.addWidget(self.groupbox7)
        layout03.addWidget(self.groupbox8)
        layout03.addWidget(self.groupbox9)
        layout03.addWidget(self.groupbox10)
        layout03.addWidget(self.groupbox11)
        layout03.addWidget(self.groupbox12)
        layout.addLayout(layout01)
        layout.addLayout(layout02)
        layout.addLayout(layout03)
        self.setLayout(layout)

        #################### System Standard Output Redirect ################# 
        sys.stdout = OutLog(self.txt3, sys.stdout)
        sys.stderr = OutLog(self.txt3, sys.stderr, QColor(255,0,0))

        # singal and slot
        # Nepre-F
        self.button1.clicked.connect(lambda:self.OpenFile(self.txt1))
        self.button3.clicked.connect(lambda:self.OpenFile(self.txt4))
        self.button4.clicked.connect(self.nepref_cal)
        # Nepre-R
        self.button2.clicked.connect(lambda:self.OpenFile(self.txt2))
        self.button5.clicked.connect(lambda:self.OpenFile(self.txt5))
        self.button6.clicked.connect(lambda:self.OpenFile(self.txt6))
        self.button7.clicked.connect(self.neprer_cal)

        # primary extract
        self.button23.clicked.connect(lambda:self.OpenFile(self.txt21))
        self.button29.clicked.connect(self.PrimaryExtract)

        # rmsd
        self.button21.clicked.connect(lambda:self.OpenFile(self.txt19))
        self.button22.clicked.connect(lambda:self.OpenFile(self.txt20))
        self.button28.clicked.connect(self.RMSD)

        # 3D View
        self.button17.clicked.connect(lambda:self.OpenFile(self.txt15))
        self.button25.clicked.connect(self.Show_3D)

        # segment extract
        self.button19.clicked.connect(lambda:self.OpenFile(self.txt17))
        self.button20.clicked.connect(lambda:self.OpenFile(self.txt18))
        self.button27.clicked.connect(self.match_segment)

        # basic plot
        self.button15.clicked.connect(lambda:self.OpenFile(self.txt13))
        self.button16.clicked.connect(lambda:self.OpenFile(self.txt14))
        self.button24.clicked.connect(self.pearson)
        

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
        self.resize(1300,900)
        self.nepre = nepreUI()
        #self.dataprocess = DataProcess()
        self.info = InstucAndInfo()
        self.addTab(self.nepre,"Structure Assessment")
        #self.addTab(self.dataprocess, "Data Analyze")
        self.addTab(self.info, "Instruction && Information")
        title = "NEPRE--Scoring Function Based on Neighbourhood Preference Statistics"
        self.setWindowTitle(title)
        



if __name__ == "__main__":
    app = QApplication(sys.argv)
    t = TabWidget()
    t.show()
    app.exec_()