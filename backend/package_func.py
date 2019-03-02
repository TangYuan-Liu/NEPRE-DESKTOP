import sys
sys.path.append("..")
import os
import Cutoff.Nepre_F as nf
import Radius.Nepre_R as nr 
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from PyQt4.QtWebKit import *

############## NEPRE-F ##################
def nepre_f(pdb_path,matrix_path,cutoff):
    if not os.path.exists(pdb_path):
        print("File or path not exists")
        return
    if os.path.isfile(pdb_path):
        matrix = nf.load_EnergyMatrix(matrix_path)
        f = open(pdb_path)
        print "Analyzing: ", pdb_path
        QApplication.processEvents()
        eng = nf.calculate_Energy(f,matrix,cutoff)
        print("*****Results*****")
        print(eng)
        print("*****************")
    
    else:
        correct_list = []
        error_list = []
        for f in os.listdir(pdb_path):
            if(f[-3:] == "pdb"):
                correct_list.append(f)
            else:
                error_list.append(f)
        if(len(correct_list) == 0):
            print("Error: All of the file format uncorrect!")
            return
        if(len(error_list) != 0):
            print("Warning: Uncorrect file detected. But system will continue calculate the correct files.")
            print("Uncorrect List:")
            print(error_list)

        matrix = nf.load_EnergyMatrix(matrix_path)
        eng_list = []
        if(pdb_path[-1] != '/'):
            pdb_path += '/'
        for decoy in correct_list:
            print "Analyzing: ",decoy
            QApplication.processEvents()
            file_path = pdb_path + decoy
            f = open(file_path)
            eng_list.append(nf.calculate_Energy(f,matrix,cutoff))
            f.close()
        
        print("*****Results*****")
        print("FileName" + ' ' + "Energy")
        for i in range(len(correct_list)):
            print correct_list[i] + ' ', eng_list[i] 
        print("*****************")
        

############### NEPRE-R ###################
def nepre_r(pdb_path,matrix_path,radius_path):
    if not os.path.exists(pdb_path):
        print("File or path not exists")
        return
    if os.path.isfile(pdb_path):
        matrix = nr.load_EnergyMatrix(matrix_path)
        radius_dict = nr.LoadRadius(radius_path)
        f = open(pdb_path)
        print "Analyzing: ", pdb_path
        QApplication.processEvents()
        eng = nr.calculate_Energy(f,matrix,radius_dict)
        print("*****Results*****")
        print(eng)
        print("*****************")
    
    else:
        correct_list = []
        error_list = []
        for f in os.listdir(pdb_path):
            if(f[-3:] == "pdb"):
                correct_list.append(f)
            else:
                error_list.append(f)
        if(len(correct_list) == 0):
            print("Error: All of the file format uncorrect!")
            return
        if(len(error_list) != 0):
            print("Warning: Uncorrect file detected. But system will continue calculate the correct files.")
            print("Uncorrect List:")
            print(error_list)
        
        matrix = nr.load_EnergyMatrix(matrix_path)
        radius_dict = nr.LoadRadius(radius_path)
        eng_list = []

        if(pdb_path[-1] != '/'):
            pdb_path += '/'
        for decoy in correct_list:
            print "Analyzing: ",decoy
            QApplication.processEvents()
            file_path = pdb_path + decoy
            f = open(file_path)
            eng_list.append(nr.calculate_Energy(f,matrix,radius_dict))
            f.close()
    
        print("*****Results*****")
        print("FileName" + ' ' + "Energy")
        for i in range(len(correct_list)):
            print correct_list[i] + ' ', eng_list[i] 
        print("*****************")
            



