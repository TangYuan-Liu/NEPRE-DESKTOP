import sys
sys.path.append("..")
import os
import Cutoff.Nepre_F as nf
import Radius.Nepre_R as nr 
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from PyQt4.QtWebKit import *
import matplotlib.pyplot as plt

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
        print("*****Prediction*****")
        print("The prediction is:")
        print correct_list[eng_list.index(min(eng_list))]
        print("********************")
        

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


def plot_bar(bar_dict):
    #data = sorted(bar_dict.items(), lambda x, y: cmp(x[1], y[1]))
    data = [(k,bar_dict[k]) for k in sorted(bar_dict.keys())] 
    keys = []
    values = []
    for i in range(len(data)):
        keys.append(data[i][0])
        values.append(data[i][1])

    plt.rcParams['savefig.dpi'] = 90
    plt.rcParams['figure.figsize'] = [12.0,8.0]
    plt.bar(keys,values,width=0.5)
    plt.xticks(keys,rotation=45,fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel("Type",fontsize=18)
    plt.ylabel("Amount",fontsize=18)
     
    plt.savefig("../cache/pics/statistic.png")
    plt.close()

def plot_scatter(eng,rmsd):
    plt.rcParams['savefig.dpi'] = 90
    plt.rcParams['figure.figsize'] = [12.0,8.0]
    plt.scatter(rmsd,eng,color='r',marker='.')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel("RMSD",fontsize=18)
    plt.ylabel("Energy",fontsize=18)
    plt.savefig("../cache/pics/pearson.png")
    plt.close()



def primary_extract(f):
    primary = []
    current_num = None
    amino_dict = {}

    for line in f.readlines():
        temp = line.strip().split()
        if(temp[0] != "ATOM"):
            continue
        line = nf.extract_Data(line)
        residue_num = line[6]
        residue_name = line[4]
        if(current_num is None):
            primary.append(residue_name)
            current_num = residue_num
        else:
            if(current_num != residue_num):
                primary.append(residue_name)
                current_num = residue_num

    for d in primary:
        if(d not in amino_dict):
            amino_dict[d] = 1
        else:
            amino_dict[d] += 1
    
    return primary,amino_dict
            



