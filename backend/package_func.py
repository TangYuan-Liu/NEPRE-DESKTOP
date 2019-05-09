import sys
sys.path.append("..")
import os
import Cutoff.Nepre_F as nf
import Radius.Nepre_R as nr 
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from PyQt4.QtWebKit import *
import matplotlib.pyplot as plt
import Cutoff.AminoAcid as AA
import numpy as np

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


def bilinear_interpolation(matrix,emp=0):
    theta = np.shape(matrix)[0]
    phi = np.shape(matrix)[1]
    res = np.zeros([theta,phi])

    for i in range(theta):
        for j in range(phi):
            if(matrix[i][j] == 0):
                val = 0
                neighbor = []
                neighbor.append(matrix[i-1][j])
                neighbor.append(matrix[(i+1)%theta][j])
                neighbor.append(matrix[i][j-1])
                neighbor.append(matrix[i][(j+1)%phi])
                while(0 in neighbor):
                    neighbor.remove(0)
                if(neighbor):
                    res[i][j] = int(float(sum(neighbor)/len(neighbor)))
                else:
                    res[i][j] = emp
            else:
                res[i][j] = matrix[i][j]

    return res

def extract_Data(line):
    """
    This part will extracted data from line according to the standard 
    PDB file format(Version 3.3.0, Nov.21, 2012)
    """
    res = []

    line = line.strip()
    #record_name
    res.append(line[0:4].strip(' '))

    #atom_serial
    res.append(line[6:11].strip(' '))

    #atom_name
    res.append(line[12:16].strip(' '))

    #alternate_indicator
    res.append(line[16])

    #residue_name
    res.append(line[17:20].strip(' '))

    #chain_id
    res.append(line[21].strip(' '))

    #residue_num
    res.append(line[22:26].strip(' '))

    #xcor
    res.append(line[30:38].strip(' '))

    #ycor
    res.append(line[38:46].strip(' '))

    #zcor
    res.append(line[46:54].strip(' '))
   
    return res

# scan a single pdb file and return the pair coordinate within cutoff
def single_pdb_scan_neprer(f,radius_dict):
    ccDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
            "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
            "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
            "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}
    csDict={"ALA":0,"VAL":0,"LEU":0,"ILE":0,"PHE":0,\
            "TRP":0,"MET":0,"PRO":0,"GLY":0,"SER":0,\
            "THR":0,"CYS":0,"TYR":0,"ASN":0,"GLN":0,\
            "HIS":0,"LYS":0,"ARG":0,"ASP":0,"GLU":0,}
     
    for amino1 in ccDict.keys():
        for amino2 in ccDict.keys():
            ccDict[amino1][amino2] = 0

    # define some useful paramter
    CurrentAANitrogen = None
    CurrentAACA = None
    Currentresidue_num = None
    CurrentAA = None
    # list of amino acids which have side chain
    UseAA_list = []
    # list of amino acids which are not used
    IgnoreAA_list = []

    # coordinate return list
    res = []    
     
    # scan pdb file line one by one
    for line in f.readlines():        
        if(line[0:4] != "ATOM"):
            continue
        # obtain information
        element_list = extract_Data(line)
        record_name = element_list[0]
        atom_name = element_list[2]
        residue_name = element_list[4]
        alternate_indicator = element_list[3]
        residue_num = element_list[-4]
        xcor = float(element_list[-3])
        ycor = float(element_list[-2])
        zcor = float(element_list[-1])
        
        # ignore hydrogen
        if(atom_name == "H"):
            continue
        # ignore amino acid out of the list
        if(residue_name not in ccDict.keys()):
            continue
        # from here start to scan useful amino acid
        # first amino acid
        if(CurrentAA is None):
            CurrentAA = AA.AminoAcid(residue_name)
            Currentresidue_num = residue_num
            if(atom_name == "N" or atom_name == "CA"):
                if(alternate_indicator == " " or alternate_indicator == "A"):
                    if(atom_name == "N"):
                        CurrentAA.InputN(np.array([xcor,ycor,zcor]))
                    else:
                        CurrentAA.InputCA(np.array([xcor,ycor,zcor]))
                else:
                    continue
                
            if(residue_name == "GLY" or atom_name not in {"N","CA","C","O","O1","02"}):
                if(alternate_indicator == " " or alternate_indicator == "A"):
                    CurrentAA.SumCenters(xcor,ycor,zcor,atom_name)
                else:
                    continue
        
        # current amino acid is not the first
        else:
            #If another amino acid begins
            if(residue_num != Currentresidue_num):
                state = CurrentAA.Check()
                # previous amino acid has no problem
                if(state == True):
                    CurrentAA.CalculateCenter()
                    UseAA_list.append(CurrentAA)
                # previous amino acid has problem
                else:
                    info = [state,Currentresidue_num]
                    IgnoreAA_list.append(info)
                
                CurrentAA = AA.AminoAcid(residue_name)
                Currentresidue_num = residue_num
                if(atom_name == "N" or atom_name == "CA"):
                    if(alternate_indicator == " " or alternate_indicator == "A"):
                        if(atom_name == "N"):
                            CurrentAA.InputN(np.array([xcor,ycor,zcor]))
                        else:
                            CurrentAA.InputCA(np.array([xcor,ycor,zcor]))
                    else:
                        continue
                if(residue_name == "GLY" or atom_name not in {"N","CA","C","O","O1","02"}):
                    if(alternate_indicator == " " or alternate_indicator == "A"):
                        CurrentAA.SumCenters(xcor,ycor,zcor,atom_name)
                    else:
                        continue
            #If still the same amino acid
            else:
                if(atom_name == "N" or atom_name == "CA"):
                    if(alternate_indicator == " " or alternate_indicator == "A"):
                        if(atom_name == "N"):
                            CurrentAA.InputN(np.array([xcor,ycor,zcor]))
                        else:
                            CurrentAA.InputCA(np.array([xcor,ycor,zcor]))
                    else:
                        continue
                if(residue_name == "GLY" or atom_name not in {"N","CA","C","O","O1","02"}):
                    if(alternate_indicator == " " or alternate_indicator == "A"):
                        CurrentAA.SumCenters(xcor,ycor,zcor,atom_name)
                    else:
                        continue
    
    state = CurrentAA.Check()
    if(state == True):
        CurrentAA.CalculateCenter()
        UseAA_list.append(CurrentAA)
        CurrentAA = AA.AminoAcid(residue_name)
        Currentresidue_num = residue_num
    
    else:
        info = [state,Currentresidue_num]
        IgnoreAA_list.append(info)
        CurrentAA = AA.AminoAcid(residue_name)
        Currentresidue_num = residue_num

    for i in range(len(UseAA_list)):
        UseAA_list[i].EstablishCoordinate()
        for j in range(len(UseAA_list)):
            if(i == j):
                continue
            else:
                dis = UseAA_list[i].DistanceBetweenAA(UseAA_list[j].center)
                radiusSum = radius_dict[UseAA_list[i].name] + radius_dict[UseAA_list[j].name]
                if(dis < radiusSum):
                    _,_,_,x,y,z = UseAA_list[i].ChangeCoordinate(UseAA_list[j].center)
                    res.append([x,y,z,UseAA_list[i].name,UseAA_list[j].name])
                    ccDict[UseAA_list[i].name][UseAA_list[j].name] += 1
    
    for i in range(len(UseAA_list)):
        csDict[UseAA_list[i].name] += 1
    return res,ccDict,csDict

def neprer_eng(dataset_path,radius_file):
    ccDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
            "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
            "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
            "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},} 
    csDict={"ALA":0.,"VAL":0.,"LEU":0.,"ILE":0.,"PHE":0.,\
            "TRP":0.,"MET":0.,"PRO":0.,"GLY":0.,"SER":0.,\
            "THR":0.,"CYS":0.,"TYR":0.,"ASN":0.,"GLN":0.,\
            "HIS":0.,"LYS":0.,"ARG":0.,"ASP":0.,"GLU":0.,}

    for amino1 in ccDict.keys():
        for amino2 in ccDict.keys():
            ccDict[amino1][amino2] = 0.

    amino_dict = {}
    f = open(radius_file)
    write_file = open("../cache/temp/coordinate.txt",'w')
    write_file.close()
    for line in f.readlines():
        temp = line.strip().split()
        if(temp[1] not in amino_dict):
            amino_dict[temp[1]] = float(temp[0])
    file_list = []
    for f in os.listdir(dataset_path):
        file_list.append(f)
    print("Dataset loading finished, start to calculate.")
    QApplication.processEvents()
    if(dataset_path[-1] != '/'):
        dataset_path += '/'
    for f in file_list:
        path = dataset_path + f
        df = open(path)
        res,tmp1,tmp2 = single_pdb_scan_neprer(df,amino_dict)
        with open("../cache/temp/coordinate.txt",'a') as write_file:
            for cor in res:
                write_file.write(str(cor[0]) + ' ' + str(cor[1]) + ' ' + str(cor[2]) + ' ' + cor[3] + ' ' + cor[4] + '\n')
        write_file.close()
        for amino1 in ccDict.keys():
            for amino2 in ccDict.keys():
                ccDict[amino1][amino2] += tmp1[amino1][amino2]
        for amino1 in ccDict.keys():
            csDict[amino1] += tmp2[amino1]
    
    pair_sum = 0
    for amino1 in ccDict.keys():
        for amino2 in ccDict.keys():
            pair_sum += ccDict[amino1][amino2]
    
    for amino1 in ccDict.keys():
        for amino2 in ccDict.keys():
            ccDict[amino1][amino2] /= pair_sum
    
    sigle_sum = sum(csDict.values())
    for amino1 in csDict.keys():
        csDict[amino1] /= sigle_sum

    for amino1 in ccDict.keys():
        for amino2 in ccDict.keys():
            ccDict[amino1][amino2] = np.log((ccDict[amino1][amino2])/(csDict[amino1]*csDict[amino2]))

    cor_file_path = "../cache/temp/coordinate.txt"
    print("AminoAcid pair scanning finished, start to generate energymatrix.")
    QApplication.processEvents()
    generate_engmatrix_neprer(cor_file_path,ccDict)
    print("Finished!")
    return True

def generate_engmatrix_neprer(coordinate_file,ccDict):
    cdDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
            "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
            "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
            "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}   
    sort_list = cdDict.keys()
    sort_list.sort()
    
    dualArray = np.zeros((20,20))
    for amino1 in cdDict.keys():
        for amino2 in cdDict.keys():
            cdDict[amino1][amino2] = dualArray
    with open(coordinate_file) as f:
        for line in f:
            line = line.strip().split()
            amino1 = line[3]
            amino2 = line[4]
            cor = np.array([float(line[0]),float(line[1]),float(line[2])])
            rho = sum(cor**2)**0.5
            theta = np.arccos(cor[2]/rho)
            phi = np.arctan2(cor[1],cor[0])
            theta = min(int(theta*20/np.pi),19)
            phi = min(int((phi+np.pi)/(2*np.pi/20)),19)
            cdDict[amino1][amino2][theta][phi] += 1
    
    for amino1 in cdDict.keys():
        for amino2 in cdDict.keys():
            min_data = float('inf')
            print("Analyzing: " + amino1 + '-' + amino2)
            QApplication.processEvents()
            for i in range(20):
                for j in range(20):
                    if cdDict[amino1][amino2][i][j] < min_data:
                        min_data = cdDict[amino1][amino2][i][j]
            
            cdDict[amino1][amino2] = bilinear_interpolation(cdDict[amino1][amino2],min_data)
            cdDict[amino1][amino2] = cdDict[amino1][amino2]/cdDict[amino1][amino2].sum()
            Integral = np.ones((20,20))
            for j in range(20):
                Integral[j] = Integral[j]*(np.cos(j*np.pi/20)-np.cos((j+1)*np.pi/20))*((1./3)*(2*np.pi/20))
            Integral = Integral / Integral.sum()
            cdDict[amino1][amino2] = cdDict[amino1][amino2] / Integral
            cdDict[amino1][amino2] = ccDict[amino1][amino2] - np.log(cdDict[amino1][amino2])
    
    f = open("../cache/temp/latest.npy",'wb')
    for amino1 in sort_list:
        for amino2 in sort_list:
            np.save(f, cdDict[amino1][amino2])
    
    return True
            
    

            
        
        


    
    
    




