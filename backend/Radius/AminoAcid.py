import numpy as np

class AminoAcid:
    def __init__(self,name):
        self.center = np.zeros([3,])
        self.SideChainAtomAmount = 0
        self.name = name
        self.xAxis = None
        self.yAxis = None
        self.zAxis = None
        self.CA = None
        self.N = None
        self.SideChain = []
        self.xlist = []
        self.ylist = []
        self.zlist = []
    #Sum all atom coordinate
    def SumCenters(self,x,y,z,atom):
        self.center[0] += x
        self.center[1] += y
        self.center[2] += z
        self.xlist.append(x)
        self.ylist.append(y)
        self.zlist.append(z)
        self.SideChainAtomAmount += 1
        self.SideChain.append(atom)

    #Calculate the center of the amino acid
    def CalculateCenter(self):
        self.center = self.center / self.SideChainAtomAmount
    
    #Calculate distance between two amino acid
    def DistanceBetweenAA(self,center):
        dis = np.sqrt(np.sum((center - self.center)**2))
        return dis

    #Input the CA and N
    def InputCA(self,CA):
        self.CA = CA
    
    def InputN(self,N):
        self.N = N

    #Establish the 3-dimension coordinate
    def EstablishCoordinate(self):
        self.xAxis = (self.N - self.center) / np.sqrt(np.dot((self.N - self.center),(self.N - self.center)))
        self.yAxis = (self.CA - self.center) - (np.dot((self.CA - self.center), self.xAxis))*self.xAxis
        self.yAxis = self.yAxis / np.sqrt((np.dot(self.yAxis,self.yAxis)))
        self.zAxis = np.cross(self.xAxis, self.yAxis)

    #Return the rho, theta, phi
    def ChangeCoordinate(self,center):
        x = np.dot((center - self.center), self.xAxis)
        y = np.dot((center - self.center), self.yAxis)
        z = np.dot((center - self.center), self.zAxis)
        rho = np.sqrt(x**2+y**2+z**2)
        theta = np.arccos(z/rho)
        phi = np.arctan2(y,x)
        return rho,theta,phi

    def Check(self):
        if(self.SideChainAtomAmount == 0):
            return "No Side Chain"
        if(self.N is None):
            return "No Nitrigen"
        if(self.CA is None):
            return "No Alpha Carbon"
        return True
