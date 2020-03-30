'''
Created on Feb 10, 2020

@author: griff
'''
import math 
class Particle:
    def __init__(self, pdg, pt, eta, phi, m):            
        self.pdg = pdg
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.m = m
    def getpt(self):
        return self.pt
    def geteta(self):
        return self.eta
    def getphi(self):
        return self.phi
    def Px(self):
        return self.pt*math.cosh(self.phi)
    def Py(self):
        return self.pt*math.sin(self.phi)
    def Pz(self):
        return self.pt*math.sinh(self.eta)
    def E(self): 
        return math.sqrt(self.m**2+self.pt**2*math.cosh(self.eta)**2)
    def Mag(self):
        x = self.Px()
        y = self.Py()
        z =self.Pz()
        E = self.E()
        return math.sqrt((E**2)+(x**2)+(y**2)+(z**2))
    def Angle_3D(self):
        return 
    def DP(self,m2):
        return self.m*self.m2*math.cos(self.Angle_3D)
    def inv_M(self, pt2, eta2, phi2):
        return math.sqrt(2*self.pt*self.pt2*(math.cosh(self.eta2-self.eta)-math.cos(self.phi2-self.phi)))
        
def muon():        
    fo = open('muons.txt','r')
    f = open('outfile.txt','w')
    count = 0
    line = [fo.readline(),fo.readline(),fo.readline(),fo.readline(),fo.readline()]
    # while count < 5:                #for few events
    while line != ['','','','','']: #for all events
        muon1 = line[2].split(' ')
        muon2 = line[3].split(' ')
        pt1,eta1,phi1,m1 = float(muon1[2]), float(muon1[3]), float(muon1[4]), float(muon1[5])
        pt2,eta2,phi2,m2 = float(muon2[2]), float(muon2[3]), float(muon2[4]), float(muon2[5])
        x = Particle(13,pt1,eta1,phi1,m1)
        y = Particle(13,pt2,eta2,phi2,m2)
        line = [fo.readline(),fo.readline(),fo.readline(),fo.readline(),fo.readline()]
        count += 1
        data1 = (x.inv_M(y.getpt(),y.geteta(),y.getphi()))
        data2 = (y.inv_M(x.getpt(),x.geteta(),x.getphi()))
        f.write(str(data1) + '\n')
        f.write(str(data2) + '\n')
    fo.close()
    f.close()

    

