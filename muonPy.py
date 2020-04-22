'''
Created on Feb 10, 2020

@author: griff
'''
import ROOT as r     #r is now the variable that will import root
r.PyConfig.IgnoreCommandLineOptions = True #root cannot steal cmd-line options
r.gROOT.SetBatch(True)   # no idea what gROOT is
r.gStyle.SetOptStat(False) #histogram
r.gROOT.ProcessLine("gErrorIgnoreLevel=3001;") #ignores root info messages
r.TCanvas.__init__._creates = False #do not create the TCanvas

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
        z = self.Pz()
        E = self.E()
        return math.sqrt((E**2)+(x**2)+(y**2)+(z**2))
    def Angle_3D(self):
        return 
    def DP(self,m2):
        return self.m*self.m2*math.cos(self.Angle_3D)
    def inv_M(self, pt2, eta2, phi2):
        return math.sqrt(2*self.pt*pt2*(math.cosh(eta2-self.eta)-math.cos(phi2-self.phi)))
        


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
	x_tlv = r.TLorentzVector()
        x_tlv.SetPtEtaPhiM(pt1,eta1,phi1,m1) 
        y = Particle(13,pt2,eta2,phi2,m2)
	y_tlv = r.TLorentzVector()
	y_tlv.SetPtEtaPhiM(pt2,eta2,phi2,m2)
	tlv_new = x_tlv + y_tlv
        line = [fo.readline(),fo.readline(),fo.readline(),fo.readline(),fo.readline()]
        count += 1
        data1 = (x.inv_M(y.getpt(),y.geteta(),y.getphi()))
	data2 = tlv_new.M()
        f.write(str(data1) + '\t' + str(data2) + '\n')
    fo.close()
    f.close()

print(muon())

'''def muon_tlv():
    fo = open('muons.txt', 'r')
    count = 0
    line = [fo.readline(),fo.readline(),fo.readline(),fo.readline(),fo.readline()]
    while line != ['','','','','']: #for all events
        muon1 = line[2].split(' ')
        muon2 = line[3].split(' ')
        pt1,eta1,phi1,m1 = float(muon1[2]), float(muon1[3]), float(muon1[4]), float(muon1[5])
        pt2,eta2,phi2,m2 = float(muon2[2]), float(muon2[3]), float(muon2[4]), float(muon2[5])
        x = Particle(13,pt1,eta1,phi1,m1)
        y = Particle(13,pt2,eta2,phi2,m2)
        line = [fo.readline(),fo.readline(),fo.readline(),fo.readline(),fo.readline()]
        count += 1
        data1 = (x.tlv_muon1(y.getpt(),y.geteta(),y.getphi()))
    print(data1)
    fo.close()

print(muon_tlv())'''




'''histogram work

import numbpy as np
import matplotlib.pyplot as plt    (but in ROOT database (TH1F) not numbpy)
tlv_new = tlv_m1 + tlv_m2
x = tlv_new
plt.hist(x, bins = 100)
plt.xlabel('TLV Inv_M')
plt.ylabel('Probability')
plt.show()'''
