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


#in order to use local python module. need to tell python where to look
WORKDIR = '/home/glbradle/workarea/simpleAna'
import sys 
sys.path.append(WORKDIR)

#Use Anyes pyUtils 
import pyUtils.plots.plot_utils as pu

import math 


class Particle(r.TLorentzVector):
    def __init__(self, pdg, pt, eta, phi, m): 
        r.TLorentzVector.__init__(self) #initializes derived class
        self.pdg = int(pdg) 
        self.pt = float(pt)  #float = 
        self.eta = float(eta)
        self.phi = float(phi)
        self.m = float(m)
        self.SetPtEtaPhiM(self.pt,self.eta,self.phi,self.m)
    def getpt(self):
        return self.pt
    def geteta(self):
        return self.eta
    def getphi(self):
        return self.phi
    def Px(self):
        return self.pt*math.cos(self.phi)
    def Py(self):
        return self.pt*math.sin(self.phi)
    def Pz(self):
        return self.pt*math.sinh(self.eta)
    def E(self): 
        return math.sqrt(pow(self.m,2)+pow(self.pt,2)*pow(math.cosh(self.eta),2))
    def Mag(self):
        x = self.Px()
        y = self.Py()
        z = self.Pz()
        E = self.E()
        return math.sqrt((E**2)+(x**2)+(y**2)+(z**2))
    def Angle_3D(self,p):
        return math.acos(self.dotProduct(p)/(self.Mag()*p.Mag()))
    def dotProduct(self,p):
        return self.Px()*p.Px()+self.Py()*p.Py()+self.Pz()*p.Pz()
    def DP(self,m2):
        return self.m*self.m2*math.cos(self.Angle_3D)
    def inv_M(self, pt2, eta2, phi2):
        return math.sqrt((2*self.pt*pt2)*(math.cosh(self.eta-eta2)-math.cos(self.phi-phi2)))
    def inv_M2(self,p):
        return math.sqrt(pow(self.m,2)+pow(p.m,2)+2*(self.E()*p.E()-self.dotProduct(p)))


def readFile():
    fo = open('muons.txt','r')
    f = open('outfile1.txt','w')
    count = 0
    events = []
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
        tlv_new = x_tlv + y_tlv #uses TLV operator for + bc my particle class has no + operator 
        newP = Particle(99,tlv_new.Pt(), tlv_new.Eta(), tlv_new.Phi(), tlv_new.M())
        line = [fo.readline(),fo.readline(),fo.readline(),fo.readline(),fo.readline()]
        count += 1
        events.append(event(0,0,x,y,newP)) #replace 0 with reading run and event
        data1 = (x.inv_M(y.getpt(),y.geteta(),y.getphi()))
        data2 = x.inv_M2(y)
        data3 = newP.M()
        f.write(str(data1) + '\t' + str(data2) + '\t' + str(data3) + '\n')  
       # print (linearray)
    fo.close()
    f.close()
    return events

class event:
    def __init__(self,run = None, evt = None, m1 = None, m2 = None, newP = None):
        self.run = int(run)
        self.evt = int(evt)
        self.mu1 = m1
        self.mu2 = m2
        self.p = newP
#print function, def function    __str(self)__

# =============================================================================
#  Running part of the script   
# =============================================================================
def main(): 
    events = readFile() #loop over events and fill histo #returns event list

    #Set the histogram style to be of the ATLAS style
    style = pu.atlasStyle()
    r.gROOT.SetStyle("ATLAS")
    r.gROOT.ForceStyle()

    #Book a histogram
    h_invM = pu.th1f("newP_invM","newP M",500,0,1000,'M [GeV]')
    #fill 100 dummy entries at the same value of 200 GeV.
    for i in range(0,100):
        h_invM.Fill(newP.M())  
    
    '''
    Ususally, one woudl loop over the events and for each event
    fill the histogram like
    h.Fill(newP.M())  #This is more what you would one to do.
    '''

    #Dump the histos to a ROOT file
    f = r.TFile('output.root','RECREATE')
    h_invM.Write();
    f.Close()





#print(muon())






'''histogram work

import numbpy as np
import matplotlib.pyplot as plt    (but in ROOT database (TH1F) not numbpy)
tlv_new = tlv_m1 + tlv_m2
x = tlv_new
plt.hist(x, bins = 100)
plt.xlabel('TLV Inv_M')
plt.ylabel('Probability')
plt.show()'''

#main function executable

if __name__=="__main__" : 
    main()


#now must run from Simple Ana b/c that is where muons.txt is
