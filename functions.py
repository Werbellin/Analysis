#!/usr/bin/env python
# Based on P.Onyisi example
import sys,time, os,string
import math
import numpy as np
#if len(sys.argv)<2:
#    print "Usage: python test.py file.root"
#    print "Exit.."
#    sys.exit()

import logging

LEVELS = {'debug': logging.DEBUG,
          'info': logging.INFO,
          'warning': logging.WARNING,
          'error': logging.ERROR,
          'critical': logging.CRITICAL}



from ROOT import TCanvas, TH1D, gSystem, TFile, TTree, TLorentzVector, TChain
TH1D.SetDefaultSumw2()
gSystem.Load('libDelphes.so')

from itertools import tee, islice, chain, izip

def update_progress(progress):
    barLength = 50 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


def previous_and_next(some_iterable):
    prevs, items, nexts = tee(some_iterable, 3)
    prevs = chain([None], prevs)
    nexts = chain(islice(nexts, 1, None), [None])
    return izip(prevs, items, nexts)

class Particle :
    def __init__(self, lorentzvector, charge) :
        self.PT  = lorentzvector.Pt()
        self.Eta = lorentzvector.Eta()
        self.Phi = lorentzvector.Phi()
        self.E   = lorentzvector.E()
        self.Charge = charge
    def __repr__(self) :
        return  "(%f, %f, %f, %f)"% (self.PT, self.Eta, self.Phi, self.E)
    def P4(self) :
        temp = TLorentzVector(0,0,0,0)
        temp.SetPtEtaPhiE(self.PT,self.Eta, self.Phi, self.E)
        return temp

def ExtractObjectsFromGenRecord(event) :
    temp = TLorentzVector(0,0,0,0)
    GenMuons = []
    GenElectrons = []
    for object in event.data.Particle :
        if object.PT <> 0.0 : # need this to avoid ROOT warning on pseudo-rapidity

            temp.SetPtEtaPhiE(object.PT, object.Eta, object.Phi, object.E)
            p = Particle(temp, object.Charge)
            if object.Status == 1 :
                if abs(object.PID) == 11 :
                    GenElectrons.append(p)
                if abs(object.PID) == 13 :
                    GenMuons.append(p)

    event.GenMuon = GenMuons
    event.GenElectron = GenElectrons

    GenJets = []

    for jet in event.data.GenJet :
        GenJets.append(jet)

    event.GenJet = GenJets


    for muon in GenMuons :
        logging.debug(muon)
    for el in GenElectrons :
        logging.debug(el)
    
   # numberOfObjectsInEvent = len(event.data.Particle.PID)
    #
    #for pnumber in xrange(numberOfObjectsInEvent) :
    #    if event.Particle.PID[pnumber] == 11 :
     #       print "Found electron!"



def dR(object1, object2) :
    dR = math.sqrt( (object1.Eta() - object2.Eta())**2  + (object1.Phi() - object2.Phi())**2)
    return dR;

def M(object1, object2) :
    M = math.sqrt(2*object1.PT*object2.PT*(math.cosh(object1.Eta - object1.Eta) - math.cos(object1.Phi - object2.Phi)))
    return M; 

def build_matrix(func,args):
    return func(*args)

def f1(A,B):
    return np.abs(A[:,np.newaxis] - B)


def SeperateByCharge(charge, particle_list):
    result = []
    for p in particle_list:
        p_charge = p.Charge
        if p_charge == charge :
            result.append(p)
    return result;
        #elif mu_charge < 0:
        #   neg_muons.append(mu)
        #else:
            #print "charge error: ", mu.Charge
def GetZMassMatrix(neg_par, pos_par, lower, upper):
    i = len(neg_par)
    j = len(pos_par)
    
    result= np.zeros(shape=(i,j))
    for n in range(0,i):
        for p in range(0,j):
            mass = M(neg_par[n],pos_par[p])
            #massDiffToZ = abs(91.1876 - mass)
            if mass > lower and mass < upper :
                result[n][p] = mass
    return result

def GetScalarPTMatrix(neg_par, pos_par, lower, upper):
    i = len(neg_par)
    j = len(pos_par)
    result= np.zeros(shape=(i,j))
    for n in range(0,i):
        for p in range(0,j):
            mass = M(neg_par[n],pos_par[p])
            if mass > lower and mass < upper :
                result[n][p] = abs(neg_par[n].PT) + abs(pos_par[p].PT)
    return result

def ZCandidatesMuEl(pos_par_type1, neg_par_type1, pos_par_type2, neg_par_type2, lower, upper):
    logging.debug('Enterered ZCandidatesMuEl')
    #print "Mass par type 1: ", M(neg_par_type1[0], pos_par_type1[0])
    #print "Mass par type 2: ", M(neg_par_type2[0], pos_par_type2[0])

    if len(pos_par_type1) >= 1 and len(neg_par_type1) >= 1 and len(pos_par_type2) >=1 and len(neg_par_type2) >= 1 :
        Z1Mass = (neg_par_type1[0].P4() + pos_par_type1[0].P4()).M()
        Z2Mass = (neg_par_type2[0].P4() + pos_par_type2[0].P4()).M()
        if Z1Mass <= upper and Z1Mass >= lower and Z2Mass <= upper and Z2Mass >= lower :
            Z1 = (neg_par_type1[0], pos_par_type1[0])
            Z2 = (neg_par_type2[0], pos_par_type2[0])

            return (Z1, Z2)
        else :
            return ()
    else :
        return ()


def ZCandidates(pos_par_type1, neg_par_type1, pos_par_type2, neg_par_type2, lower, upper):
    logging.debug('Enterered ZCandidates')
    #print "Mass par type 1: ", M(neg_par_type1[0], pos_par_type1[0])
    #print "Mass par type 2: ", M(neg_par_type2[0], pos_par_type2[0])

    Type1MassMatrix = GetZMassMatrix(neg_par_type1, pos_par_type1, lower, upper) 
    Type2MassMatrix = GetZMassMatrix(neg_par_type2, pos_par_type2, lower, upper)


    if len(neg_par_type1) <> 0 and len(pos_par_type1) <> 0 :
        Type1min = np.absolute(91.1876 - Type1MassMatrix).min()
        Type1SPTMatrix = GetScalarPTMatrix(neg_par_type1, pos_par_type1, lower, upper)

    else :
        Type1min = 10e10
        Type1SPTMatrix = np.zeros(1)
    if len(neg_par_type2) <> 0 and len(pos_par_type2) <> 0 :
        Type2min = np.absolute(91.1876 - Type2MassMatrix).min()    
        Type2SPTMatrix = GetScalarPTMatrix(neg_par_type2, pos_par_type2, lower, upper)
    else:
        Type2min = 10e10
        Type2SPTMatrix = np.zeros(1)

    if Type1min < Type2min and np.sum(Type1MassMatrix,dtype=np.int32) > 0:
        #print "found a Z1 candidate of type 1 particles!"
        Z1index = np.unravel_index(Type1MassMatrix.argmin(),Type1MassMatrix.shape)
        Z1 = (neg_par_type1[Z1index[0]] , pos_par_type1[Z1index[1]])
        Type1SPTMatrix.itemset(Z1index,0)
         #Type1DiffMatrix.itemset(Z1index,0)
    elif np.sum(Type2MassMatrix,dtype=np.int32) > 0:
        #print "found a Z1 candidate of type 2 particles!"
        Z1index = np.unravel_index(Type2MassMatrix.argmin(),Type2MassMatrix.shape)
        Z1 = (neg_par_type2[Z1index[0]] , pos_par_type2[Z1index[1]])
        Type2SPTMatrix.itemset(Z1index,0)
    else : 
        return ();

    Type1SPTmax = Type1SPTMatrix.max()
    Type2SPTmax = Type2SPTMatrix.max()
    if Type1SPTmax > Type2SPTmax and Type1SPTmax <> 0 :
        #print "found a Z2 candidate of type 1 particles!"
        Z2index = np.unravel_index(Type1SPTMatrix.argmin(),Type1SPTMatrix.shape)                                                   
        Z2 = (neg_par_type1[Z2index[0]] , pos_par_type1[Z2index[1]])
        return (Z1, Z2)
    elif Type2SPTmax > Type1SPTmax and Type2SPTmax <> 0 :                                                                                   
        #print "found a Z2 candidate of type 2 particles!"                                                                            
        Z2index = np.unravel_index(Type2SPTMatrix.argmin(),Type2SPTMatrix.shape)                                                     
        Z2 = (neg_par_type2[Z2index[0]] , pos_par_type2[Z2index[1]])
        return (Z1, Z2)
    else :
        return ();  


