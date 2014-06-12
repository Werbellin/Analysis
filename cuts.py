#!/usr/local/bin/python$
# -*- coding: utf-8 -*-$

import logging
from functions import *
from ROOT import TH2D
import numpy as np
from step import Step
class Cut(Step) :
    def __init__(self, step_name) :
        super(Cut, self).__init__(step_name)
        self.NumberEventsPassedCut = 0.

    def Initialize(self, Histos, data_name, cut_flow) :
        super(Cut, self).Initialize(data_name, cut_flow)

class Start(Cut) :
    def __init__(self, step_name) :
        super(Start, self).__init__(step_name)
        self.CutAbbreviation = ""
    def Initialize(self, Histos, data_name, cut_flow) :
        super(Start,self).Initialize(Histos, data_name, cut_flow)

    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        self.NumberEventsPassedCut += 1.0
        return True


class ZeppenfeldVariableCut(Cut) :
    def __init__(self, step_name) :
        super(ZeppenfeldVariableCut, self).__init__(step_name)
        self.CutAbbreviation = "YSTAR"
    def Initialize(self, Histos, data_name, cut_flow) :
        super(ZeppenfeldVariableCut,self).Initialize(Histos, data_name, cut_flow)

    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )

        if abs(event.Y1Star) >= 0.2 :
            self.NumberEventsPassedCut += 1.0
            return True
        else :
            return False

class JetVetoCut(Cut) :
    def __init__(self, step_name) :
        super(JetVetoCut, self).__init__(step_name)
        self.CutAbbreviation = "JVeto"
    def Initialize(self, Histos, data_name, cut_flow) :
        super(JetVetoCut,self).Initialize(Histos, data_name, cut_flow)

    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )

        TJ_max = max(event.TaggingJet1.Eta(), event.TaggingJet2.Eta())
        TJ_min = min(event.TaggingJet1.Eta(), event.TaggingJet2.Eta())
        Tracker_eta_max = 2.5
        Tracker_eta_min = -2.5
        if TJ_max < Tracker_eta_max:
            eta_max = TJ_max
        else :
            eta_max = Tracker_eta_max
        if TJ_min < Tracker_eta_min :
            eta_min = Tracker_eta_min
        else:
            eta_min = TJ_min
        #print eta_min , " max: ", eta_max
 
        if data_type == "SIM" :
            for jet in event.data.Jet :
                if jet <> event.TaggingJet1Pointer and jet <> event.TaggingJet2Pointer :
                    #print "eta: ", jet.Eta, " pt ", jet.PT
                    if jet.PT > 25. and jet.Eta > eta_min and jet.Eta < eta_max :
                        #print "Event failed cut"
                        return False
        self.NumberEventsPassedCut += 1.0
        return True


class TaggingJetY1Y2Cut(Cut) :
    def __init__(self, step_name) :
        super(TaggingJetY1Y2Cut, self).__init__(step_name)
        self.CutAbbreviation = "TJY1Y2"
    def Initialize(self, Histos, data_name, cut_flow) :
        super(TaggingJetY1Y2Cut,self).Initialize(Histos, data_name, cut_flow)

    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )

        if event.TaggingJetY1Y2 <= 0. :
            self.NumberEventsPassedCut += 1.0
            return True
        else :
            return False


class LeptonsBetweenTaggingJetsEtaCut(Cut) :
    def __init__(self, step_name) :
        super(LeptonsBetweenTaggingJetsEtaCut, self).__init__(step_name)
        self.CutAbbreviation = "LbTJEta"

    def Initialize(self, Histos, data_name, cut_flow) :
        super(LeptonsBetweenTaggingJetsEtaCut, self).Initialize(Histos, data_name, cut_flow)

    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        eta_max = max(event.TaggingJet1.Eta(), event.TaggingJet2.Eta())
        eta_min = min(event.TaggingJet1.Eta(), event.TaggingJet2.Eta())

        result = True
        leptonsFromZ = []
        leptonsFromZ.extend(event.Z1Particles)
        leptonsFromZ.extend(event.Z2Particles)
        #print leptonsFromZ
        #print eta_min, "max: ", eta_max
        for lepton in leptonsFromZ :
            if lepton.Eta() <= eta_max and lepton.Eta() >= eta_min :
                result *= True
            else :
                result *= False
        if result == True :
            self.NumberEventsPassedCut += 1.0
            return True
        else :
            return False


class LeptonIsolationCut(Cut) :
    def __init__(self, step_name) :
        super(LeptonIsolationCut, self).__init__(step_name)
        self.CutAbbreviation = "LISO"

    def Initialize(self, Histos, data_name, cut_flow) :
        super(LeptonIsolationCut, self).Initialize(Histos, data_name, cut_flow)

    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        for lepton1 in event.goodLeptons :
            for lepton2 in event.goodLeptons :
                if lepton1 <> lepton2 and dR(lepton1, lepton2) <= 0.3 :
                    event.Cuts[self.name] = False
                    return False

        event.Cuts[self.name] = True
        self.NumberEventsPassedCut += 1.
        return True


class DefineTaggingJetsCut(Cut) :
    def __init__(self, step_name) :
        super(DefineTaggingJetsCut, self).__init__(step_name)
        self.CutAbbreviation = "TJDEF"

    def Initialize(self, Histos, data_name, cut_flow) :
        super(DefineTaggingJetsCut, self).Initialize(Histos, data_name, cut_flow)

    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" + "in mode " + data_type )
        event.Cuts[self.name] = False

        goodJets = []

        if data_type == "SIM" :
            for jet in event.data.Jet :
                if jet.PT > 30. and abs(jet.Eta) < 5.2 :
                    goodJets.append(jet)
        if data_type == "GEN" :
            for jet in event.GenJet :
                 if jet.PT > 30. and abs(jet.Eta) < 5.2 :
                    goodJets.append(jet)
        goodJets.sort(key=lambda x: x.PT, reverse=True)

        if len(goodJets) >= 2 :
            event.Cuts[self.name] = True
            TaggingJet1 = goodJets[0]
            TaggingJet2 = goodJets[1]
            event.TaggingJet1 = TaggingJet1.P4()
            event.TaggingJet1Pointer = TaggingJet1
            event.TaggingJet2 = TaggingJet2.P4()
            event.TaggingJet2Pointer = TaggingJet2
            event.goodJets = goodJets
            self.NumberEventsPassedCut += 1.0

            rapidity_product = TaggingJet1.P4().Rapidity() * TaggingJet2.P4().Rapidity()
            gap = abs(TaggingJet1.P4().Rapidity() - TaggingJet2.P4().Rapidity())

            DPhi = abs(TaggingJet1.Phi - TaggingJet2.Phi)
            EtaSum = TaggingJet1.Eta + TaggingJet2.Eta
            event.TaggingJetEtaSum = EtaSum
            event.TaggingJetDPhi = DPhi
            event.TaggingJetY1Y2 = rapidity_product
            event.TaggingJetRapidityGap = gap
            return True
        return False

class TaggingJetRapidityGapCut(Cut) :
    def __init__(self, step_name) :
        super(TaggingJetRapidityGapCut, self).__init__(step_name)
        self.CutAbbreviation = "TJDY"

    def Initialize(self, Histos, data_name, cut_flow) :
        super(TaggingJetRapidityGapCut, self).Initialize(Histos, data_name, cut_flow)

    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        if event.TaggingJetRapidityGap > 3.0 :
            event.Cuts[self.name] = True
            self.NumberEventsPassedCut += 1.
            return True
        else :
            event.Cuts[self.name] = False
            return False

        return True

class TaggingJetSumEtaCut(Cut) :
    def __init__(self, step_name) :
        super(TaggingJetSumEtaCut, self).__init__(step_name)
        self.CutAbbreviation = "TJSETA"

    def Initialize(self, Histos, data_name, cut_flow) :
        super(TaggingJetSumEtaCut, self).Initialize(Histos, data_name, cut_flow)

    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        etaSum = event.TaggingJet1.Eta() + event.TaggingJet2.Eta()
        if etaSum > -380. :
            event.Cuts[self.name] = True
            self.NumberEventsPassedCut += 1.
            return True
        else :
            event.Cuts[self.name] = False
            return False

        return True

class TrueZCut(Cut) :
    def __init__(self, step_name) :
        super(TrueZCut, self).__init__(step_name)
        self.CutAbbreviation = "TrueZ"

    def Initialize(self, Histos, data_name, cut_flow) :
        super(TrueZCut, self).Initialize(Histos, data_name, cut_flow)

    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        if data_type == "GEN" :
            if len(event.goodMuons) >= 2 and len(event.goodElectrons) >=2 :
                Z1Mass = (event.goodMuons[0].P4() + event.goodMuons[1].P4()).M()
                Z2Mass = (event.goodElectrons[0].P4() + event.goodElectrons[1].P4()).M()
                l = 66.
                u = 116.
                if Z1Mass >= l and Z1Mass <= u and Z2Mass >= l and Z2Mass <= u:
                    self.NumberEventsPassedCut += 1.
                    return True
        return False

class TaggingJetInvariantMassCut(Cut) :
    def __init__(self, step_name) :
        super(TaggingJetInvariantMassCut, self).__init__(step_name)
        self.CutAbbreviation = "TJM"

    def Initialize(self, Histos, data_name, cut_flow) :
        super(TaggingJetInvariantMassCut, self).Initialize(Histos, data_name, cut_flow)

    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        mass = (event.TaggingJet1 + event.TaggingJet2).M()
        if mass > 380. :
            event.Cuts[self.name] = True
            self.NumberEventsPassedCut += 1.
            return True
        else :
            event.Cuts[self.name] = False
            return False

        return True

class LeptonAcceptanceAnalysisCut(Cut) :
    def __init__(self, step_name) :
        super(LeptonAcceptanceAnalysisCut, self).__init__(step_name)
        self.CutAbbreviation = "LAccAna"

    def Initialize(self, Histos, data_name, cut_flow) :
        super(LeptonAcceptanceAnalysisCut, self).Initialize(Histos, data_name, cut_flow)

        eta_bin = 5
        eta_max  = 2.5
        eta_min  = -2.5

        PT_bin = 20
        PT_max  = 200.
        PT_min  = 0.

        PT_ETA_bin = 10
        PT_ETA_max = 200

        self.ElectronEta = TH1D(self.DataName + self.Cutflow + "El_Eta","Number of electrons", eta_bin, eta_min, eta_max)
        self.MuonEta = TH1D(self.DataName + self.Cutflow + "Muon_Eta","Number of muons",eta_bin, eta_min, eta_max)
        self.ElectronPtEta = TH2D(self.DataName + self.Cutflow + "El_Pt_Eta", "Electron efficency", PT_ETA_bin, PT_min, PT_ETA_max, eta_bin, eta_min, eta_max)
        Histos.append(self.ElectronPtEta)
        self.MuonPtEta = TH2D(self.DataName + self.Cutflow + "Muon_Pt_Eta", "Muon efficency", PT_ETA_bin, PT_min, PT_ETA_max, eta_bin, eta_min, eta_max)
        
        Histos.append(self.MuonPtEta)
        Histos.append(self.ElectronEta)
        Histos.append(self.MuonEta)

        self.ElectronPT = TH1D(self.DataName + self.Cutflow + "El_PT","PT of leading lepton", PT_bin, PT_min, PT_max)
        self.MuonPT = TH1D(self.DataName +self.Cutflow +"Muon_PT","PT of subleading lepton",PT_bin, PT_min, PT_max)

        Histos.append(self.ElectronPT)
        Histos.append(self.MuonPT)

        self.LeptonMinPT = TH1D(self.DataName + self.Cutflow + "Lepton_Min_PT","PT of 4th lepton", 30, PT_min, 30.)
        Histos.append(self.LeptonMinPT)

        R_bin = 100
        R_min = 0.
        R_max = 7.
        self.LeptonR = TH1D(self.DataName + self.Cutflow + "-RLiLj","Mass of tagging jets",R_bin, R_min, R_max)
        Histos.append(self.LeptonR)

    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        event.Cuts[self.name] = False
        goodMuons       = []
        goodElectrons   = []

        if data_type == "GEN" :
            ExtractObjectsFromGenRecord(event)
            for mu in event.GenMuon :
                if abs(mu.Eta) < 2.4 and mu.PT > 10. :
                    goodMuons.append(mu)

            for el in event.GenElectron :
                if abs(el.Eta) < 2.5 and el.PT > 4. :
                    goodElectrons.append(el)

        if data_type == "SIM" :
            for mu in event.data.Muon :
                if abs(mu.Eta) < 2.4 and mu.PT > 10. :
                    goodMuons.append(mu)

            for el in event.data.Electron :
                if abs(el.Eta) < 2.5 and el.PT > 4. :
                    goodElectrons.append(el)

        goodLeptons = []
        goodLeptons.extend(goodMuons)
        goodLeptons.extend(goodElectrons)
        goodLeptons.sort(key=lambda x: x.PT, reverse=True)
        goodMuons.sort(key=lambda x: x.PT, reverse=True)
        goodElectrons.sort(key=lambda x: x.PT, reverse=True)

        for mu in goodMuons :
            self.MuonPtEta.Fill(mu.PT, mu.Eta)
            self.MuonEta.Fill(mu.Eta)
            self.MuonPT.Fill(mu.PT)

        for el in goodElectrons :
            self.ElectronPtEta.Fill(el.PT, el.Eta)
            self.ElectronEta.Fill(el.Eta)
            self.ElectronPT.Fill(el.PT)

        for l1 in goodLeptons :
            for l2 in goodLeptons :
                if dR(l1.P4(),l2.P4()) <> 0. :
                    self.LeptonR.Fill(dR(l1.P4(), l2.P4()))

        if len(goodLeptons) >= 4 :
            self.NumberEventsPassedCut += 1
            self.LeptonMinPT.Fill(goodLeptons[3].PT)
        return True



class LeptonDefinitionCut(Cut) :
    def __init__(self, step_name, el_pt = 7., el_eta = 2.5, mu_pt = 7., mu_eta = 2.4) :
        super(LeptonDefinitionCut, self).__init__(step_name)
        self.CutAbbreviation = "LDEF"
        self.ElPtCut    = el_pt
        self.ElEtaCut   = el_eta
        self.MuPtCut    = mu_pt
        self.MuEtaCut   = mu_eta

    def __str__(self) :
        name = self.name + "[pTe>" + str(self.ElPtCut) + ", |ηe|<" + str(self.ElEtaCut) \
                        + ", pTμ>" + str(self.MuPtCut) + ", |ημ|<" + str(self.MuEtaCut) + "]"
        return name

    def Initialize(self, Histos, data_name, cut_flow) :
        super(LeptonDefinitionCut, self).Initialize(Histos, data_name, cut_flow)

    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        event.Cuts[self.name] = False
        goodMuons       = []
        goodElectrons   = []

        if data_type == "GEN" :
            ExtractObjectsFromGenRecord(event)
            for mu in event.GenMuon :
                if abs(mu.Eta) < 2.4 and mu.PT > 7. :
                    goodMuons.append(mu)

            for el in event.GenElectron :
                if abs(el.Eta) < 2.5 and el.PT > 7. :
                    goodElectrons.append(el)

        if data_type == "SIM" :
            for mu in event.data.Muon :
                if abs(mu.Eta) < 2.4 and mu.PT > 7. :
                    goodMuons.append(mu)

            for el in event.data.Electron :
                if abs(el.Eta) < 2.5 and el.PT > 7. :
                    goodElectrons.append(el)

        goodLeptons = []
        goodLeptons.extend(goodMuons)
        goodLeptons.extend(goodElectrons)
        goodLeptons.sort(key=lambda x: x.PT, reverse=True)
        goodMuons.sort(key=lambda x: x.PT, reverse=True)
        goodElectrons.sort(key=lambda x: x.PT, reverse=True)

        if (len(goodMuons) + len(goodElectrons)) >= 4 :
            logging.debug("Found 4 good leptons")
            event.Cuts[self.name] = True
            self.NumberEventsPassedCut += 1
            event.goodMuons = goodMuons
            event.goodElectrons = goodElectrons

            event.goodLeptons = goodLeptons

            return True
        else:
               # for mu in event.data.Muon :
               #     print "Muon pt: ", mu.PT , "eta: ", mu.Eta
               # for el in event.data.Electron :
               #     print "Elec PT: " , el.PT , "eta ", el.Eta
            return False

class ZPairDefinitionCut(Cut) :
    def __init__(self, step_name) :
        super(ZPairDefinitionCut, self).__init__(step_name)
        self.CutAbbreviation = "2ZDEF"

    def Initialize(self, Histos, data_name, cut_flow) :
        super(ZPairDefinitionCut, self).Initialize(Histos, data_name, cut_flow)

    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut")
        event.Cuts[self.name] = False

        pos_muons       = SeperateByCharge(1,event.goodMuons)
        pos_electrons   = SeperateByCharge(1,event.goodElectrons)
        neg_muons       = SeperateByCharge(-1,event.goodMuons)
        neg_electrons   = SeperateByCharge(-1,event.goodElectrons)

        if ((len(neg_muons) >= 2 and len(pos_muons) >= 2) or
            (len(neg_electrons) >= 2 and len(pos_electrons) >= 2) or
            (len(neg_muons) >= 1 and len(pos_muons) >= 1 and
             len(neg_electrons) >= 1 and len(pos_electrons) >=1 )):
            Z_candidates = ZCandidatesMuEl(neg_muons, pos_muons,neg_electrons, pos_electrons, 80., 100.)
            #print Z_candidates
#print "length: ", len(Z_candidates)
            if len(Z_candidates) == 2 :
                logging.debug("Two Z candidates found")
                event.Cuts[self.name] = True
                self.NumberEventsPassedCut += 1.

                Z1 = Z_candidates[0]
                Z1_p1 = Z1[0].P4()
                Z1_p2 = Z1[1].P4()
                Z1Particles = [Z1_p1, Z1_p2]
                Z1 = Z1_p1 + Z1_p2
                Z2 = Z_candidates[1]
                Z2_p1 = Z2[0].P4()
                Z2_p2 = Z2[1].P4()
                Z2 = Z2_p1 + Z2_p2
                Z2Particles = [Z2_p1, Z1_p2]
                event.Z1Particles = Z1Particles
                event.Z1 = Z1
                event.Z2 = Z2
                event.Z2Particles = Z2Particles
                return True

        else :
            #print "no 2 Z in event"
            return False;
