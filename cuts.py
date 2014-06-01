import logging
from functions import *
from ROOT import TH2D
class Start :
    def __init__(self, cut_name, Histos, data_name) :
        logging.debug("Start constructor")
        self.name = cut_name
        self.DataName = data_name
        self.NumberEventsPassedCut = 0.

    def IsCut(self) :
        return True

    def ApplyCut(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        event.Cuts[self.name] = True
        self.NumberEventsPassedCut += 1.
        return True

class LeptonIsolation :
    def __init__(self, cut_name, Histos, data_name) :
        logging.debug("LeptonIsolation constructor")
        self.name = cut_name
        self.DataName = data_name
        self.NumberEventsPassedCut = 0.

    def IsCut(self) :
        return True

    def ApplyCut(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        for lepton1 in event.goodLeptons :
            for lepton2 in event.goodLeptons :
                if lepton1 <> lepton2 and dR(lepton1, lepton2) <= 0.5 :
                    event.Cuts[self.name] = False
                    return False

        event.Cuts[self.name] = True
        self.NumberEventsPassedCut += 1.
        return True


class ZJetsKinematics :
    def __init__(self, cut_name, Histos, data_name) :
        self.DataName = data_name
        logging.debug("ZJetsKinematics constructor")
        self.name = cut_name
        self.NumberEventsPassedCut = 0.
        R_bin = 40
        R_max  = 8
        R_min  = 0.
        
        self.RJ1L1 = TH1D(self.DataName + "R_J1_L1","R difference of leading tagging jet and leading lepton", R_bin, R_min, R_max)
        self.RJ1L2 = TH1D(self.DataName + "R_J1_L2","R difference of leading tagging jet and subleading lepton",R_bin, R_min, R_max)

        Histos.append(self.RJ1L1)
        Histos.append(self.RJ1L2)
    def IsCut(self) :
        return False


    def ApplyCut(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        event.Cuts[self.name] = False
        self.NumberEventsPassedCut += 1.

        J1 = event.goodJets[0]
        L1 = event.goodLeptons[0]
        L2 = event.goodLeptons[1]

        self.RJ1L1.Fill(dR(J1, L1))
        self.RJ1L2.Fill(dR(J1, L2))

        return True

class DefineTaggingJets :
    def __init__(self, cut_name, Histos, data_name) :
        self.DataName = data_name
        logging.debug("Tagging jets constructor")
        self.name = cut_name
        self.NumberEventsPassedCut = 0.

        M_bin = 100
        M_min = 0.
        M_max = 2000.
        self.TaggingJetMass = TH1D(self.DataName + "TaggingJetMass","Mass of tagging jets",M_bin, M_min, M_max)

        DY_bin = 100
        DY_min = 0.0
        DY_max = 9.0
        self.TaggingJetRapidityGap = TH1D(self.DataName + "TaggingJetRapidityGap","Rapidity of tagging jets",DY_bin, DY_min, DY_max)
        H_bin = 40
        H_min = -20.0
        H_max = 20.0
        self.TaggingJetY1Y2 = TH1D(self.DataName + "TaggingJetY1Y2"," y_1 * y_2 of  tagging jets",H_bin, H_min, H_max)
        Histos.append(self.TaggingJetY1Y2)

        PT1_bin = 100
        PT1_min = 0.0
        PT1_max = 500.0
        self.TaggingJetLeadingPT = TH1D(self.DataName + "TaggingJetLeadingPT","PT of leading tagging jets",PT1_bin, PT1_min, PT1_max)
        Histos.append(self.TaggingJetLeadingPT)

        PT2_bin = 100
        PT2_min = 0.0
        PT2_max = 500.0
        self.TaggingJetSubleadingPT = TH1D(self.DataName + "TaggingJetSubleadingPT","PT of subleading tagging jets",PT2_bin, PT2_min, PT2_max)
        Histos.append(self.TaggingJetSubleadingPT)

        DPhi_bin = 100
        DPhi_min = 0.0
        DPhi_max = 7.0
        self.TaggingJetDPhi = TH1D(self.DataName + "TaggingJetDPhi","Delta phi of tagging jets",DPhi_bin, DPhi_min, DPhi_max)
        Histos.append(self.TaggingJetDPhi)

        Histos.append(self.TaggingJetRapidityGap)
        Histos.append(self.TaggingJetMass)


    def IsCut(self) :
        return True

    def ApplyCut(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" + "in mode " + data_type )
        event.Cuts[self.name] = False

        goodJets = []

        if data_type == "SIM" :
            for jet in event.data.Jet :
                goodJets.append(jet)
        if data_type == "GEN" :
            for jet in event.GenJet :
                goodJets.append(jet)
        goodJets.sort(key=lambda x: x.PT, reverse=True)

        if len(goodJets) >= 2 :
            event.Cuts[self.name] = True
            TaggingJet1 = goodJets[0]
            TaggingJet2 = goodJets[1]
            event.TaggingJet1 = TaggingJet1.P4()
            event.TaggingJet2 = TaggingJet2.P4()
            event.goodJets = goodJets
            self.NumberEventsPassedCut += 1.0

            rapidity_product = TaggingJet1.P4().Rapidity() * TaggingJet2.P4().Rapidity()
            gap = abs(TaggingJet1.P4().Rapidity() - TaggingJet2.P4().Rapidity())

            DPhi = abs(TaggingJet1.Phi - TaggingJet2.Phi)
            self.TaggingJetDPhi.Fill(DPhi)

            self.TaggingJetLeadingPT.Fill(TaggingJet1.PT)
            self.TaggingJetSubleadingPT.Fill(TaggingJet2.PT)

            self.TaggingJetY1Y2.Fill(rapidity_product)
            self.TaggingJetMass.Fill((TaggingJet1.P4() + TaggingJet2.P4()).M())
            self.TaggingJetRapidityGap.Fill(gap)

            return True
        return False

class TaggingJetInvariantMass :
    def __init__(self, cut_name, Histos, data_name) :
        self.DataName = data_name
        logging.debug("TaggingJetInvariantMass constructor")
        self.name = cut_name
        self.NumberEventsPassedCut = 0.

    def IsCut(self) :
        return True
    def ApplyCut(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        
        mass = (event.TaggingJet1 + event.TaggingJet2).M()
        if mass > 800. :
            event.Cuts[self.name] = True
            self.NumberEventsPassedCut += 1.
            return True
        else :
            event.Cuts[self.name] = False
            return False

        return True




class ZeppenfeldVar :
    def __init__(self, cut_name, Histos, data_name) :
        self.DataName = data_name
        logging.debug("Zeppenfeld constructor")
        self.name = cut_name
        self.NumberEventsPassedCut = 0.

        Y1Star_bin = 30
        Y1Star_max  = 4.8
        Y1Star_min  = 4.8
        self.Y1Star = TH1D(self.DataName + "Y1Star","y* of Z_1", Y1Star_bin, Y1Star_min, Y1Star_max)

        Histos.append(self.Y1Star)

        Y2Star_bin = 30
        Y2Star_max  = 4.8
        Y2Star_min  = 4.8
        self.Y2Star = TH1D(self.DataName + "Y2Star","y* of Z_2", Y2Star_bin, Y2Star_min, Y2Star_max)

        Histos.append(self.Y2Star)


    def IsCut(self) :
        return True
    def ApplyCut(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        event.Cuts[self.name] = True
        self.NumberEventsPassedCut += 1.

        Y1Star = event.Z1.Rapidity() - 0.5*(event.TaggingJet1.Rapidity() + event.TaggingJet2.Rapidity())
        Y2Star = event.Z2.Rapidity() - 0.5*(event.TaggingJet1.Rapidity() + event.TaggingJet2.Rapidity())

        self.Y1Star.Fill(Y1Star)
        self.Y2Star.Fill(Y2Star)
        return True




class ZKinematics :
    def __init__(self, cut_name, Histos, data_name) :
        self.DataName = data_name
        logging.debug("Zkinematics constructor")
        self.name = cut_name
        self.NumberEventsPassedCut = 0.

        eta_bin = 30
        eta_max  = 4.8
        eta_min  = 4.8
        self.Z1Eta = TH1D(self.DataName + "Z1Eta","Eta of Z_1", eta_bin, eta_min, eta_max)
        self.Z2Eta = TH1D(self.DataName + "Z2Eta","Eta of Z_1",eta_bin, eta_min, eta_max)

        Histos.append(self.Z1Eta)
        Histos.append(self.Z2Eta)

        PT_bin = 100
        PT_max  = 100.
        PT_min  = 0.

        self.Z1PT = TH1D(self.DataName + "Z1PT","PT of Z_1", PT_bin, PT_min, PT_max)
        self.Z2PT = TH1D(self.DataName + "Z2PT","PT of Z_2",PT_bin, PT_min, PT_max)

        Histos.append(self.Z1PT)
        Histos.append(self.Z2PT)

    def IsCut(self) :
        return False
    def ApplyCut(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        event.Cuts[self.name] = False
        self.NumberEventsPassedCut += 1.

        self.Z1Eta.Fill(event.Z1.Eta())

        self.Z1PT.Fill(event.Z1.Eta())
        return True



#class AcceptanceCutSet :
#    def __init__(self, PTcut, ETAcut)

class LeptonAcceptanceAnalysis :
    def __init__(self, cut_name, Histos, data_name) :
        self.DataName = data_name
        logging.debug("LeptonAcceptance constructor")
        self.name = cut_name
        self.NumberEventsPassedCut = 0.

        eta_bin = 5
        eta_max  = 2.5
        eta_min  = -2.5

        PT_bin = 10
        PT_max  = 300.
        PT_min  = 0.

        self.ElectronEta = TH1D(self.DataName + "ElectronEta","Number of electrons", eta_bin, eta_min, eta_max)
        self.MuonEta = TH1D(self.DataName + "MuonEta","Number of muons",eta_bin, eta_min, eta_max)
        self.ElectronPtEta = TH2D(self.DataName + "Electron_Pt_Eta", "Electron efficency", PT_bin, PT_min, PT_max, eta_bin, eta_min, eta_max)
        Histos.append(self.ElectronPtEta)
        self.MuonPtEta = TH2D(self.DataName + "Muon_Pt_Eta", "Muon efficency", PT_bin, PT_min, PT_max, eta_bin, eta_min, eta_max)
        
        Histos.append(self.MuonPtEta)
        Histos.append(self.ElectronEta)
        Histos.append(self.MuonEta)

        self.ElectronPT = TH1D(self.DataName + "ElectronPT","PT of leading lepton", PT_bin, PT_min, PT_max)
        self.MuonPT = TH1D(self.DataName + "MuonPT","PT of subleading lepton",PT_bin, PT_min, PT_max)

        Histos.append(self.ElectronPT)
        Histos.append(self.MuonPT)

    def IsCut(self) :
        return True




    def ApplyCut(self, event, Histos, data_type) :
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

        for mu in goodMuons :
            self.MuonPtEta.Fill(mu.PT, mu.Eta)
            self.MuonEta.Fill(mu.Eta)
            self.MuonPT.Fill(mu.PT)

        for el in goodElectrons :
            self.ElectronPtEta.Fill(el.PT, el.Eta)
            self.ElectronEta.Fill(el.Eta)
            self.ElectronPT.Fill(el.PT)
        return True



class LeptonAcceptance :
    def __init__(self, cut_name, Histos, data_name) :
        self.DataName = data_name
        logging.debug("LeptonAcceptance constructor")
        self.name = cut_name
        self.NumberEventsPassedCut = 0.

        eta_bin = 56
        eta_max  = 2.8
        eta_min  = -2.8
        self.LeadingLeptonEta = TH1D(self.DataName + "LeadingLeptonEta","Eta of leading lepton", eta_bin, eta_min, eta_max)
        self.SubleadingLeptonEta = TH1D(self.DataName + "SubleadingLeptonEta","Eta of subleading lepton",eta_bin, eta_min, eta_max)
        self.ThirdleadingLeptonEta = TH1D(self.DataName + "ThirdleadingLeptonEta","Eta of third lepton",eta_bin, eta_min, eta_max)
        self.FourthleadingLeptonEta = TH1D(self.DataName + "FourthleadingLeptonEta","Eta of fourth lepton",eta_bin, eta_min, eta_max)

        Histos.append(self.LeadingLeptonEta)
        Histos.append(self.SubleadingLeptonEta)
        Histos.append(self.ThirdleadingLeptonEta)
        Histos.append(self.FourthleadingLeptonEta)

        PT_bin = 400
        PT_max  = 400.
        PT_min  = 0.

        self.LeadingLeptonPT = TH1D(self.DataName + "LeadingLeptonPT","PT of leading lepton", PT_bin, PT_min, PT_max)
        self.SubleadingLeptonPT = TH1D(self.DataName + "SubleadingLeptonPT","PT of subleading lepton",PT_bin, PT_min, PT_max)
        self.ThirdleadingLeptonPT = TH1D(self.DataName + "ThirdleadingLeptonPT","PT of third lepton",PT_bin, PT_min, PT_max)
        self.FourthleadingLeptonPT = TH1D(self.DataName + "FourthleadingLeptonPT","PT of fourth lepton",PT_bin, PT_min, PT_max)

        Histos.append(self.LeadingLeptonPT)
        Histos.append(self.SubleadingLeptonPT)
        Histos.append(self.ThirdleadingLeptonPT)
        Histos.append(self.FourthleadingLeptonPT)

        M_bin = 40
        M_min = 0.
        M_max = 800.
        self.FourLeptonMass = TH1D(self.DataName + "FourLeptonMass","Mass of 4 lepton system",M_bin, M_min, M_max)

        Histos.append(self.FourLeptonMass)

    def IsCut(self) :
        return True




    def ApplyCut(self, event, Histos, data_type) :
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

            fourLepton = TLorentzVector(0,0,0,0)
            for lepton in goodLeptons :
                fourLepton += lepton.P4()
            Mass = fourLepton.M()

            self.LeadingLeptonEta.Fill(goodLeptons[0].Eta)
            self.SubleadingLeptonEta.Fill(goodLeptons[1].Eta)
            self.ThirdleadingLeptonEta.Fill(goodLeptons[2].Eta)
            self.FourthleadingLeptonEta.Fill(goodLeptons[3].Eta)

            self.LeadingLeptonPT.Fill(goodLeptons[0].PT)
            self.SubleadingLeptonPT.Fill(goodLeptons[1].PT)
            self.ThirdleadingLeptonPT.Fill(goodLeptons[2].PT)
            self.FourthleadingLeptonPT.Fill(goodLeptons[3].PT)
            
            self.FourLeptonMass.Fill(Mass)
            
            return True
        else:
                for mu in event.data.Muon :
                    print "Muon pt: ", mu.PT , "eta: ", mu.Eta
                for el in event.data.Electron :
                    print "Elec PT: " , el.PT , "eta ", el.Eta
        return False

class ZPair :
    def __init__(self, cut_name, Histos, data_name) :
        self.DataName = data_name
        self.name = cut_name
        self.NumberEventsPassedCut = 0.
    def IsCut(self) :
        return True


    def ApplyCut(self, event, Histos, data_type) :
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
            Z_candidates = ZCandidates(neg_muons, pos_muons,neg_electrons, pos_electrons, 1., 1000.)
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
