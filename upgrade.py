import logging

from functions import *
from ROOT import TH1D, TH2D

class ZJetsKinematics :
    def __init__(self, step_name) :
        logging.debug("ZJetsKinematics constructor")
        self.name = step_name
        self.NumberEventsPassedCut = 0.
    def Initialize(self, Histos, data_name, cut_flow) :
        self.DataName = data_name

        R_bin  = 40
        R_max  = 8
        R_min  = 0.

        self.RJ1L1 = TH1D(self.DataName + "R_J1_L1","R difference of leading tagging jet and leading lepton", R_bin, R_min, R_max)
        self.RJ1L2 = TH1D(self.DataName + "R_J1_L2","R difference of leading tagging jet and subleading lepton",R_bin, R_min, R_max)

        Histos.append(self.RJ1L1)
        Histos.append(self.RJ1L2)
    def IsCut(self) :
        return False


    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        event.Cuts[self.name] = False
        self.NumberEventsPassedCut += 1.

        J1 = event.goodJets[0].P4()
        L1 = event.goodLeptons[0].P4()
        L2 = event.goodLeptons[1].P4()

        self.RJ1L1.Fill(dR(J1, L1))
        self.RJ1L2.Fill(dR(J1, L2))

        return True

class ZeppenfeldVar :
    def __init__(self, step_name) :
        logging.debug("Zeppenfeld constructor")
        self.name = step_name
        self.NumberEventsPassedCut = 0.

    def Initialize(self, Histos, data_name, cut_flow) :
        self.DataName = data_name
        self.Cutflow = cut_flow
        Y1Star_bin  = 30
        Y1Star_max  = 4.8
        Y1Star_min  = 4.8
        self.Y1Star = TH1D(self.DataName +self.Cutflow+ "Y1Star","y* of Z_1", Y1Star_bin, Y1Star_min, Y1Star_max)
        Histos.append(self.Y1Star)

        Y2Star_bin  = 30
        Y2Star_max  = 4.8
        Y2Star_min  = 4.8
        self.Y2Star = TH1D(self.DataName +self.Cutflow+ "Y2Star","y* of Z_2", Y2Star_bin, Y2Star_min, Y2Star_max)

        Histos.append(self.Y2Star)


    def IsCut(self) :
        return True
    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        event.Cuts[self.name] = True
        self.NumberEventsPassedCut += 1.

        Y1Star = event.Z1.Rapidity() - 0.5*(event.TaggingJet1.Rapidity() + event.TaggingJet2.Rapidity())
        Y2Star = event.Z2.Rapidity() - 0.5*(event.TaggingJet1.Rapidity() + event.TaggingJet2.Rapidity())

        self.Y1Star.Fill(Y1Star)
        self.Y2Star.Fill(Y2Star)
        return True


class ZKinematics :
    def __init__(self, step_name) :
        logging.debug("Zkinematics constructor")
        self.name = step_name
        self.NumberEventsPassedCut = 0.

    def Initialize(self, Histos, data_name, cut_flow) :
        self.DataName = data_name
        self.Cutflow = cut_flow
        eta_bin  = 30
        eta_max  = 4.8
        eta_min  = 4.8
        self.Z1Eta = TH1D(self.DataName +self.Cutflow+ "Z1Eta","Eta of Z_1", eta_bin, eta_min, eta_max)
        self.Z2Eta = TH1D(self.DataName +self.Cutflow+ "Z2Eta","Eta of Z_1",eta_bin, eta_min, eta_max)

        Histos.append(self.Z1Eta)
        Histos.append(self.Z2Eta)

        PT_bin = 100
        PT_max  = 200.
        PT_min  = 0.

        self.Z1PT = TH1D(self.DataName +self.Cutflow+ "Z1PT","PT of Z_1", PT_bin, PT_min, PT_max)
        self.Z2PT = TH1D(self.DataName +self.Cutflow +"Z2PT","PT of Z_2",PT_bin, PT_min, PT_max)

        Histos.append(self.Z1PT)
        Histos.append(self.Z2PT)

    def IsCut(self) :
        return False
    def PerformStep(self, event, Histos, data_type) :
        logging.debug("Called ApplyCut of " + self.name + " cut" )
        event.Cuts[self.name] = False
        self.NumberEventsPassedCut += 1.

        self.Z1Eta.Fill(event.Z1.Eta())
        self.Z1PT.Fill(event.Z1.Pt())
        self.Z2Eta.Fill(event.Z2.Eta())
        self.Z2PT.Fill(event.Z2.Pt())

        return True

