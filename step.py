from functions import *
from ROOT import TH1D
class Step(object) :
    def __init__(self, step_name) :
        self.name = step_name
        self.NumberEventsPassedCut = 0.
    def Initialize(self, data_name, cut_flow) :
        self.DataName = data_name
        self.Cutflow = cut_flow

class AddPlot(Step) :
    def __init__(self, step_name) :
        super(AddPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(AddPlot, self).Initialize(data_name, cut_flow)
class LeptonIsolationPlot(AddPlot) :
    def __init__(self, step_name) :
        super(LeptonIsolationPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(LeptonIsolationPlot,self).Initialize(Histos, data_name, cut_flow)
        R_bin = 100
        R_min = 0.
        R_max = 20.
        self.Z1LeptonR = TH1D(self.DataName + self.Cutflow + "-Z1RL","Mass of tagging jets",R_bin, R_min, R_max)
        Histos.append(self.Z1LeptonR)
        self.Z2LeptonR = TH1D(self.DataName + self.Cutflow + "-Z2RL","Mass of tagging jets",R_bin, R_min, R_max)
        Histos.append(self.Z2LeptonR)

    def PerformStep(self, event, Histos, data_type) :
        self.Z1LeptonR.Fill(dR(event.Z1Particles[0], event.Z1Particles[1]))
        return True

class TaggingJetMassPlot(AddPlot) :
    def __init__(self, step_name) :
        super(TaggingJetMassPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(TaggingJetMassPlot,self).Initialize(Histos, data_name, cut_flow)
        M_bin = 100
        M_min = 0.
        M_max = 2000.
        self.TaggingJetMass = TH1D(self.DataName + self.Cutflow + "-TJMass","Mass of tagging jets",M_bin, M_min, M_max)
        Histos.append(self.TaggingJetMass)
        DY_bin = 100
        DY_min = 0.0
        DY_max = 9.0
        self.TaggingJetRapidityGap = TH1D(self.DataName + self.Cutflow + "-TJDY","Rapidity of tagging jets",DY_bin, DY_min, DY_max)
        H_bin = 60
        H_min = -15.0
        H_max = 15.0
        self.TaggingJetY1Y2 = TH1D(self.DataName + self.Cutflow + "-TJY1Y2"," y_1 * y_2 of  tagging jets",H_bin, H_min, H_max)
        Histos.append(self.TaggingJetY1Y2)

        PT1_bin = 100
        PT1_min = 0.0
        PT1_max = 500.0
        self.TaggingJetLeadingPT = TH1D(self.DataName + self.Cutflow + "-TJ1PT","PT of leading tagging jets",PT1_bin, PT1_min, PT1_max)
        Histos.append(self.TaggingJetLeadingPT)

        PT2_bin = 100
        PT2_min = 0.0
        PT2_max = 500.0
        self.TaggingJetSubleadingPT = TH1D(self.DataName + self.Cutflow + "-TJ2PT","PT of subleading tagging jets",PT2_bin, PT2_min, PT2_max)
        Histos.append(self.TaggingJetSubleadingPT)

        DPhi_bin = 100
        DPhi_min = 0.0
        DPhi_max = 7.0
        self.TaggingJetDPhi = TH1D(self.DataName + self.Cutflow + "TJDPhi","Delta phi of tagging jets",DPhi_bin, DPhi_min, DPhi_max)
        Histos.append(self.TaggingJetDPhi)

        Histos.append(self.TaggingJetRapidityGap)

    def PerformStep(self, event, Histos, data_type) :
        self.TaggingJetMass.Fill((event.TaggingJet1 + event.TaggingJet2).M())
        self.TaggingJetDPhi.Fill(event.TaggingJetDPhi)

        self.TaggingJetLeadingPT.Fill(event.TaggingJet1.Pt())
        self.TaggingJetSubleadingPT.Fill(event.TaggingJet2.Pt())
        self.TaggingJetY1Y2.Fill(event.TaggingJetY1Y2)
        self.TaggingJetRapidityGap.Fill(event.TaggingJetRapidityGap)

        return True
