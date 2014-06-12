from functions import *
from ROOT import TH1D, TH2D
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
        self.histos = {}
    def Initialize(self, external_histos, data_name, cut_flow) :
        super(AddPlot, self).Initialize(data_name, cut_flow)
        self._external_histos = external_histos
    def AddHistogram(self, local_name, ROOT_name, title, x_min, x_max, x_bin) :

        if local_name in self.histos :
            logging.error("Trying to add histgram with exact name")
        else :
            self.histos[local_name] = TH1D(self.DataName + self.Cutflow + "-" + ROOT_name,
                                           title, x_bin, x_min, x_max)
            self._external_histos.append(self.histos[local_name])

class ZPairMassPlot(AddPlot) :
    def __init__(self, step_name) :
        super(ZPairMassPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(ZPairMassPlot,self).Initialize(Histos, data_name, cut_flow)
        Jmul_bin  = 40
        Jmul_max  = 400.
        Jmul_min  = 0.
        self.ZPairMass = TH2D(self.DataName +self.Cutflow+ "-ZPairMass","Z pair mass distribution", Jmul_bin, Jmul_min, Jmul_max, Jmul_bin, Jmul_min, Jmul_max)
        self.ZPairMass.SetTitle("Z pair mass distribution;m_{Z_{1}};m_{Z_{2}}")
        Histos.append(self.ZPairMass)

    def PerformStep(self, event, Histos, data_type) :
        Z1Mass = 0.
        Z2Mass = 0.

        if data_type == "GEN" :
            if len(event.goodMuons) >= 2 and len(event.goodElectrons) >=2 :
                Z1Mass = (event.goodMuons[0].P4() + event.goodMuons[1].P4()).M()
                Z2Mass = (event.goodElectrons[0].P4() + event.goodElectrons[1].P4()).M()
        
        self.ZPairMass.Fill(Z1Mass, Z2Mass)
        return True


class JetMultiplicityPlot(AddPlot) :
    def __init__(self, step_name) :
        super(JetMultiplicityPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(JetMultiplicityPlot,self).Initialize(Histos, data_name, cut_flow)
        Jmul_bin  = 10
        Jmul_max  = 10
        Jmul_min  = 0
        self.JetMultiplicity = TH1D(self.DataName +self.Cutflow+ "-JetMultiplicity","y* of Z_1", Jmul_bin, Jmul_min, Jmul_max)
        Histos.append(self.JetMultiplicity)

    def PerformStep(self, event, Histos, data_type) :
        jetMul = 0
        if data_type == "SIM" :
            for jet in event.data.Jet : jetMul += 1
        if data_type == "Gen" :
            for jet in event.data.GenJet : jetMul+= 1
        self.JetMultiplicity.Fill(jetMul)
        return True

class ZeppenfeldVariablesPlot(AddPlot) :
    def __init__(self, step_name) :
        super(ZeppenfeldVariablesPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(ZeppenfeldVariablesPlot,self).Initialize(Histos, data_name, cut_flow)
        self.Cutflow = cut_flow
        Y1Star_bin  = 96
        Y1Star_max  = 4.8
        Y1Star_min  = 4.8
        self.Y1Star = TH1D(self.DataName +self.Cutflow+ "Y1Star","y* of Z_1", Y1Star_bin, Y1Star_min, Y1Star_max)
        Histos.append(self.Y1Star)

        Y2Star_bin  = 96
        Y2Star_max  = 4.8
        Y2Star_min  = 4.8
        self.Y2Star = TH1D(self.DataName +self.Cutflow+ "Y2Star","y* of Z_2", Y2Star_bin, Y2Star_min, Y2Star_max)

        Histos.append(self.Y2Star)

    def PerformStep(self, event, Histos, data_type) :
        Y1Star = event.Z1.Rapidity() - 0.5*(event.TaggingJet1.Rapidity() + event.TaggingJet2.Rapidity())
        Y2Star = event.Z2.Rapidity() - 0.5*(event.TaggingJet1.Rapidity() + event.TaggingJet2.Rapidity())
        event.Y1Star = Y1Star
        event.Y2Star = Y2Star
        self.Y1Star.Fill(Y1Star)
        self.Y2Star.Fill(Y2Star)
        return True

class LeptonIsolationPlot(AddPlot) :
    def __init__(self, step_name) :
        super(LeptonIsolationPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(LeptonIsolationPlot,self).Initialize(Histos, data_name, cut_flow)
        R_bin = 100
        R_min = 0.
        R_max = 9.
        self.Z1LeptonR = TH1D(self.DataName + self.Cutflow + "-Z1RL","Mass of tagging jets",R_bin, R_min, R_max)
        Histos.append(self.Z1LeptonR)
        self.Z2LeptonR = TH1D(self.DataName + self.Cutflow + "-Z2RL","Mass of tagging jets",R_bin, R_min, R_max)
        Histos.append(self.Z2LeptonR)
        self.LeptonR = TH1D(self.DataName + self.Cutflow + "-RLiLj","Mass of tagging jets",R_bin, R_min, R_max)
        Histos.append(self.LeptonR)
    def PerformStep(self, event, Histos, data_type) :
        self.Z1LeptonR.Fill(dR(event.Z1Particles[0], event.Z1Particles[1]))
        self.Z2LeptonR.Fill(dR(event.Z2Particles[0], event.Z2Particles[1]))
        for l1 in event.goodLeptons :
            for l2 in event.goodLeptons :
                if dR(l1.P4(),l2.P4()) <> 0. :
                    self.LeptonR.Fill(dR(l1.P4(), l2.P4()))
        return True

class TaggingJetKinematicsPlot(AddPlot) :
    def __init__(self, step_name) :
        super(TaggingJetKinematicsPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(TaggingJetKinematicsPlot,self).Initialize(Histos, data_name, cut_flow)
        ES_bin = 100
        ES_min = -10.
        ES_max = 10.
        self.TaggingJetEtaSum = TH1D(self.DataName + self.Cutflow + "-TJETASUM","Mass of tagging jets",ES_bin, ES_min, ES_max)
        Histos.append(self.TaggingJetEtaSum)
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

        self.AddHistogram("TJEta1Eta2", "Eta1Eta2", "#eta_{j_1}#eta_{j_2}", -10., 10., 10)


    def PerformStep(self, event, Histos, data_type) :
        self.TaggingJetMass.Fill((event.TaggingJet1 + event.TaggingJet2).M())
        self.TaggingJetDPhi.Fill(event.TaggingJetDPhi)

        self.TaggingJetLeadingPT.Fill(event.TaggingJet1.Pt())
        self.TaggingJetSubleadingPT.Fill(event.TaggingJet2.Pt())
        self.TaggingJetY1Y2.Fill(event.TaggingJetY1Y2)
        self.TaggingJetRapidityGap.Fill(event.TaggingJetRapidityGap)

        self.TaggingJetEtaSum.Fill(event.TaggingJetEtaSum)
        self.histos["TJEta1Eta2"].Fill(event.TaggingJet1.Eta() * event.TaggingJet2.Eta())
        return True

class AllLeptonPtEtaPlot(AddPlot) :
    def __init__(self, step_name) :
        super(AllLeptonPtEtaPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(AllLeptonPtEtaPlot,self).Initialize(Histos, data_name, cut_flow)

        eta_bin = 56
        eta_max  = 2.8
        eta_min  = -2.8
        self.LeadingLeptonEta = TH1D(self.DataName +self.Cutflow+ "LeadingLeptonEta","Eta of leading lepton", eta_bin, eta_min, eta_max)
        self.SubleadingLeptonEta = TH1D(self.DataName +self.Cutflow+ "SubleadingLeptonEta","Eta of subleading lepton",eta_bin, eta_min, eta_max)
        self.ThirdleadingLeptonEta = TH1D(self.DataName +self.Cutflow+ "ThirdleadingLeptonEta","Eta of third lepton",eta_bin, eta_min, eta_max)
        self.FourthleadingLeptonEta = TH1D(self.DataName +self.Cutflow+ "FourthleadingLeptonEta","Eta of fourth lepton",eta_bin, eta_min, eta_max)
        Histos.append(self.LeadingLeptonEta)
        Histos.append(self.SubleadingLeptonEta)
        Histos.append(self.ThirdleadingLeptonEta)
        Histos.append(self.FourthleadingLeptonEta)

        PT_bin = 400
        PT_max  = 400.
        PT_min  = 0.

        self.LeadingLeptonPT = TH1D(self.DataName +self.Cutflow+ "LeadingLeptonPT","PT of leading lepton", PT_bin, PT_min, PT_max)
        self.SubleadingLeptonPT = TH1D(self.DataName +self.Cutflow+ "SubleadingLeptonPT","PT of subleading lepton",PT_bin, PT_min, PT_max)
        self.ThirdleadingLeptonPT = TH1D(self.DataName +self.Cutflow+ "ThirdleadingLeptonPT","PT of third lepton",PT_bin, PT_min, PT_max - 200.)
        self.FourthleadingLeptonPT = TH1D(self.DataName +self.Cutflow+ "FourthleadingLeptonPT","PT of fourth lepton",PT_bin, PT_min, PT_max - 250.)

        Histos.append(self.LeadingLeptonPT)
        Histos.append(self.SubleadingLeptonPT)
        Histos.append(self.ThirdleadingLeptonPT)
        Histos.append(self.FourthleadingLeptonPT)

        M_bin = 40
        M_min = 0.
        M_max = 800.
        self.FourLeptonMass = TH1D(self.DataName +self.Cutflow+ "FourLeptonMass","Mass of 4 lepton system",M_bin, M_min, M_max)

        Histos.append(self.FourLeptonMass)


    def PerformStep(self, event, Histos, data_type) :

        fourLepton = TLorentzVector(0,0,0,0)
        leptons = []
        if data_type == "SIM" :

            for mu in event.data.Muon :
                leptons.append(mu)
            for el in event.data.Electron :
                leptons.append(el)
        leptons.sort(key=lambda x: x.PT, reverse=True)
        Mass = 0.
        if len(leptons) >= 4 :
            for lepton in leptons[:3] :
                fourLepton += lepton.P4()
                Mass = fourLepton.M()
        if len(leptons) >= 1 :
            self.LeadingLeptonEta.Fill(leptons[0].Eta)
            self.LeadingLeptonPT.Fill(leptons[0].PT)
        if len(leptons) >= 2 :
            self.SubleadingLeptonEta.Fill(leptons[1].Eta)
            self.SubleadingLeptonPT.Fill(leptons[1].PT)
        if len(leptons) >= 3 :
            self.ThirdleadingLeptonEta.Fill(leptons[2].Eta)
            self.ThirdleadingLeptonPT.Fill(leptons[2].PT)
        if len(leptons) >=4 :
            self.FourthleadingLeptonEta.Fill(leptons[3].Eta)
            self.FourthleadingLeptonPT.Fill(leptons[3].PT)
        self.FourLeptonMass.Fill(Mass)
        return True



class GoodLeptonPtEtaPlot(AddPlot) :
    def __init__(self, step_name) :
        super(GoodLeptonPtEtaPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(GoodLeptonPtEtaPlot,self).Initialize(Histos, data_name, cut_flow)

        eta_bin = 56
        eta_max  = 2.8
        eta_min  = -2.8
        self.LeadingGoodLeptonEta = TH1D(self.DataName +self.Cutflow+ "LeadingGoodLeptonEta","Eta of leading lepton", eta_bin, eta_min, eta_max)
        self.SubleadingGoodLeptonEta = TH1D(self.DataName +self.Cutflow+ "SubleadingGoodLeptonEta","Eta of subleading lepton",eta_bin, eta_min, eta_max)
        self.ThirdleadingGoodLeptonEta = TH1D(self.DataName +self.Cutflow+ "ThirdleadingGoodLeptonEta","Eta of third lepton",eta_bin, eta_min, eta_max)
        self.FourthleadingGoodLeptonEta = TH1D(self.DataName +self.Cutflow+ "FourthleadingGoodLeptonEta","#eta of the fourth lepton;#eta;leptons",eta_bin, eta_min, eta_max)
        Histos.append(self.LeadingGoodLeptonEta)
        Histos.append(self.SubleadingGoodLeptonEta)
        Histos.append(self.ThirdleadingGoodLeptonEta)
        Histos.append(self.FourthleadingGoodLeptonEta)

        PT_bin = 400
        PT_max  = 400.
        PT_min  = 0.

        self.LeadingGoodLeptonPT = TH1D(self.DataName +self.Cutflow+ "LeadingGoodLeptonPT","PT of leading lepton", PT_bin, PT_min, PT_max)
        self.SubleadingGoodLeptonPT = TH1D(self.DataName +self.Cutflow+ "SubleadingGoodLeptonPT","PT of subleading lepton",PT_bin, PT_min, PT_max)
        self.ThirdleadingGoodLeptonPT = TH1D(self.DataName +self.Cutflow+ "ThirdleadingGoodLeptonPT","PT of third lepton",PT_bin, PT_min, PT_max)
        self.FourthleadingGoodLeptonPT = TH1D(self.DataName +self.Cutflow+ "FourthleadingGoodLeptonPT","p_{T} of the fourth lepton;p_{T};leptons",75, PT_min, 150.)

        Histos.append(self.LeadingGoodLeptonPT)
        Histos.append(self.SubleadingGoodLeptonPT)
        Histos.append(self.ThirdleadingGoodLeptonPT)
        Histos.append(self.FourthleadingGoodLeptonPT)

        M_bin = 40
        M_min = 0.
        M_max = 800.
        self.FourGoodLeptonMass = TH1D(self.DataName +self.Cutflow+ "FourGoodLeptonMass","Mass of 4 lepton system",M_bin, M_min, M_max)

        Histos.append(self.FourGoodLeptonMass)


    def PerformStep(self, event, Histos, data_type) :

        fourLepton = TLorentzVector(0,0,0,0)
        for lepton in event.goodLeptons[:3] :
            fourLepton += lepton.P4()
            Mass = fourLepton.M()

        self.LeadingGoodLeptonEta.Fill(event.goodLeptons[0].Eta)
        self.SubleadingGoodLeptonEta.Fill(event.goodLeptons[1].Eta)
        self.ThirdleadingGoodLeptonEta.Fill(event.goodLeptons[2].Eta)
        self.FourthleadingGoodLeptonEta.Fill(event.goodLeptons[3].Eta)

        self.LeadingGoodLeptonPT.Fill(event.goodLeptons[0].PT)
        self.SubleadingGoodLeptonPT.Fill(event.goodLeptons[1].PT)
        self.ThirdleadingGoodLeptonPT.Fill(event.goodLeptons[2].PT)
        self.FourthleadingGoodLeptonPT.Fill(event.goodLeptons[3].PT)
        self.FourGoodLeptonMass.Fill(Mass)
        return True

class TaggingJetZKinematicsPlot(AddPlot) :
    def __init__(self, step_name) :
        super(TaggingJetZKinematicsPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(TaggingJetZKinematicsPlot,self).Initialize(Histos, data_name, cut_flow)
        DPhi_bin = 40
        DPhi_max  = 7.
        DPhi_min  = -7.


        self.DPhiZ1j1 = TH1D(self.DataName +self.Cutflow+ "DPhiZ1j1","PT of leading lepton", DPhi_bin, DPhi_min, DPhi_max)
        self.DPhiZ1j2 = TH1D(self.DataName +self.Cutflow+ "DPhiZ1j2","PT of leading lepton", DPhi_bin, DPhi_min, DPhi_max)
        self.DPhi4lj1 = TH1D(self.DataName +self.Cutflow+ "DPhi4lj1","PT of leading lepton", DPhi_bin, DPhi_min, DPhi_max)
        Histos.append(self.DPhiZ1j1)
        Histos.append(self.DPhiZ1j2)
        Histos.append(self.DPhi4lj1)

    def PerformStep(self, event, Histos, data_type) :
        DPhiZ1j1 = event.Z1.Phi() - event.TaggingJet1.Phi()
        DPhiZ1j2 = event.Z1.Phi() - event.TaggingJet2.Phi()
        
        
        DPhi4lj1 = (event.Z1 + event.Z2).Phi() - event.TaggingJet1.Phi()
        self.DPhiZ1j1.Fill(DPhiZ1j1)
        self.DPhiZ1j2.Fill(DPhiZ1j2)
        self.DPhi4lj1.Fill(DPhi4lj1)
         
        return True
