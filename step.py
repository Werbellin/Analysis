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

    def Add2DHistogram(self, local_name, ROOT_name, title, x_min, x_max, x_bin, y_min, y_max, y_bin) :
        if local_name in self.histos :
            logging.error("Trying to add histgram with exact name")
        else :
            self.histos[local_name] = TH2D(self.DataName + self.Cutflow + "-" + ROOT_name,
                                           title, x_bin, x_min, x_max, y_bin, y_min, y_max)
            self._external_histos.append(self.histos[local_name])


class ZPairMassPlot(AddPlot) :
    def __init__(self, step_name) :
        super(ZPairMassPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(ZPairMassPlot,self).Initialize(Histos, data_name, cut_flow)
        Jmul_bin  = 20
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


class ZZRapidityPlot(AddPlot) :
    def __init__(self, step_name) :
        super(ZZRapidityPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(ZZRapidityPlot,self).Initialize(Histos, data_name, cut_flow)
        self.AddHistogram("ZZRapidity", "ZZRapidity", "Rapidity of ZZ system;y_{ZZ};Events/0.50 a.u.", 0., 5., 10)

    def PerformStep(self, event, Histos, data_type) :
        ZZRapidity = (event.Z1Particles[0] + event.Z1Particles[1] + event.Z2Particles[0] + event.Z2Particles[1]).Rapidity()
        ZZRapidity = abs(ZZRapidity)
        event.ZZRapidity = ZZRapidity
        self.histos["ZZRapidity"].Fill(ZZRapidity)
        return True

class YEtaPlot(AddPlot) :
    def __init__(self, step_name) :
        super(YEtaPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(YEtaPlot,self).Initialize(Histos, data_name, cut_flow)
        self.Add2DHistogram("YEta", "YEta", "y vs eta;y_{ZZ};Events/0.50 a.u.", -5., 5., 50, -5., 5., 50)

    def PerformStep(self, event, Histos, data_type) :
        self.histos["YEta"].Fill(event.TaggingJet1.Rapidity(), event.TaggingJet1.Eta())
        return True

class YStarZZPlot(AddPlot) :
    def __init__(self, step_name) :
        super(YStarZZPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(YStarZZPlot,self).Initialize(Histos, data_name, cut_flow)
        self.AddHistogram("YStarZZ", "YStarZZ", "Zeppenfeld variable of ZZ system;y*_{ZZ};Events/0.50 a.u.", 0., 5., 10)

    def PerformStep(self, event, Histos, data_type) :
        ZZRapidity = (event.Z1Particles[0] + event.Z1Particles[1] + event.Z2Particles[0] + event.Z2Particles[1]).Rapidity()
        YZZ = ZZRapidity - 0.5 * (event.TaggingJet1.Rapidity() + event.TaggingJet2.Rapidity())
        YZZ = abs(YZZ)
        self.histos["YStarZZ"].Fill(YZZ)
        return True

class Z1LeptonPtPlot(AddPlot) :
    def __init__(self, step_name) :
        super(Z1LeptonPtPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(Z1LeptonPtPlot,self).Initialize(Histos, data_name, cut_flow)
        self.AddHistogram("Z1LeptonPt", "-Z1LeptonPt", "p_{T} of Z_{1};p_{T}", 0., 140., 40)

    def PerformStep(self, event, Histos, data_type) :
        Z1Mass = 0.
        self.histos["Z1LeptonPt"].Fill(event.Z1Particles[0].Pt())
        self.histos["Z1LeptonPt"].Fill(event.Z1Particles[1].Pt())
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
            for jet in event.data.Jet : 
                if jet.PT > 20.0 :
                    jetMul += 1
        if data_type == "Gen" :
            for jet in event.data.GenJet : 
                if jet.PT > 20.0 :
                    jetMul+= 1
        self.JetMultiplicity.Fill(jetMul)
        return True

class ZeppenfeldVariablesPlot(AddPlot) :
    def __init__(self, step_name) :
        super(ZeppenfeldVariablesPlot, self).__init__(step_name)

    def Initialize(self, Histos, data_name, cut_flow) :
        super(ZeppenfeldVariablesPlot,self).Initialize(Histos, data_name, cut_flow)
        self.Cutflow = cut_flow
        Y1Star_bin  = 12
        Y1Star_max  = 4.8
        Y1Star_min  = 4.8
        self.Y1Star = TH1D(self.DataName +self.Cutflow+ "Y1Star","y* of Z_1", Y1Star_bin, Y1Star_min, Y1Star_max)
        Histos.append(self.Y1Star)

        Y2Star_bin  = 12
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
        self.AddHistogram("TJ1E", "TJ1E", "E of leading jet;E_{j_{1}};Events/10 GeV", 0., 800., 80)
        self.AddHistogram("TJ2E", "TJ2E", "E of subleading jet;E_{j_{2}};Events/10 GeV", 0., 800., 80)
        self.AddHistogram("TJETASUM", "TJETASUM", "ll_{j_{1}} + ta_{j_{2}};#eta_{j1} + #eta_{j2};Events/1 a.u.", -10., 10., 20)
        self.AddHistogram("TJMASS", "TJMASS", "Invariant mass of tagging jets;m_{jj};Events/20 GeV", 0., 2000., 100)
        self.AddHistogram("TJDY", "TJDY", "Rapidity gap of tagging jets;#delta(j,j);Events/0.3 a.u. ", 0., 9., 30)
        self.AddHistogram("TJY1Y2", "TJY1Y2", "y_{1} * y_{2} of  tagging jets;y_{j_{1}}*y_{j_{2}};Events/0.5 a.u.", -15., 15., 60)
        self.AddHistogram("TJ1PT", "TJ1PT", "p_T of  leading tagging jet;p_{T}^{j_{1}};Events/10 GeV", 0., 500., 50)
        self.AddHistogram("TJ2PT", "TJ2PT", "p_T of  subleading tagging jet;p_{T}^{j_{2}};Events/10 GeV", 0., 500., 50)
        self.AddHistogram("TJDPhi", "TJDPhi", "#phi_{1}-#phi_{2} for  tagging jets;#Delta#phi_{jj};Events/0.07 a.u.", 0., 7., 100)
        self.AddHistogram("TJEta1Eta2", "Eta1Eta2", "ll_{j_1}ta_{j_2};Events/2.0 a.u.", -10., 10., 10)
        self.AddHistogram("TJDEtaDY", "TJDEtDaY", "#deltall - #delta y;#deltata - #delta y;Events/1.0 a.u.",-10., 10., 20)

    def PerformStep(self, event, Histos, data_type) :
        self.histos["TJMASS"].Fill((event.TaggingJet1 + event.TaggingJet2).M())
        self.histos["TJDPhi"].Fill(event.TaggingJetDPhi)
        self.histos["TJ1PT"].Fill(event.TaggingJet1.Pt())
        self.histos["TJ2PT"].Fill(event.TaggingJet2.Pt())
        self.histos["TJY1Y2"].Fill(event.TaggingJetY1Y2)
        self.histos["TJDY"].Fill(event.TaggingJetRapidityGap)
        self.histos["TJETASUM"].Fill(event.TaggingJetEtaSum)
        self.histos["TJEta1Eta2"].Fill(event.TaggingJet1.Eta() * event.TaggingJet2.Eta())
        self.histos["TJDEtaDY"].Fill(abs(event.TaggingJetRapidityGap - event.TaggingJetEtaGap))
        self.histos["TJ1E"].Fill(event.TaggingJet1.E())
        self.histos["TJ2E"].Fill(event.TaggingJet2.E())
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

        self.AddHistogram("GL1ETA", "GL1ETA", "Leading good lepton #eta;\eta_{\ell_{1}};Events/0.1 a.u.", -2.8, 2.8, 56)
        self.AddHistogram("GL2ETA", "GL2ETA", "Subleading good lepton #eta;\eta_{\ell _{2}};Events/0.1 a.u.", -2.8, 2.8, 56)
        self.AddHistogram("GL3ETA", "GL3ETA", "3rd good lepton #eta;\eta_{\ell_{3}};Events/0.1 a.u.", -2.8, 2.8, 56)
        self.AddHistogram("GL4ETA", "GL4ETA", "4th good lepton #eta;\eta_{\ell_4};Events/0.1 a.u.", -2.8, 2.8, 56)

        self.AddHistogram("GL1PT", "GL1PT", "Leading good lepton p_{T};p_{T}^{l_{1}};Events/1 GeV", 0., 400., 400)
        self.AddHistogram("GL2PT", "GL2PT", "Subleading good lepton p_{T};p_{T}^{l2};Events/1 GeV", 0., 400., 400)
        self.AddHistogram("GL3PT", "GL3PT", "3rd good lepton p_{T};p_{T}^{l3};Events/1 GeV", 0., 400., 400)
        self.AddHistogram("GL4PT", "GL4PT", "4th good lepton p_{T};p_{T}^{l4};Events/2 GeV", 0., 150., 75)

        self.AddHistogram("GLMASS", "GLMASS", "Good 4 lepton mass;m_{4\ell};Events/20 GeV", 0., 800., 40)

    def PerformStep(self, event, Histos, data_type) :
        fourLepton = TLorentzVector(0,0,0,0)
        for lepton in event.goodLeptons[:3] :
            fourLepton += lepton.P4()
            Mass = fourLepton.M()

        self.histos["GL1ETA"].Fill(event.goodLeptons[0].Eta)
        self.histos["GL2ETA"].Fill(event.goodLeptons[1].Eta)
        self.histos["GL3ETA"].Fill(event.goodLeptons[2].Eta)
        self.histos["GL4ETA"].Fill(event.goodLeptons[3].Eta)

        self.histos["GL1PT"].Fill(event.goodLeptons[0].PT)
        self.histos["GL2PT"].Fill(event.goodLeptons[1].PT)
        self.histos["GL3PT"].Fill(event.goodLeptons[2].PT)
        self.histos["GL4PT"].Fill(event.goodLeptons[3].PT)

        self.histos["GLMASS"].Fill(Mass)
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
