class Event :
    def __init__(self, leaf) :
        self.data = leaf
        self.Cuts = {'start':True}
        self.CutHistory = ""
        for event in self.data.Event :
            self.weight = event.Weight

    def Calculate(self, qty_name) :
        if qty_name == "ZZSPT" :
            self.ZZSPT = (self.Z1Particles[0] + self.Z1Particles[1]).Pt() + (self.Z2Particles[0] + self.Z2Particles[1]).Pt()
            return self.ZZSPT

        if qty_name == "ZZPT" :
            self.ZZPT = (self.Z1Particles[0] + self.Z1Particles[1] + self.Z2Particles[0] + self.Z2Particles[1]).Pt()
            return self.ZZPT

        if qty_name == "ZZY" :
            self.ZZY = (self.Z1Particles[0] + self.Z1Particles[1] + self.Z2Particles[0] + self.Z2Particles[1]).Rapidity()
            return self.ZZY

        if qty_name == "TJDEta" :
            self.TJDEta = self.TaggingJet1.Eta() - self.TaggingJet2.Eta()
            return self.TJDEta


        if qty_name == "TJDY" :
            self.TJDY = self.TaggingJet1.Rapidity() - self.TaggingJet2.Rapidity()
            return self.TJDY

        if qty_name == "TJM" :
            self.TJM = (self.TaggingJet1 + self.TaggingJet2).M()
            return self.TJM

        if qty_name == "YStar1" :
            self.YStar1 = self.Z1.Rapidity() - 0.5 * (self.TaggingJet1.Rapidity() + self.TaggingJet2.Rapidity())
            return self.YStar1

        if qty_name == "YStar2" :
            self.YStar2 = self.Z2.Rapidity() - 0.5 * (self.TaggingJet1.Rapidity() + self.TaggingJet2.Rapidity())
            return self.YStar2


        if qty_name == "TJM" :
            self.TJM = (self.TaggingJet1 + self.TaggingJet2).M()
            return self.TJM

        if qty_name == "TJDPhi" :
            self.TJDPhi = self.TaggingJet1.Phi() - self.TaggingJet2.Phi()
            return self.TJDPhi

        if qty_name == "YStarZZ" :
            self.YStarZZ = self.getQuantity("ZZY") - 0.5 * (self.TaggingJet1.Rapidity() + self.TaggingJet2.Rapidity())
            return self.YStarZZ

        else :
            return None

    def getQuantity(self, qty_name) :
        if hasattr(self, qty_name) :
            return getattr(self, qty_name)
        else :
            return self.Calculate(qty_name)
