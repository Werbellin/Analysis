import logging, copy

LEVELS = {'debug': logging.DEBUG,
          'info': logging.INFO,
          'warning': logging.WARNING,
          'error': logging.ERROR,
          'critical': logging.CRITICAL}

from collections import Counter#, OrderedDict

from ROOT import TCanvas, TH1D, TH2D, gSystem, TFile, TTree, TLorentzVector, TChain, THStack, TColor
TH1D.SetDefaultSumw2()
gSystem.Load('libDelphes.so')

from cuts import *
from Event import *
from upgrade import *
logger = logging.getLogger('log')


class EventSelection :
    def __init__(self, dataset, step_list) :
        self.dataset = dataset
        logger.debug("EventSelection constructor called")
        self.ROOTFile = TFile("EventSelection.root","RECREATE");
        self.Histos = {}
        self.currentDataset = None
        self.luminosity = 100.
        self._StepList = {}
        for data in dataset :
            self.Histos[data.name] = []
            newStepList = copy.deepcopy(step_list)
            for step in newStepList :
                step.Initialize(self.Histos[data.name], data.name)

            self._StepList[data.name] = newStepList

    def ChangeDataSet(self, dataset_name) :
        self.currentDataset = dataset_name
        #new_dir = self.ROOTFile.mkdir(dataset_name);
        #new_dir.cd() 

    def SaveHistos(self) :
        for data in self.dataset :
            for hist in self.Histos :
                hist.Write()

    def ProcessEvent(self, leaf, data_type) :
        self.event = Event(leaf)
        self.BeginSelection(data_type)

    def BeginSelection(self, data_type) :
        for step in self._StepList[self.currentDataset] :
            if step.PerformStep(self.event, self.Histos[self.currentDataset], data_type) :
                logger.debug("Event passed " + step.name)
            else :
                logger.debug("Event failed cut " + step.name)
                return True

    def GetEfficencies(self) :
        print "Efficiency of cuts on Background suppression"

        data = ("", 0.)
        signal = {}
        background = {}
        significance = []
        for step in self._StepList[self.currentDataset] :
            data = (step.name, 0.)
            signal[step.name]= 0.
            background[step.name] = 0.
            #significance.append(data)

        SignalData = []
        BackgroundData = []
        for data in self.dataset :
            if data.category == "Signal" :
                SignalData.append(data)
            if data.category == "Background" :
                BackgroundData.append(data)

        print "The following datasets are considered as signal:"
        for item in SignalData :
            print item.name
        print "The following datasets are considered as backgrounds: "
        for item in BackgroundData :
            print item.name

        for data in SignalData :
            for step in self._StepList[data.name] :
                #print "Sig ", step.name, step.NumberEventsPassedCut / data.totalNumberEvents
                #for k,v in signal :
                 #   if k == step.name :
                  #      v += step.NumberEventsPassedCut / data.totalNumberEvents * data.xsection
                
                signal[step.name] += step.NumberEventsPassedCut / data.totalNumberEvents * data.xsection
                    
        for data in BackgroundData :
            for step in self._StepList[data.name] :
                #for k,v in background :
                 #   if k == step.name :
                #        v += step.NumberEventsPassedCut / data.totalNumberEvents * data.xsection               

                #item for item in a if item[0] == step.number
                #print "BG ", step.name, step.NumberEventsPassedCut / data.totalNumberEvents
                background[step.name] += step.NumberEventsPassedCut / data.totalNumberEvents * data.xsection

        for step in self._StepList[self.currentDataset] :
            s = signal[step.name]
            b = math.sqrt(background[step.name]) * math.sqrt(self.luminosity)
            if b <> 0. :
                significance.append( (step.name, s / b))
            else:
                significance.append((step.name, -1))
            
        for k in significance :
            print k

    def PrintCutflow(self) :    
        for data in self.dataset :
            print "Cutflow for the ", data.name, " dataset:"
            print "========================================================="
            template = "{0:40}|{1:10}|{2:10}|{3:10}" # column widths: 8, 10, 15, 7, 10
            print template.format("Cut name", "Passed", "abs.diff", "rel. diff") # header

            for previous, item, nxt in previous_and_next(self._StepList[data.name]):
                if isinstance(item, Cut) :
                    diff = 0.
                    rel_diff = 0.
                    if previous is not None : 
                        diff = - (previous.NumberEventsPassedCut - item.NumberEventsPassedCut)
                        if previous.NumberEventsPassedCut <> 0. : 
                            rel_diff = diff/previous.NumberEventsPassedCut * 100.

                    tuple = (item.name, item.NumberEventsPassedCut, "{:10.2f}".format(diff) , "{:10.1f}".format(rel_diff))
                    print template.format(*tuple)
        print "========================================================="

    def Finalize(self) :    
        for data in self.dataset :
        #self.currentDataset = dataset_name
            data_dir = self.ROOTFile.mkdir(data.name)
            data_dir.cd()
            #print data.name, self.Histos
            for hist in self.Histos[data.name] :
                hist.Write()

        data_dir = self.ROOTFile.mkdir("Combined")
        data_dir.cd()
        stackdir = {}
        for hist in self.Histos[self.currentDataset] :
            stack = THStack(hist.GetName()[len(self.currentDataset):], hist.GetName()[len(self.currentDataset):])
            stackdir[hist.GetName()[len(self.currentDataset):]] = stack
        for data in self.dataset :
            for hist in self.Histos[data.name] :
                hist.SetTitle(data.name)
                if data.name.find("QCD")<> -1 :
                    hist.SetLineColor(2)
                if data.name.find("EWK")<> -1 :
                    hist.SetLineColor(3)
                if data.name.find("VBF")<> -1 :
                    hist.SetLineColor(6)
                #print "data.name: ", data.name
                #print "Trying to retrieve",hist.GetName()[len(data.name):]
                #print "Content of stackdir: ", stackdir
                stack = stackdir[hist.GetName()[len(data.name):]]
                #if hist.Integral() <> 0. :
                #    hist.Scale(1.0/hist.Integral())#, "width")
                stack.Add(hist)


        for stack in stackdir.itervalues() :
            canvas = TCanvas()
            stack.Draw("nostack")
            canvas.BuildLegend()
            #stack.Draw("nostack")
            canvas.SaveAs(stack.GetName() + ".pdf")
            
            stack.Write()



        self.ROOTFile.Close()
        return

