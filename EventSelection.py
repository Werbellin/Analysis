#!/usr/local/bin/python
# -*- coding: utf-8 -*-
import logging, copy

LEVELS = {'debug': logging.DEBUG,
          'info': logging.INFO,
          'warning': logging.WARNING,
          'error': logging.ERROR,
          'critical': logging.CRITICAL}

from collections import Counter#, OrderedDict

from ROOT import TCanvas, TH1D, TH2D, gSystem, TFile, TTree, TLorentzVector, TChain, THStack, TColor, TText
TH1D.SetDefaultSumw2()
gSystem.Load('libDelphes.so')

from cuts import *
from Event import *
from upgrade import *
from step import *
log = logging.getLogger(__name__)
logging.getLogger(__name__).setLevel(logging.INFO)


class EventSelection :
    def __init__(self, dataset,run_name, step_list) :
        self.dataset = dataset
        self._RunName = run_name

        log.debug("EventSelection constructor called")
        self.ROOTFile = TFile("output/" + self._RunName + "/A-Result.root","RECREATE");
        self.Histos = {}
        self.currentDataset = None
        self.luminosity = 100.
        self._StepList = {}
        self._CutList = {}

        for data in dataset :
            self.Histos[data.name] = []
            newStepList = copy.deepcopy(step_list)
            cutflow = ""
            newCutList = []
            for step in newStepList :
                if isinstance(step, Cut) :
                    cutflow += step.CutAbbreviation + "_"
                    newCutList.append(step)
                step.Initialize(self.Histos[data.name], data.name, cutflow)

            self._StepList[data.name] = newStepList
            self._CutList[data.name] = newCutList
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
                log.debug("Event passed " + step.name)
            else :
                log.debug("Event failed cut " + step.name)
                return True

    def GetEfficencies(self) :
        width = 80
        str_double = "=" * width
        str_single = "-" * width
        log.info("%s", str_double)
        log.info("Cutflow efficency analysis")
        print str_single
        data = ("", 0.)
        signal = {}
        background = {}
        significance = []
        for step in self._StepList[self.currentDataset] :
            data = (step.name, 0.)
            signal[str(step)]= 0.
            background[str(step)] = 0.
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
            for step in self._CutList[data.name] :
                signal[str(step)] += step.NumberEventsPassedCut / data.totalNumberEvents * data.xsection

        for data in BackgroundData :
            for step in self._CutList[data.name] :
                background[str(step)] += step.NumberEventsPassedCut / data.totalNumberEvents * data.xsection

        for step in self._CutList[self.currentDataset] :
            s = signal[str(step)]
            b = math.sqrt(background[str(step)])
            if b <> 0. :
                significance.append( (str(step), s / b * math.sqrt(self.luminosity) ) )
            else:
                significance.append((str(step), -1))

        template = "{0:55}|{1:25}"
        print str_single
        print template.format("Cut name", "S/√B @ " + "{:3.2f}".format(self.luminosity) + "fb⁻¹")
        print str_single
        significanceFormated = []
        for tuple in significance :
            temp = (tuple[0], "{:4.2f}".format(tuple[1]))
            significanceFormated.append(temp)
        for tuple in significanceFormated :
            print template.format(*tuple)
        print str_double

    def PrintCutflow(self) :
        width = 100
        str_double = "=" * width
        str_single = "-" * width
        print str_double
        for data in self.dataset :
            print "Cutflow for the ", data.name, " dataset:"
            template = "{0:60}|{1:15}|{2:15}|{3:15}|{4:15}" # column widths: 8, 10, 15, 7, 10
            log.info("%s",str_single)
            print template.format("Cut name", "Passed","@"+ "{:3.1f}".format(self.luminosity) + "fb-1", "abs.diff", "efficency") # header
            print str_single
            for previous, item, nxt in previous_and_next(self._CutList[data.name]):
                diff     = 0.
                rel_diff = 0.
                eff      = 0.
                if previous is not None :
                    diff = - (previous.NumberEventsPassedCut - item.NumberEventsPassedCut)
                    if previous.NumberEventsPassedCut <> 0. :
                       rel_diff = diff/previous.NumberEventsPassedCut * 100.
                       eff = item.NumberEventsPassedCut/previous.NumberEventsPassedCut * 100.
                evtNumberMCFile = item.NumberEventsPassedCut
                evtNumberActual = item.NumberEventsPassedCut / data.totalNumberEvents * data.xsection * self.luminosity
                tuple = (item.__str__(), evtNumberMCFile, "{:10.3f}".format(evtNumberActual) , "{:10.2f}".format(diff) , "{:10.1f}".format(eff))
                print template.format(*tuple)
            print str_single


    def Finalize(self, plot_normalization = "LUMI") :
        selection = ""
        self.PlotNormalization = plot_normalization
        for data in self.dataset :
        #self.currentDataset = dataset_name
            data_dir = self.ROOTFile.mkdir(data.name)
            data_dir.cd()
            #print data.name, self.Histos
            for hist in self.Histos[data.name] :
                selection = hist.GetName()
                hist.Write()
                hist.SetLineColor(data.LineColor)
                hist.SetMarkerColor(data.LineColor)
                hist.SetLineWidth(3)
                print dir(hist)
                print "norm: ", getattr(hist, "_normalization")
                #if hist._normalization == "ONE" and hist.Integral() <> 0.:
                #    hist.Scale(1.0/hist.Integral())
                #if hist._drawOption <> "" :
                #    hist.Draw(hist._drawOption)
                #    canvas.Print("output//" + self._RunName + "//" + hist.GetName() + ".png")

        data_dir = self.ROOTFile.mkdir("Combined")
        data_dir.cd()
        stackdir = {}
        for hist in self.Histos[self.currentDataset] :
            stack = THStack(hist.GetName()[len(self.currentDataset):], hist.GetTitle())
            stackdir[hist.GetName()[len(self.currentDataset):]] = stack
        for data in self.dataset :
            for hist in self.Histos[data.name] :
                hist.SetTitle(data.name)
                stack = stackdir[hist.GetName()[len(data.name):]]
                if hist.Integral() <> 0. :
                    if self.PlotNormalization == "LUMI" :
                        hist.Scale( data.xsection * self.luminosity / data.totalNumberEvents )
                    if self.PlotNormalization == "PDF_WIDTH" :
                        hist.Scale(1.0/hist.Integral(), "width")
                stack.Add(hist)
                stack.Draw()
                stack.GetHistogram().GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
                stack.GetHistogram().GetYaxis().SetTitle(hist.GetYaxis().GetTitle())



        for stack in stackdir.itervalues() :
            canvas = TCanvas()
            stack.Draw("nostack")
            #canvas.Divide(1,2)
            #canvas.cd(1)
            #canvas.BuildLegend()
            #canvas.SetLogy()
            #if stack.GetName().find("LeptonMass") <> -1 :
                #print "found stack"
                #stack.GetXaxis().SetTitle("m_{4l}")
                #stack.GetYaxis().SetTitle("Events/20GeV")
            #stack.Draw("nostack")
            #stack.SetTitle("")
            #canvas.cd(2)
            text = TText()
            #text.DrawText(0, 0, selection)
            canvas.Print("output//" + self._RunName + "//" + stack.GetName() + ".png")
            stack.Write()



        self.ROOTFile.Close()
        return

