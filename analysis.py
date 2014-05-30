#!/usr/bin/env python
# Based on P.Onyisi example
import sys,os,string
import math
import numpy as np
#import OrderedSet

import xml.etree.ElementTree as ET

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

from collections import Counter#, OrderedDict

from ROOT import TCanvas, TH1D, gSystem, TFile, TTree, TLorentzVector, TChain, THStack, TColor
TH1D.SetDefaultSumw2()
gSystem.Load('libDelphes.so')

from cuts import *

if len(sys.argv) > 1:
    level_name = sys.argv[1]
    level = LEVELS.get(level_name, logging.NOTSET)
    logging.basicConfig(level=level)

logger = logging.getLogger('log')
class Event :
    def __init__(self, leaf) :
        self.data = leaf
        self.Cuts = {'start':True}
    

class EventSelection :
    def __init__(self, dataset) : # step_list) :
        self.dataset = dataset
        logger.debug("EventSelection constructor called")
        self.ROOTFile = TFile("EventSelection.root","RECREATE");
        self.cutflow = {}
        self.Histos = {}
        self.currentDataset = None
        self.luminosity = 100.




        for data in dataset :
            self.Histos[data.name] =[]
            #newStepList = copy.deepcopy(step_list)
            #for step in newSte :
            #    step.Initialize(self.Histos[data.name], data.name)

            #self._StepList[data.name] = newStepList
            
            cut0 = Start("Total events processed", self.Histos[data.name], data.name)
            cut1 = LeptonAcceptance("Lepton Acceptance Eta PT", self.Histos[data.name], data.name)
            cut2 = ZPair("Two Z candidates", self.Histos[data.name], data.name)
            cut3 = DefineTaggingJets("DefineTaggingJets", self.Histos[data.name], data.name)
            cut4 = ZKinematics("Kinematics of Z bosons", self.Histos[data.name],data.name)
            cut5 = ZJetsKinematics("Kinematics of TJ1 and L1/2", self.Histos[data.name], data.name)
            cut6 = ZeppenfeldVar("Zeppenfeld variables", self.Histos[data.name], data.name)
            cut7 = TaggingJetInvariantMass("Invariant mass of tagging jets", self.Histos[data.name], data.name)
            cut8 = LeptonIsolation("Lepton Isolation", self.Histos[data.name], data.name)
            newcutflow = [cut0, cut1, cut2, cut3, cut4, cut5, cut6, cut8, cut7]
            self.cutflow[data.name] = newcutflow
     
            #print self.cutflow

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
        for cut in self.cutflow[self.currentDataset] :
            if cut.ApplyCut(self.event, self.Histos[self.currentDataset], data_type) :
                logger.debug("Event passed " + cut.name)
            else :
                logger.debug("Event failed cut " + cut.name)
                return True
        

    def GetEfficencies(self) :
        print "Efficiency of cuts on Background suppression"

        data = ("", 0.)
        signal = {}
        background = {}
        significance = []
        for step in self.cutflow[self.currentDataset] :
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
            for step in self.cutflow[data.name] :
                #print "Sig ", step.name, step.NumberEventsPassedCut / data.totalNumberEvents
                #for k,v in signal :
                 #   if k == step.name :
                  #      v += step.NumberEventsPassedCut / data.totalNumberEvents * data.xsection
                
                signal[step.name] += step.NumberEventsPassedCut / data.totalNumberEvents * data.xsection
                    
        for data in BackgroundData :
            for step in self.cutflow[data.name] :
                #for k,v in background :
                 #   if k == step.name :
                #        v += step.NumberEventsPassedCut / data.totalNumberEvents * data.xsection               

                #item for item in a if item[0] == step.number
                #print "BG ", step.name, step.NumberEventsPassedCut / data.totalNumberEvents
                background[step.name] += step.NumberEventsPassedCut / data.totalNumberEvents * data.xsection

        for step in self.cutflow[self.currentDataset] :
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

            for previous, item, nxt in previous_and_next(self.cutflow[data.name]):
                if True : #item.IsCut() :
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
                if hist.Integral() <> 0. :
                    hist.Scale(1.0/hist.Integral())#, "width")
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

logging.debug("Opening XML datafile")
tree = ET.parse('data.xml')
root = tree.getroot()

#for child in root :
#    print child.tag, child.attrib

class Data :
    def __init__(self, name, xsection, files, data_type, category, total_events) :
        self.name = name
        self.xsection = xsection
        self.files = files
        self.type = data_type
        self.totalNumberEvents = total_events
        self.category = category

RunData = []
for dataset in root.findall('dataset') :
    xsection = float(dataset.find('xsection').text)
    category = dataset.get('category')
    name      = dataset.get('name')
    type = dataset.get('type')
    print "For dataset ", name, " with cross section of ", xsection, "fb the following files will be processed: "
    files = root.findall(".//*[@name='" + name + "']/file")
    #print files
    filelist =[]
    for file in files :
        print file.text
        filelist.append(file.text)
    t_chain = TChain("Delphes")
    for file in filelist :
        t_chain.Add(file[:])

    totalNumberEvents = t_chain.GetEntries()

    d = Data(name, xsection, filelist, type, category, totalNumberEvents)
    RunData.append(d)

print "Finished reading data.xml"
mEventSelection = EventSelection(RunData)
for data in RunData:
    chain = TChain("Delphes")
    print data.name
    for file in data.files :
        print file
        chain.Add(file[:])

    mEventSelection.ChangeDataSet(data.name)
    print "Number of events in chain for dataset " + data.name +":  ", chain.GetEntries()
 
    for eventNumber in range(chain.GetEntries()) : #chain:   # loop over all event      
        logging.debug("Starting new event")                                                                                                     
        if eventNumber%1000==0: 
            print "Event=",eventNumber                                                                                                  
        chain.GetEntry(eventNumber)
        event = Event(chain)     
        mEventSelection.ProcessEvent(chain, data.type)

    #mEventSelection.SaveHistos()
mEventSelection.PrintCutflow()
mEventSelection.GetEfficencies()
mEventSelection.Finalize()
