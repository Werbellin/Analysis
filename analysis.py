#!/usr/bin/env python
# Based on P.Onyisi example
# -*- coding: utf-8 -*-
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

from ROOT import TCanvas, TH1D, TH2D, gROOT, gSystem, TFile, TTree, TLorentzVector, TChain, THStack, TColor
TH1D.SetDefaultSumw2()
gSystem.Load('libDelphes.so')
gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")
gROOT.SetBatch(True);
from cuts import *
from EventSelection import *
from step import *
if len(sys.argv) > 1:
    level_name = sys.argv[1]
    level = LEVELS.get(level_name, logging.NOTSET)
    logging.basicConfig(level=level)

logger = logging.getLogger('log')

logging.debug("Opening XML datafile")
#tree = ET.parse('SeperateFinalStates.xml')
#tree = ET.parse('test.xml')
#tree = ET.parse('LHE.xml')
tree = ET.parse('data.xml')
root = tree.getroot()

class Data :
    def __init__(self, name, xsection, files, data_type, category, total_events) :
        self.name               = name
        self.xsection           = xsection
        self.files              = files
        self.type               = data_type
        #self.identifier         = name + "(" + data_type + ")"
        self.totalNumberEvents  = total_events
        self.category           = category

RunData = [] # stores the data objects which include all the info about a dataset and the associated files

for dataset in root.findall('dataset') :
    xsection = float(dataset.find('xsection').text)
    category = dataset.get('category')
    name     = dataset.get('name')
    type     = dataset.get('type')
    #name     = name + "(" + type + ")"
    files    = root.findall(".//*[@name='" + name + "'][@type='" + type + "']/file")

    print "For dataset ", name, " with cross section of ", xsection, "fb the following files will be processed: "
    #print files

    filelist = []

    for file in files :
        print file.text
        filelist.append(file.text)

    t_chain = TChain("Delphes")

    for file in filelist :
        t_chain.Add(file[:])

    totalNumberEvents = t_chain.GetEntries()

    data = Data(name, xsection, filelist, type, category, totalNumberEvents)
    RunData.append(data)

print "Finished reading data.xml"

newStepList =[  Start("Total events processed"),
                #LeptonAcceptanceAnalysisCut("Delphes")]
                AllLeptonPtEtaPlot("Lepton plots"),
                JetMultiplicityPlot("JetMul"),
                LeptonDefinitionCut("Lepton Acceptance"),
                ZPairMassPlot("ZMass"),
                #TrueZCut("trueZ"),
                GoodLeptonPtEtaPlot("GoodLeptons"),
                ZPairDefinitionCut("Two Z candidates"),
                ZPairMassPlot("ZMass"),
                GoodLeptonPtEtaPlot("Z leptons"),
                LeptonIsolationPlot("R of Z Leptons"),
                DefineTaggingJetsCut("DefineTaggingJets"),
                TaggingJetKinematicsPlot("TJ kinematics"),
                TaggingJetZKinematicsPlot("TJZ kinematics"),
                ZKinematics("Kinematics of Z bosons"),
                ZJetsKinematics("Kinematics of TJ1 and L1/2"),
                ZeppenfeldVariablesPlot("YStar1"),
                TaggingJetRapidityGapCut("Tagging jet rapidity gap > 1.8"),
                #JetVetoCut("Jet veto"),
                TaggingJetInvariantMassCut("Invariant mass of tagging jets"),
                TaggingJetKinematicsPlot("TJ kinematics"),
                LeptonsBetweenTaggingJetsEtaCut("Leptons between tagging jets"),
                TaggingJetY1Y2Cut("Cut on y1y2"),
                JetVetoCut("Jet veto"),
                LeptonIsolationPlot("R of Z Leptons"),
                GoodLeptonPtEtaPlot("Leptons between jets")]


mEventSelection = EventSelection(RunData, newStepList)

for data in RunData:
    chain = TChain("Delphes")
    print data.name
    for file in data.files :
        print file
        chain.Add(file[:])

    mEventSelection.ChangeDataSet(data.name)
    print "Number of events in chain for dataset " + data.name +":  ", data.totalNumberEvents

    for eventNumber in range(data.totalNumberEvents) :
        logging.debug("Starting new event")
        if eventNumber % 500 == 0:
            update_progress(float(eventNumber)/data.totalNumberEvents)
            #print "Event=",eventNumber, "Fraction: ", eventNumber/chain.GetEntries()
        chain.GetEntry(eventNumber)
        #event = Event(chain)
        mEventSelection.ProcessEvent(chain, data.type)



mEventSelection.PrintCutflow()
mEventSelection.GetEfficencies()
mEventSelection.Finalize()
