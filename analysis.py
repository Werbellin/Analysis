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

from  time import strftime
RunName = strftime("%Y-%m-%d %H:%M:%S")
#if not os.path.exists(RunName): 
os.makedirs("output//" + RunName)

class Tee(object):
     def __init__(self, *files):
         self.files = files
     def write(self, obj):
         for f in self.files:
             f.write(obj)

f = open('out.txt', 'w')
original = sys.stdout
#sys.stdout = Tee(sys.stdout, f)

logging.debug("Opening XML datafile")
#tree = ET.parse('SeperateFinalStates.xml')
#tree = ET.parse('NPAC_talk.xml')
#tree = ET.parse('eff.xml')
tree = ET.parse('data_Hv1.xml')
root = tree.getroot()

class Data :
    def __init__(self, name, xsection, files, data_type, category, total_events, line_color) :
        self.name               = name
        self.xsection           = xsection
        self.LineColor          = line_color
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
    line_color = int(dataset.find('color').text)
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

    data = Data(name, xsection, filelist, type, category, totalNumberEvents, line_color)
    RunData.append(data)

print "Finished reading data.xml"

newStepList =[  Start("Total events processed"),
                LeptonTriggerCut("Trigger RunI"),
                LeptonTriggerCut("Trigger RunII", l_leading_pt = 22., l_subleading_pt = 10.),
                #LeptonAcceptanceAnalysisCut("Delphes", el_pt = 6., el_eta = 2.5, mu_pt = 6., mu_eta = 2.5)]
                AllLeptonPtEtaPlot("Lepton plots"),
                JetMultiplicityPlot("JetMul"),
                LeptonDefinitionCut("Lepton Acceptance", el_eta = 3.0),
                GoodLeptonPtEtaPlot("GoodLeptons"),
                ZPairMassPlot("ZMass"),
                #TrueZCut("trueZ"),
                ZPairDefinitionCut("Two Z candidates"),
                ZPairMassPlot("ZMass"),
                ZZRapidityPlot("ZZ system boost"),
                #ZZRapidityCut("ZZ y cut"),
                Z1LeptonPtPlot("Z1 Lepton pT"),
                GoodLeptonPtEtaPlot("Z leptons"),
                LeptonIsolationPlot("R of Z Leptons"),
                DefineTaggingJetsCut("DefineTaggingJets", jet_pt = 30.0, jet_eta = 4.7),
                YStarZZPlot("y_ZZ"),
                YEtaPlot("yeta"),
                TaggingJetKinematicsPlot("TJ kinematics"),
                TaggingJetZKinematicsPlot("TJZ kinematics"),
                ZKinematics("Kinematics of Z bosons"),
                ZJetsKinematics("Kinematics of TJ1 and L1/2"),
                ZeppenfeldVariablesPlot("YStar1"),
                TaggingJetRapidityGapCut("Tagging jet rapidity gap 2"),
                ZZRapidityPlot("ZZ system boost"),
                YStarZZPlot("y_ZZ"),
                GoodLeptonPtEtaPlot("Z leptons 2"),
                #TaggingJetKinematicsPlot("TJ kinematics 2"),
                #JetVetoCut("Jet veto"),
                TaggingJetKinematicsPlot("TJ kinematics 3"),
                TaggingJetInvariantMassCut("Invariant mass of tagging jets"),
                TaggingJetKinematicsPlot("TJ kinematics 2"),
                TaggingJetInvariantMassCut("Invariant mass of tagging jets 2", mass_cut = 600.),
                TaggingJetInvariantMassCut("Invariant mass of tagging jets 3", mass_cut = 700.),
                ZeppenfeldVariablesPlot("YStar1"),
                TaggingJetKinematicsPlot("TJ kinematics"),
                GoodLeptonPtEtaPlot("Z leptons 3"),
                #TaggingJetRapidityGapCut("Tagging jet rapidity gap"),
                LeptonsBetweenTaggingJetsEtaCut("Leptons between tagging jets"),
                #TaggingJetY1Y2Cut("Cut on y1y2"),
                GoodLeptonPtEtaPlot("Z leptons"),
                JetMultiplicityPlot("Jet multiplicity"),
                JetVetoCut("Jet veto"),
                JetVetoCut("Jet veto 2", pt_cut = 20.),
                LeptonIsolationPlot("R of Z Leptons"),
                YStarZZPlot("y_ZZ"),
                ZZRapidityPlot("ZZ system boost"),
                GoodLeptonPtEtaPlot("Leptons between jets")]


mEventSelection = EventSelection(RunData,RunName, newStepList)

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
mEventSelection.Finalize(plot_normalization = "PDF_WIDTH")
