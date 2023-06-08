#!/usr/bin/env python
# -*-coding:Latin-1 -*
import ROOT
import ROOT as root
import matplotlib.pyplot as plt
import atlasplots as aplt
import numpy as np
import sys
import os
import time
import gc

from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TH1
from source_plotting.source  		import *
from source_plotting.source_Xs  	import *
from source_plotting.source_Optim  	import *
from source_plotting.source_diff  	import *

from math import *


""" ****************************************************************** Define the input ********************************************************************* """

MC_directory = "/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[3]+"_MC_"+sys.argv[6]+"TeV/Nominal/"
Directory    = os.listdir(MC_directory)
for i in Directory:
	if "Sherpa_225" in i:
		MC_Sherpa = ROOT.TFile.Open(MC_directory+i)
	if "PhH7EG_225" in i:
		MC_PhH7EG = ROOT.TFile.Open(MC_directory+i)

print(MC_Sherpa)
print(MC_PhH7EG)
pT_Sherpa  = MC_Sherpa.Get("TruthSelection/Lepton_pT_Truth_cut4")
Eta_Sherpa = MC_Sherpa.Get("TruthSelection/Lepton_Eta_Truth_cut4")
pT_PhH7EG  = MC_PhH7EG.Get("TruthSelection/Lepton_pT_Truth_cut4")
Eta_PhH7EG = MC_PhH7EG.Get("TruthSelection/Lepton_Eta_Truth_cut4")

path 	   = "/sps/atlas/h/hatmani/W_xs/RooUnfold/Output/output_"
Summarize  = ROOT.TFile.Open(path+sys.argv[2]+str(sys.argv[6])+"/"+sys.argv[4]+"/Summarize_"+sys.argv[2]+str(sys.argv[6])+".root")
Bias       = ROOT.TFile.Open(path+sys.argv[2]+str(sys.argv[6])+"/"+sys.argv[4]+"/Bias_"		+sys.argv[2]+str(sys.argv[6])+".root")

IsoSF      = ROOT.TFile.Open(path+sys.argv[2]+str(sys.argv[6])+"/"+sys.argv[4]+"/Syst_IsoSys.root")
RecoSF     = ROOT.TFile.Open(path+sys.argv[2]+str(sys.argv[6])+"/"+sys.argv[4]+"/Syst_RecoSys.root")
TrigSF     = ROOT.TFile.Open(path+sys.argv[2]+str(sys.argv[6])+"/"+sys.argv[4]+"/Syst_TrigSys.root")
Recoil     = ROOT.TFile.Open(path+sys.argv[2]+str(sys.argv[6])+"/"+sys.argv[4]+"/Recoil_Syst.root")
Calib      = ROOT.TFile.Open(path+sys.argv[2]+str(sys.argv[6])+"/"+sys.argv[4]+"/Calib_Syst.root")

if "enu" in sys.argv[2]:
	IdSF       = ROOT.TFile.Open(path+sys.argv[2]+str(sys.argv[6])+"/"+sys.argv[4]+"/Syst_ElIDSys.root")

if "munu" in sys.argv[2]:
	IdSF       = ROOT.TFile.Open(path+sys.argv[2]+str(sys.argv[6])+"/"+sys.argv[4]+"/Syst_MuTTVASys.root")


NIter = 5

""" ********************************************************************* Nominal Plot ********************************************************************** """
"""
GetAcceptanceFactors(		Summarize, 							sys.argv[2]+sys.argv[6], sys.argv[4], sys.argv[8])
ShowNominalDistribution(	Summarize, 							sys.argv[2]+sys.argv[6], sys.argv[4], sys.argv[8])
CompareBias(				Summarize, 	Bias, 	   1,  NIter,   sys.argv[2]+sys.argv[6], sys.argv[4], sys.argv[8])
CompareStatError(			Summarize, 			   1,  NIter,   sys.argv[2]+sys.argv[6], sys.argv[4], sys.argv[8])
BiasProcedure( 				Summarize, 	Bias,    			sys.argv[2]+sys.argv[6], sys.argv[4], sys.argv[8])
"""
""" ********************************************************************* Syst Plot ********************************************************************** """
"""
if "enu" in sys.argv[2]:
	CompareSystId(    IdSF,   1, 		NIter, 	sys.argv[2]+sys.argv[6],   sys.argv[4], sys.argv[8])
	CompareSystCalib( Calib, 1, 		NIter, 	sys.argv[2]+sys.argv[6],   sys.argv[4], sys.argv[8])

CompareSystIso(  	  IsoSF,  1, 		NIter, 	sys.argv[2]+sys.argv[6],   sys.argv[4], sys.argv[8])
CompareSystReco( 	  RecoSF, 1, 		NIter, 	sys.argv[2]+sys.argv[6],   sys.argv[4], sys.argv[8])
CompareSystTrig( 	  TrigSF, 1, 		NIter, 	sys.argv[2]+sys.argv[6],   sys.argv[4], sys.argv[8])
CompareSystRecoil(   Recoil, 1, 		NIter, 	sys.argv[2]+sys.argv[6],   sys.argv[4], sys.argv[8])

if "enu" in sys.argv[2]:
	CompareSystEl( Summarize, Bias, IdSF, IsoSF, RecoSF, TrigSF, Recoil, Calib, 1, NIter, 	sys.argv[2]+sys.argv[6],   sys.argv[4], sys.argv[8])
if "munu" in sys.argv[2]:
	CompareSystMu( Summarize, Bias, IsoSF, RecoSF, TrigSF, Recoil, 1, NIter, 	sys.argv[2]+sys.argv[6],   sys.argv[4], sys.argv[8])
"""
	
""" *************************************************************** fiducial Xs Plot ****************************************************************** """

#GetSummaringTableEl(		Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, sys.argv[1], sys.argv[8], sys.argv[2]+sys.argv[6], float(sys.argv[7]), sys.argv[4])
#GetFiducialXsEl(			Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, sys.argv[1], sys.argv[8], sys.argv[2]+sys.argv[6], float(sys.argv[7]), sys.argv[4])

""" *************************************************************** differential Xs Plot ****************************************************************** """

#if "enu" in sys.argv[2]:
#	GetDiffernetialXsEl(			Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, sys.argv[1], sys.argv[8], sys.argv[2]+sys.argv[6], float(sys.argv[7]), sys.argv[4])
#if "munu" in sys.argv[2]:
#	GetDiffernetialXsMu(			Summarize, Bias, TrigSF, RecoSF, IsoSF,	IdSF, Calib, Recoil, sys.argv[1], sys.argv[8], sys.argv[2]+sys.argv[6], float(sys.argv[7]), sys.argv[4])

#GetDiffernetialXsPlot(				Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, sys.argv[1], sys.argv[8], sys.argv[2]+sys.argv[6], float(sys.argv[7]), sys.argv[4], pT_Sherpa, Eta_Sherpa, pT_PhH7EG, Eta_PhH7EG)
GetDiffernetialXsPlotNormalised(	Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, sys.argv[1], sys.argv[8], sys.argv[2]+sys.argv[6], float(sys.argv[7]), sys.argv[4], pT_Sherpa, Eta_Sherpa, pT_PhH7EG, Eta_PhH7EG)





#if "enu" in sys.argv[2]:
#	GetDiffernetialXsEl(		Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, sys.argv[1], sys.argv[8], sys.argv[2]+sys.argv[6], float(sys.argv[7]), sys.argv[4])
#	GetDiffernetialXsPlot(		Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, sys.argv[1], sys.argv[8], sys.argv[2]+sys.argv[6], float(sys.argv[7]), sys.argv[4], pT_Sherpa, Eta_Sherpa, pT_PhH7EG, Eta_PhH7EG)

#if "munu" in sys.argv[2]:
#	GetDiffernetialXsMu(		Summarize, Bias, TrigSF, RecoSF, IsoSF,	Recoil, sys.argv[1], sys.argv[8], sys.argv[2]+sys.argv[6], float(sys.argv[7]), sys.argv[4])
#	GetDiffernetialXsPlot(		Summarize, Bias, TrigSF, RecoSF, IsoSF, TrigSF, Recoil, Recoil, sys.argv[1], sys.argv[8], sys.argv[2]+sys.argv[6], float(sys.argv[7]), sys.argv[4], pT_Sherpa, Eta_Sherpa, pT_PhH7EG, Eta_PhH7EG)


#DiffCrossSectionX.GetDiffernetialXsComp(Summarize_plusenu5, Summarize_plusmunu5, Summarize_minusenu5, Summarize_minusmunu5, list1, "5TeV", "W$^{+}$ $\\rightarrow$ e$^{+} \\nu $, 5TeV, Uncertainties in (\%)", "Wplusenu5", 256.827)
#GetDiffernetialXs(Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, sys.argv[1], sys.argv[8], sys.argv[2]+sys.argv[6], int(sys.argv[7])
#GetDiffernetialXsPlot(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, "5TeV", "W^{+}#rightarrow e^{+}#nu, 5TeV", "Wplusenu5", 256.827)
#CrossSectionDeter.GetDiffernetialXsPlotN(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, "5TeV", "W^{+}#rightarrow e^{+}#nu, 5TeV", "Wplusenu5", 256.827)

print("end")
os._exit(0)



""" ********************************************************************************************************************************************************* """
""" ****************************************************************** Define the input ********************************************************************* """
""" ********************************************************************************************************************************************************* """
"""
# plus enu 5 TeV
MCsamples_plusenu5  = ROOT.TFile.Open("/eos/user/h/hatmani/PostDoc/W_Xs/DataSets_5TeV/WCrossSections_wplusenu_MC_5TeV/Nominal/merge/MC_5TeV_wplusenu.root")
Summarize_plusenu5  = ROOT.TFile.Open("/afs/cern.ch/work/h/hatmani/PostDoc/W_Xs/RooUnfold/Output/output_Wplusenu5/Lepton_pT/Summarize_Wplusenu5.root")
Bias_plusenu5       = ROOT.TFile.Open("/afs/cern.ch/work/h/hatmani/PostDoc/W_Xs/RooUnfold/Output/output_Wplusenu5/Lepton_pT/Bias_Wplusenu5.root")
#IdSF_plusenu5       = ROOT.TFile.Open("/afs/cern.ch/work/h/hatmani/HistMaker_v1/Plotting/RooUnfold/output_Wplusenu5/elEtaSF/Syst_ElIDSys.root")
#IsoSF_plusenu5      = ROOT.TFile.Open("/afs/cern.ch/work/h/hatmani/HistMaker_v1/Plotting/RooUnfold/output_Wplusenu5/elEtaSF/Syst_ElIsoSys.root")
#RecoSF_plusenu5     = ROOT.TFile.Open("/afs/cern.ch/work/h/hatmani/HistMaker_v1/Plotting/RooUnfold/output_Wplusenu5/elEtaSF/Syst_ElRecoSys.root")
#TrigSF_plusenu5     = ROOT.TFile.Open("/afs/cern.ch/work/h/hatmani/HistMaker_v1/Plotting/RooUnfold/output_Wplusenu5/elEtaSF/Syst_ElTrigSys.root")
#Recoil_plusenu5     = ROOT.TFile.Open("/afs/cern.ch/work/h/hatmani/HistMaker_v1/Plotting/RooUnfold/output_Wplusenu5/elEtaSF/Recoil_Syst.root")
#Calib_plusenu5      = ROOT.TFile.Open("/afs/cern.ch/work/h/hatmani/HistMaker_v1/Plotting/RooUnfold/output_Wplusenu5/elEtaSF/Calib_Syst.root")
"""

#list1 = []
#list1.append(IdSF_plusenu5)
#list1.append(IsoSF_plusenu5)
#list1.append(RecoSF_plusenu5)
#list1.append(TrigSF_plusenu5)
#list1.append(Recoil_plusenu5)
#list1.append(Calib_plusenu5)
#list1.append(Bias_plusenu5) 

""" ********************************************************************************************************************************************************* """
""" ****************************************************************** Get Nominal Plot ********************************************************************* """
""" ********************************************************************************************************************************************************* """

''' Wplusenu 5TeV '''
"""
NominalPlots.'GetEpsilonFactors'(Summarize_plusenu5,                "Wplusenu5",   "W^{+}#rightarrow e^{+}#nu,   5TeV")
NominalPlots.GetAcceptanceFactors(Summarize_plusenu5,             	"Wplusenu5",   "W^{+}#rightarrow e^{+}#nu,   5TeV")
NominalPlots.MigrationMatrix(Summarize_plusenu5,                  	"Wplusenu5",   "W^{+}#rightarrow e^{+}#nu,   5TeV")
NominalPlots.ShowNominalDistribution(Summarize_plusenu5,          	"Wplusenu5",   "W^{+}#rightarrow e^{+}#nu,   5TeV")
NominalPlots.CompareBias(Bias_plusenu5, 1,  9,                   	"Wplusenu5",   "W^{+}#rightarrow e^{+}#nu,   5TeV")
NominalPlots.CompareStatError(Summarize_plusenu5, 1,  9,         	"Wplusenu5",   "W^{+}#rightarrow e^{+}#nu,   5TeV")
NominalPlots.BiasProcedure( Summarize_plusenu5, Bias_plusenu5,    	"Wplusenu5",   "W^{+}#rightarrow e^{+}#nu,   5TeV" )
"""
""" ********************************************************************************************************************************************************* """
""" ****************************************************************** Get System Plots ********************************************************************* """
""" ********************************************************************************************************************************************************* """

''' Wplusenu 5TeV '''

#SystematicsStudy.CompareSystId(   IdSF_plusenu5,   1, 10, "Wplusenu5",   "W^{+}#rightarrow e^{+}#nu,   5TeV")
#SystematicsStudy.CompareSystIso(  IsoSF_plusenu5,  1, 10, "Wplusenu5",   "W^{+}#rightarrow e^{+}#nu,   5TeV")
#SystematicsStudy.CompareSystReco( RecoSF_plusenu5, 1, 10, "Wplusenu5",   "W^{+}#rightarrow e^{+}#nu,   5TeV")
#SystematicsStudy.CompareSystTrig( TrigSF_plusenu5, 1, 10, "Wplusenu5",   "W^{+}#rightarrow e^{+}#nu,   5TeV")

#SystematicsStudy.CompareSystRecoil(  Summarize_plusenu5, Recoil_plusenu5, 1, 10, "Wplusenu5",   "W^{+}#rightarrow e^{+}#nu,   5TeV")
#SystematicsStudy.CompareSystCalib(  Summarize_plusenu5,  Calib_plusenu5,  1, 10, "Wplusenu5",   "W^{+}#rightarrow e^{+}#nu,   5TeV")
#SystematicsStudy.CompareSyst( Summarize_plusenu5, IdSF_plusenu5, IsoSF_plusenu5, RecoSF_plusenu5, TrigSF_plusenu5, Recoil_plusenu5, Calib_plusenu5, "Wplusenu5", "W^{+}#rightarrow e^{+}#nu,   5TeV")



""" ********************************************************************************************************************************************************* """
""" ************************************************************** Optimisation Study 1 bin ***************************************************************** """
""" ********************************************************************************************************************************************************* """

#Optimisation1B.StatStudy(Summarize_plusenu5, 1, 20, 12, "Wplusenu5", "$p^{T}_{W}$")    # define the number of iterations for the study
#Optimisation1B.StatStudy(Summarize_plusenu5, 1, 20, 13, "Wplusenu5", "$p^{T}_{W}$")    # define the number of iterations for the study
#Optimisation1B.StatStudy(Summarize_plusenu5, 1, 20, 14, "Wplusenu5", "$p^{T}_{W}$")    # define the number of iterations for the study
#Optimisation1B.StatStudy(Summarize_plusenu5, 1, 20, 15, "Wplusenu5", "$p^{T}_{W}$")    # define the number of iterations for the study
#Optimisation1B.StatStudy(Summarize_plusenu5, 1, 20, 16, "Wplusenu5", "$p^{T}_{W}$")    # define the number of iterations for the study
#Optimisation1B.StatStudy(Summarize_plusenu5, 1, 20, 17, "Wplusenu5", "$p^{T}_{W}$")    # define the number of iterations for the study
#Optimisation1B.StatStudy(Summarize_plusenu5, 1, 20, 20, "Wplusenu5", "$p^{T}_{W}$")    # define the number of iterations for the study
#Optimisation1B.StatStudy(Summarize_plusenu5, 1, 20, 21, "Wplusenu5", "$p^{T}_{W}$")    # define the number of iterations for the study
#Optimisation1B.StatStudy(Summarize_plusenu5, 1, 20, 22, "Wplusenu5", "$p^{T}_{W}$")    # define the number of iterations for the study

#Optimisation1B.BiasStudy(Summarize_plusenu5, Bias_plusenu5, 1, 20,  12,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.BiasStudy(Summarize_plusenu5, Bias_plusenu5, 1, 20,  13,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.BiasStudy(Summarize_plusenu5, Bias_plusenu5, 1, 20,  14,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.BiasStudy(Summarize_plusenu5, Bias_plusenu5, 1, 20,  15,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.BiasStudy(Summarize_plusenu5, Bias_plusenu5, 1, 20,  16,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.BiasStudy(Summarize_plusenu5, Bias_plusenu5, 1, 20,  17,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.BiasStudy(Summarize_plusenu5, Bias_plusenu5, 1, 20,  20,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.BiasStudy(Summarize_plusenu5, Bias_plusenu5, 1, 20,  21,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.BiasStudy(Summarize_plusenu5, Bias_plusenu5, 1, 20,  22,  "Wplusenu5", "$p^{T}_{W}$")

#Optimisation1B.EffSystematicStudy(Summarize_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, 1, 19, 12,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.EffSystematicStudy(Summarize_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, 1, 19, 13,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.EffSystematicStudy(Summarize_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, 1, 19, 14,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.EffSystematicStudy(Summarize_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, 1, 19, 15,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.EffSystematicStudy(Summarize_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, 1, 19, 16,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.EffSystematicStudy(Summarize_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, 1, 19, 17,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.EffSystematicStudy(Summarize_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, 1, 19, 20,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.EffSystematicStudy(Summarize_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, 1, 19, 21,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.EffSystematicStudy(Summarize_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, 1, 19, 22,  "Wplusenu5", "$p^{T}_{W}$")

#Optimisation1B.CalibSystematicStudy(Summarize_plusenu5, Calib_plusenu5, 1, 19, 12,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.CalibSystematicStudy(Summarize_plusenu5, Calib_plusenu5, 1, 19, 13,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.CalibSystematicStudy(Summarize_plusenu5, Calib_plusenu5, 1, 19, 14,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.CalibSystematicStudy(Summarize_plusenu5, Calib_plusenu5, 1, 19, 15,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.CalibSystematicStudy(Summarize_plusenu5, Calib_plusenu5, 1, 19, 16,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.CalibSystematicStudy(Summarize_plusenu5, Calib_plusenu5, 1, 19, 17,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.CalibSystematicStudy(Summarize_plusenu5, Calib_plusenu5, 1, 19, 20,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.CalibSystematicStudy(Summarize_plusenu5, Calib_plusenu5, 1, 19, 21,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.CalibSystematicStudy(Summarize_plusenu5, Calib_plusenu5, 1, 19, 22,  "Wplusenu5", "$p^{T}_{W}$")

#Optimisation1B.RecoilSystematicStudy(Summarize_plusenu5, Recoil_plusenu5, 1, 19, 12,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.RecoilSystematicStudy(Summarize_plusenu5, Recoil_plusenu5, 1, 19, 13,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.RecoilSystematicStudy(Summarize_plusenu5, Recoil_plusenu5, 1, 19, 14,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.RecoilSystematicStudy(Summarize_plusenu5, Recoil_plusenu5, 1, 19, 15,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.RecoilSystematicStudy(Summarize_plusenu5, Recoil_plusenu5, 1, 19, 16,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.RecoilSystematicStudy(Summarize_plusenu5, Recoil_plusenu5, 1, 19, 17,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.RecoilSystematicStudy(Summarize_plusenu5, Recoil_plusenu5, 1, 19, 20,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.RecoilSystematicStudy(Summarize_plusenu5, Recoil_plusenu5, 1, 19, 21,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.RecoilSystematicStudy(Summarize_plusenu5, Recoil_plusenu5, 1, 19, 22,  "Wplusenu5", "$p^{T}_{W}$")

#Optimisation1B.TotalSystematicStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 12,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalSystematicStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 13,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalSystematicStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 14,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalSystematicStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 15,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalSystematicStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 16,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalSystematicStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 17,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalSystematicStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 20,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalSystematicStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 21,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalSystematicStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 22,  "Wplusenu5", "$p^{T}_{W}$")

#Optimisation1B.TotalStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 12,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 13,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 14,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 15,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 16,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 17,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 20,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 21,  "Wplusenu5", "$p^{T}_{W}$")
#Optimisation1B.TotalStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 22,  "Wplusenu5", "$p^{T}_{W}$")


""" ********************************************************************************************************************************************************* """
""" ***************************************************************** Optimisation Study ******************************************************************** """
""" ********************************************************************************************************************************************************* """

#Optimisation.StatStudy(Summarize_plusenu5, 1, 20, 25, 60,  "Wplusenu5", "$p^{T}_{l}$")    # define the number of iterations for the study
#Optimisation.StatStudy(Summarize_plusenu5, 1, 20, 30, 55,  "Wplusenu5", "$p^{T}_{l}$")    # define the number of iterations for the study
#Optimisation.StatStudy(Summarize_plusenu5, 1, 20, 30, 50,  "Wplusenu5", "$p^{T}_{l}$")    # define the number of iterations for the study
#Optimisation.StatStudy(Summarize_plusenu5, 1, 20, 35, 45,  "Wplusenu5", "$p^{T}_{l}$")    # define the number of iterations for the study

#Optimisation.BiasStudy(Summarize_plusenu5, Bias_plusenu5, 1, 20, 25, 60,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.BiasStudy(Summarize_plusenu5, Bias_plusenu5, 1, 20, 30, 55,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.BiasStudy(Summarize_plusenu5, Bias_plusenu5, 1, 20, 30, 50,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.BiasStudy(Summarize_plusenu5, Bias_plusenu5, 1, 20, 35, 45,  "Wplusenu5", "$p^{T}_{l}$")

#Optimisation.EffSystematicStudy(Summarize_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, 1, 19, 25,  60,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.EffSystematicStudy(Summarize_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, 1, 19, 30,  55,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.EffSystematicStudy(Summarize_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, 1, 19, 30,  50,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.EffSystematicStudy(Summarize_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, 1, 19, 35,  45,  "Wplusenu5", "$p^{T}_{l}$")

#Optimisation.CalibSystematicStudy(Summarize_plusenu5, Calib_plusenu5, 1, 19, 25,   60,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.CalibSystematicStudy(Summarize_plusenu5, Calib_plusenu5, 1, 19, 30,   55,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.CalibSystematicStudy(Summarize_plusenu5, Calib_plusenu5, 1, 19, 30,   50,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.CalibSystematicStudy(Summarize_plusenu5, Calib_plusenu5, 1, 19, 35,   45,  "Wplusenu5", "$p^{T}_{l}$")

#Optimisation.RecoilSystematicStudy(Summarize_plusenu5, Recoil_plusenu5, 1, 19, 25,   60,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.RecoilSystematicStudy(Summarize_plusenu5, Recoil_plusenu5, 1, 19, 30,   55,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.RecoilSystematicStudy(Summarize_plusenu5, Recoil_plusenu5, 1, 19, 30,   50,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.RecoilSystematicStudy(Summarize_plusenu5, Recoil_plusenu5, 1, 19, 35,   45,  "Wplusenu5", "$p^{T}_{l}$")

#Optimisation.TotalSystematicStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 25,   60,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.TotalSystematicStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 30,   55,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.TotalSystematicStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 30,   50,  "Wplusenu5", "$p^{T}_{l}$")
#Optimisation.TotalSystematicStudy(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, 1, 19, 35,   45,  "Wplusenu5", "$p^{T}_{l}$")

""" ********************************************************************************************************************************************************* """
""" ******************************************************************* Nominal 2D Plot ********************************************************************* """
""" ********************************************************************************************************************************************************* """

#MatrixPlots.StatCovarianceMatrix(Summarize_plusenu5, 2,  "Wplusenu5",   "W^{+}#rightarrow e^{+}#nu, 5TeV", "5TeV")
#MatrixPlots.BiasCovarianceMatrix(Bias_plusenu5,      2,  "Wplusenu5",   "W^{+}#rightarrow e^{+}#nu, 5TeV", "5TeV")

""" ********************************************************************************************************************************************************* """
""" ******************************************************************* Differential Xs ********************************************************************* """
""" ********************************************************************************************************************************************************* """
 
#Summarize_plusenu5    = ROOT.TFile.Open("/afs/cern.ch/work/h/hatmani/HistMaker_v1/Plotting/RooUnfold/output_Wplusenu5/elEtaSF/Summarize_Wplusenu5.root")
#Summarize_plusmunu5   = ROOT.TFile.Open("/afs/cern.ch/work/h/hatmani/HistMaker_v1/Plotting/RooUnfold/output_Wplusmunu5/muEtaSF/Summarize_Wplusmunu5.root")
#Summarize_minusenu5   = ROOT.TFile.Open("/afs/cern.ch/work/h/hatmani/HistMaker_v1/Plotting/RooUnfold/output_Wminusenu5/elEtaSF/Summarize_Wminusenu5.root")
#Summarize_minusmunu5  = ROOT.TFile.Open("/afs/cern.ch/work/h/hatmani/HistMaker_v1/Plotting/RooUnfold/output_Wminusmunu5/muEtaSF/Summarize_Wminusmunu5.root")

#DiffCrossSectionX.GetDiffernetialXsComp(Summarize_plusenu5, Summarize_plusmunu5, Summarize_minusenu5, Summarize_minusmunu5, list1, "5TeV", "W$^{+}$ $\\rightarrow$ e$^{+} \\nu $, 5TeV, Uncertainties in (\%)", "Wplusenu5", 256.827)

#CrossSectionDeter.GetDiffernetialXs(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, "5TeV", "W$^{+}$ $\\rightarrow$ e$^{+} \\nu $, 5TeV, Uncertainties in (\%)", "Wplusenu5", 256.827)
#CrossSectionDeter.GetDiffernetialXsPlot(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, "5TeV", "W^{+}#rightarrow e^{+}#nu, 5TeV", "Wplusenu5", 256.827)
#CrossSectionDeter.GetDiffernetialXsPlotN(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, "5TeV", "W^{+}#rightarrow e^{+}#nu, 5TeV", "Wplusenu5", 256.827)


""" ********************************************************************************************************************************************************* """
""" ********************************************************************* fiducial Xs *********************************************************************** """
""" ********************************************************************************************************************************************************* """

#CrossSectionDeter.GetFiducialXs(Summarize_minusenu5, Bias_minusenu5, TrigSF_minusenu5, RecoSF_minusenu5, IsoSF_minusenu5, IdSF_minusenu5, Calib_minusenu5, Recoil_minusenu5, "5TeV", "W^{-} \\rightarrow e^{-} \\nu, 5TeV", "Wminusenu5", 256.827)
#CrossSectionDeter.GetSummaringTable(Summarize_plusenu5, Bias_plusenu5, TrigSF_plusenu5, RecoSF_plusenu5, IsoSF_plusenu5, IdSF_plusenu5, Calib_plusenu5, Recoil_plusenu5, "5TeV", "W$^{+}$ $\\rightarrow$ e$^{+} \\nu $, 5TeV, Uncertainties in (\%)", "Wplusenu5", 256.827, 0.460)

""" ********************************************************************************************************************************************************* """
""" ******************************************************************** Background Plot ******************************************************************** """
""" ********************************************************************************************************************************************************* """
'''
channel             = "Wminusenu"
Energy              = "5TeV"
Indice              = "W^{-}#rightarrow e^{-}#nu, 5TeV"
data                = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/pT_Input/data17_minusenu_5TeV.root")
Signal              = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/pT_Input/mc16_5TeV.361103.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminusenu.e4916_s3238_r10243_r10210_p3665.root")
Background_Top      = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wminusenu5/Background_Top.root")
Background_diboson  = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wminusenu5/Background_dilepton.root")
Background_W        = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wminusenu5/Background_W.root")
Background_Z        = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wminusenu5/Background_Z.root")
Background_MiltiJet = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wminusenu5/Background_MiltiJet.root")
BackgroundPlot.BackgroundPlotspTlepton(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)
BackgroundPlot.BackgroundPlotsetalepton(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)
BackgroundPlot.BackgroundPlotsmTW(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)
BackgroundPlot.BackgroundPlotspTW(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)

channel             = "Wplusenu"
Energy              = "5TeV"
Indice              = "W^{+}#rightarrow e^{+}#nu, 5TeV"
data                = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/pT_Input/data17_plusenu_5TeV.root")
Signal              = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/pT_Input/mc16_5TeV.361100.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu.e4916_s3238_r10243_r10210_p3665.root")
Background_Top      = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wplusenu5/Background_Top.root")
Background_diboson  = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wplusenu5/Background_dilepton.root")
Background_W        = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wplusenu5/Background_W.root")
Background_Z        = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wplusenu5/Background_Z.root")
Background_MiltiJet = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wplusenu5/Background_MiltiJet.root")
BackgroundPlot.BackgroundPlotspTlepton(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)
BackgroundPlot.BackgroundPlotsetalepton(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)
BackgroundPlot.BackgroundPlotsmTW(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)
BackgroundPlot.BackgroundPlotspTW(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)

channel             = "Wminusmunu"
Energy              = "5TeV"
Indice              = "W^{-}#rightarrow #mu^{-}#nu, 5TeV"
data                = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/pT_Input/data17_minusmunu_5TeV.root")
Signal              = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/pT_Input/mc16_5TeV.361104.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminusmunu.e4916_s3238_r10243_r10210_p3665.root")
Background_Top      = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wminusmunu5/Background_Top.root")
Background_diboson  = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wminusmunu5/Background_dilepton.root")
Background_W        = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wminusmunu5/Background_W.root")
Background_Z        = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wminusmunu5/Background_Z.root")
Background_MiltiJet = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wminusmunu5/Background_MiltiJet.root")
BackgroundPlot.BackgroundPlotspTlepton(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)
BackgroundPlot.BackgroundPlotsetalepton(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)
BackgroundPlot.BackgroundPlotsmTW(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)
BackgroundPlot.BackgroundPlotspTW(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)

channel             = "Wplusmunu"
Energy              = "5TeV"
Indice              = "W^{+}#rightarrow #mu^{+}#nu, 5TeV"
data                = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/pT_Input/data17_plusmunu_5TeV.root")
Signal              = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/pT_Input/mc16_5TeV.361101.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu.e4916_s3238_r10243_r10210_p3665.root")
Background_Top      = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wplusmunu5/Background_Top.root")
Background_diboson  = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wplusmunu5/Background_dilepton.root")
Background_W        = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wplusmunu5/Background_W.root")
Background_Z        = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wplusmunu5/Background_Z.root")
Background_MiltiJet = ROOT.TFile.Open("/Users/hatmani/pTw_Unfolding/Inputs/Background/Wplusmunu5/Background_MiltiJet.root")
BackgroundPlot.BackgroundPlotspTlepton(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)
BackgroundPlot.BackgroundPlotsetalepton(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)
BackgroundPlot.BackgroundPlotsmTW(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)
BackgroundPlot.BackgroundPlotspTW(data,  Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy)
'''
