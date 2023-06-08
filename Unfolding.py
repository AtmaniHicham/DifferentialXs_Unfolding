#!/usr/bin/env python
# -*-coding:Latin-1 -*

# ==============================================================================
#  Description:
#       Unfolding using for pT(W), mT(W), pT(lepton) and eta(lepton)
#  Author:  ATMANI Hicham, based on RooUnfold framework
# ==============================================================================
import sys
import os
import time
import gc
import csv
import numpy as np
import ROOT
from ROOT import TFile, TH1D, TH1F, TH2F, TCanvas, TPad, TLegend, gStyle, gROOT, gPad, gDirectory, TVector2, TPaveStats, TStyle, TLatex
from ROOT import TColor, kBlack, kRed, kBlue, kMagenta, kYellow, kCyan, kGreen, kOrange, kTeal, kPink, kGray
from ROOT import TArrayD, TAxis, TMath, TVectorF, TMatrixF
from ROOT import kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal

from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
#from root_numpy import *
from source.source import *
import yaml


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def CrossSectionNormalisa(reco_hist, channel, NormalisationPath):

    # initialisation
    lepton = "ee"
    if "munu" in channel: lepton = "mumu"
    NameVector = ["ttbar", "singletop", "Wtop", "Wminusenu", "Wplusenu", "Wminusmunu", "Wplusmunu", "Wplustaunu", "Wminustaunu", "Z"+lepton, "Ztautau", "ZZ", "WZ", "WW", "multijet"]

    # remove all files starting with "._"
    for filename in os.listdir(NormalisationPath):
        if filename.startswith("."):
            os.remove(os.path.join(NormalisationPath, filename))

    # read the root files and make one list of the files:
    VariationFiles   = []
    files            = os.listdir(NormalisationPath)
    for i in files:
        input = ROOT.TFile.Open(NormalisationPath+"/"+i)
        VariationFiles.append(input)

    # read the root files and make list of histogram:
    VariationHisto   = []
    DireName         = channel + "Selection_NormVariations"
    for i in VariationFiles:
        for j in NameVector:
            HistVaried   = i.Get(DireName+"/Lepton_pT_Reco_cut7_Xsec_up_"+j)
            if (HistVaried != None):
                VariationHisto.append(HistVaried)


    #print(VariationHisto)
    #print(" hello ", VariationHisto[0])
    return VariationHisto
    
# ============================================================================================================================================================
#  				                                           Read the config file
# ============================================================================================================================================================

print(sys.argv[1])

# ====================================================================== read config file ====================================================================================

print("************ read config file ************")
print("Energy           : ", sys.argv[1])
print("Channel          : ", sys.argv[2])
print("Charge 	        : ", sys.argv[3])
print("Variable         : ", sys.argv[4])
print("Niterations      : ", sys.argv[5])
print("CentreOfMass     : ", sys.argv[6])
print("Luminosity       : ", sys.argv[7])

if "plusenu"   in sys.argv[2]: MJChannel = "Wep"
if "plusmunu"  in sys.argv[2]: MJChannel = "Wmup"
if "minusenu"  in sys.argv[2]: MJChannel = "Wem"
if "minusmunu" in sys.argv[2]: MJChannel = "Wmum"

print("************ end read config file ************")

# ====================================================================== read input file =====================================================================================

print("************ read input files ************")
Data                = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[3]+"_DATA_"+sys.argv[1]+"/Nominal/DATA.root")
Signal              = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[3]+"_MC_"  +sys.argv[1]+"/Nominal/MC.root")
Background_Top      = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[3]+"_MC_"  +sys.argv[1]+"/Nominal/merge/Nominal_Top.root")
Background_diboson  = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[3]+"_MC_"  +sys.argv[1]+"/Nominal/merge/Nominal_dilepton.root")
Background_W        = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[3]+"_MC_"  +sys.argv[1]+"/Nominal/merge/Nominal_W.root")
Background_Z        = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[3]+"_MC_"  +sys.argv[1]+"/Nominal/merge/Nominal_Z.root")
Background_MiltiJet = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/MJbkg/"+sys.argv[1]+"/"+MJChannel+".root")

if "munu" in sys.argv[2]:
    Signal_Iso    = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w" + sys.argv[3] + "_MC_" + sys.argv[1]+"/MuonSF_iso/MC.root")
    Signal_Trig   = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w" + sys.argv[3] + "_MC_" + sys.argv[1]+"/MuonSF_trig/MC.root")
    Signal_recoSF = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w" + sys.argv[3] + "_MC_" + sys.argv[1]+"/MuonSF_reco/MC.root")
    Signal_ttva   = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w" + sys.argv[3] + "_MC_" + sys.argv[1]+"/MuonSF_ttva/MC.root")

if "enu" in sys.argv[2]:
    Signal_Iso    = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w" + sys.argv[3] + "_MC_" + sys.argv[1]+"/ElecSF_iso/MC.root")
    Signal_Trig   = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w" + sys.argv[3] + "_MC_" + sys.argv[1]+"/ElecSF_trig/MC.root")
    Signal_recoSF = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w" + sys.argv[3] + "_MC_" + sys.argv[1]+"/ElecSF_reco/MC.root")
    Signal_Id     = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w" + sys.argv[3] + "_MC_" + sys.argv[1]+"/ElecSF_id/MC.root")

print(" Data                 : ",	Data			)
print(" Signal               : ",	Signal			)
print(" Background_Top       : ",	Background_Top		)
print(" Background_diboson   : ",	Background_diboson	)
print(" Background_W         : ",	Background_W		)
print(" Background_Z         : ",	Background_Z		)
print(" Background_MiltiJet  : ",	Background_MiltiJet	)

print("************ end read input files ************")

# ====================================================================== calculate the Total Background ======================================================================

Background_Total = SumBackground( Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, sys.argv[2],  sys.argv[4]+"_Reco")

# ====================================================================== read the truth histograms ===========================================================================

truth_hist       = TruthDistribution( Signal,  sys.argv[4], sys.argv[6])

# ====================================================================== read the migration Matrix ===========================================================================

mig_hist         = MigrationMatrix( Signal, sys.argv[2], sys.argv[6], sys.argv[4], truth_hist)

# ====================================================================== read the reco histograms of MC and data =============================================================

data_hist  = DataDistribution( Data, 	sys.argv[2], sys.argv[4]) 
reco_hist  = RecoDistribution( Signal,  sys.argv[2], sys.argv[4])

# ====================================================================== Control table =======================================================================================

print("Background       bins: ",Background_Total.GetNbinsX(),           	"       Integral: ",Background_Total.Integral())
print("truth_hist       bins: ",truth_hist.GetNbinsX(),                  "       Integral: ",truth_hist.Integral())
print("data_hist        bins: ",data_hist.GetNbinsX(),                  	"       Integral: ",data_hist.Integral())
print("reco_hist        bins: ",reco_hist.GetNbinsX(),                  	"       Integral: ",reco_hist.Integral())
print("B",              Background_Total  )
print("M",              mig_hist          )
print("T",              truth_hist        )
print("D",              data_hist         )
print("R",              reco_hist         )

# ====================================================================== get the Acceptance and Efficiency ===================================================================

Acceptance_hist = GetAcceptance( mig_hist, truth_hist )
Efficiency_hist = GetEfficieny(  mig_hist, reco_hist  )

# ====================================================================== apply the efficieny and subrtact the background from data ===========================================

dataCorrected = CorrectData( data_hist, reco_hist, Background_Total, Efficiency_hist)

# ====================================================================== define the Response Matrix ==========================================================================

responseM       = RooUnfoldResponse(0,0,mig_hist,"UNFOLD","UNFOLD");
unfolded_data   = []
Covarinace_data = []
for i in range(1, 1+int(sys.argv[5])):
    unfoldDATA = RooUnfoldBayes (responseM, dataCorrected, i)
    unfolded_data.append( unfoldDATA.Hreco())
    Covarinace_data.append( unfoldDATA.Ereco(2) )

# ====================================================================== closure test for MC =================================================================================

unfolded_MC   = []
Covarinace_MC = []
for i in range(1, 1+int(sys.argv[5])):
    unfoldMC = RooUnfoldBayes (responseM, mig_hist.ProjectionX(), i)
    unfolded_MC.append( unfoldMC.Hreco())
    Covarinace_MC.append( unfoldMC.Ereco(2) )

# ====================================================================== Background uncertainties =================================================================================
"""
# -------- Normalisation path
NormalisationPath   = "/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[3]+"_MC_"+sys.argv[1]+"/Normalization/"

# -------- initialisation
lepton = "ee"
if "munu" in sys.argv[2]: lepton = "mumu"
NameVector = ["ttbar", "singletop", "Wtop", "Wminusenu", "Wplusenu", "Wminusmunu", "Wplusmunu", "Wplustaunu", "Wminustaunu", "Z"+lepton, "Ztautau", "ZZ", "WZ", "WW", "multijet"]

# -------- remove all files starting with "._"
for filename in os.listdir(NormalisationPath):
    if filename.startswith("."):
        os.remove(os.path.join(NormalisationPath, filename))

# -------- read the root files and make one list of the files:
VariationFiles   = []
files            = os.listdir(NormalisationPath)
for i in files:
    input = ROOT.TFile.Open(NormalisationPath+"/"+i)
    VariationFiles.append(input)

# -------- read the root files and make list of histogram:
VariationHisto   = []
DireName         = sys.argv[2] + "Selection_NormVariations"
for i in VariationFiles:
    for j in NameVector:
        HistVaried   = i.Get(DireName+"/Lepton_pT_Reco_cut7_Xsec_up_"+j)
        if (HistVaried != None):
            VariationHisto.append(HistVaried)

print(VariationHisto[0].GetBinContent(3))
# -------- redo the unfolding for each background one:

NormalisationUnfolded       = []
for i in VariationHisto:
    Background_TotalVaried  = Background_Total.Clone()
    Background_TotalVaried.Add(i)
    CorrectedData           = CorrectData( data_hist, reco_hist, Background_TotalVaried, Efficiency_hist)
    unfoldDataNormalisation = RooUnfoldBayes (responseM, dataCorrected, 4)
    NormalisationUnfolded.append( unfoldDataNormalisation.Hreco())

OutputFile = TFile("Output/output_"+sys.argv[2]+str(sys.argv[6])+"/"+sys.argv[4]+"/Normalisation_"+sys.argv[2]+str(sys.argv[6])+".root",'RECREATE')
k=0
for i in NormalisationUnfolded:
    i.Write("Normalisation"+str(k))
    k=k+1

# ======================================================================  define the bias for Unfolding ======================================================================
"""
if(sys.argv[4] == "Lepton_Eta"):
    VarMin = -2.4
    VarMax =  2.4
if(sys.argv[4] == "Lepton_pT"):
    VarMin =  0
    VarMax =  100
"""
GetTheBias(responseM, dataCorrected, mig_hist, int(sys.argv[5]), VarMin, VarMax, sys.argv[2], sys.argv[6], sys.argv[4])

# ====================================================================== define the background uncertainties =================================================================

MultiJetUncertai = MultiJetUncerta(             Background_Top, Background_MiltiJet, sys.argv[2],  sys.argv[4]+"_Reco")
LuminosityUncert = LuminosityUncer( reco_hist,  Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, sys.argv[2],  sys.argv[4]+"_Reco", sys.argv[1])
BackgrodUncertai = BackgroundUncer(             Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, sys.argv[2],  sys.argv[4]+"_Reco", sys.argv[1])

# ====================================================================== Calculate the Sf Systematics ========================================================================
"""

if "enu" in sys.argv[2]:
	Id_Systematics   = GetSFSystematic( Signal_Id,    reco_hist, mig_hist, sys.argv[2], sys.argv[1], sys.argv[4],   "ElIDSys",      sys.argv[5], VarMin, VarMax,  sys.argv[6]) 


if "munu" in sys.argv[2]:
	Id_Systematics   = GetSFSystematic( Signal_ttva,   reco_hist, mig_hist, sys.argv[2], sys.argv[1], sys.argv[4],  "MuTTVASys",    sys.argv[5], VarMin, VarMax,  sys.argv[6]) 

Iso_Systematics  = GetSFSystematic(     Signal_Iso,   reco_hist, mig_hist, sys.argv[2], sys.argv[1], sys.argv[4],   "IsoSys",       sys.argv[5], VarMin, VarMax,  sys.argv[6])
Reco_Systematics = GetSFSystematic(     Signal_recoSF,  reco_hist, mig_hist, sys.argv[2], sys.argv[1], sys.argv[4], "RecoSys",      sys.argv[5], VarMin, VarMax,  sys.argv[6])
Trig_Systematics = GetSFSystematic(     Signal_Trig,  reco_hist, mig_hist, sys.argv[2], sys.argv[1], sys.argv[4],   "TrigSys",      sys.argv[5], VarMin, VarMax,  sys.argv[6])

# ====================================================================== Calculate the Calib Systematics =====================================================================

if  "enu" in sys.argv[2]:
	CalibCovarianceItera  = []
	CalibVariation        = GetTheCalibVariation(sys.argv[3], sys.argv[1])
	GetPrinSystematic( CalibVariation, reco_hist, mig_hist, sys.argv[2],  sys.argv[1], sys.argv[4], sys.argv[5], VarMin, VarMax, 1, "Calib", sys.argv[6])

if  "munu" in sys.argv[2]:
	CalibCovarianceItera  = []
	CalibVariation        = GetTheCalibVariationMu(sys.argv[3], sys.argv[1])
	GetPrinSystematic( CalibVariation, reco_hist, mig_hist, sys.argv[2],  sys.argv[1], sys.argv[4], sys.argv[5], VarMin, VarMax, 1, "Calib", sys.argv[6])

# ====================================================================== Calculate the Recoil Systematics ====================================================================

RecoilCovarianceItera  = []
RecoilVariation        = GetTheRecoilVariation(sys.argv[3], sys.argv[1])
GetPrinSystematic( RecoilVariation, reco_hist, mig_hist, sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[5], VarMin, VarMax, 1, "Recoil", sys.argv[6])

# ====================================================================== Save Output File ====================================================================================

OutputFile = TFile("Output/output_"+sys.argv[2]+str(sys.argv[6])+"/"+sys.argv[4]+"/Summarize_"+sys.argv[2]+str(sys.argv[6])+".root",'RECREATE')

mig_hist.Write("mig_hist")
truth_hist.Write("truth_hist")
reco_hist.Write("reco_hist")

Background_Total.Write("Background_Total")

(mig_hist.ProjectionX()).Write("reco_prjection")
(mig_hist.ProjectionY()).Write("truth_prjection")

Efficiency_hist.Write("Efficiency_hist")
Acceptance_hist.Write("Acceptance_hist")

data_hist.Write("data_hist")
dataCorrected.Write("dataCorrected")

for i in range(0, 0+int(sys.argv[5])):
    unfolded_data[i].Write("unfolded_data"+str(i+1))
    unfolded_MC[i].Write("unfolded_MC"+str(i+1))
for i in range(0, len(Covarinace_data) ):
    Covarinace_data[i].Write("CovarianceMatrix_Data_Iter"+str(i+1))       
    Covarinace_MC[i].Write("CovarianceMatrix_MC_Iter"+str(i+1))

# ------------------------------------------- Input for the combination -------------------------------------------------------
# The unfolded distribution:
unfolded_data_subtr_pur_eff = unfolded_data[1].Clone("UnfoldedData")
for i in range(0, (unfolded_data_subtr_pur_eff.GetNbinsX())):
    if Acceptance_hist.GetBinContent(i+1) != 0:
        unfolded_data_subtr_pur_eff.SetBinContent(i+1, unfolded_data_subtr_pur_eff.GetBinContent(i+1)/Acceptance_hist.GetBinContent(i+1))
    if Acceptance_hist.GetBinContent(i+1) == 0:
        unfolded_data_subtr_pur_eff.SetBinContent(i+1, 0)

# statistical uncertainty of data and MC
StatData = unfolded_data[1].Clone("StatData")
StatMC   = unfolded_MC[1].Clone("StatMC")

for i in range(0, (StatData.GetNbinsX())):
    StatData.SetBinContent(i+1, sqrt(Covarinace_data[1](i,i))   )
    StatMC.SetBinContent(i+1,   sqrt(Covarinace_MC[1](i,i))     )

# correlation matrix for data and simulation
corrMatrix_statdata_pur_eff = mig_hist.Clone("corrMatrix_statdata_pur_eff")
corrMatrix_statMC_pur_eff = mig_hist.Clone("corrMatrix_statMC_pur_eff")


for i in range(0, corrMatrix_statdata_pur_eff.GetNbinsX()):
    for j in range(0, corrMatrix_statdata_pur_eff.GetNbinsY()):
        Multi = Covarinace_data[1](i, i)*Covarinace_data[1](j, j)
        if Multi != 0:
            corrMatrix_statdata_pur_eff.SetBinContent(i+1, j+1,  Covarinace_data[1](i,j) / (sqrt(Covarinace_data[1](i,i))*sqrt(Covarinace_data[1](j,j))))
        if Multi == 0:
            corrMatrix_statdata_pur_eff.SetBinContent(i+1, j+1,   0)

for i in range(0, corrMatrix_statdata_pur_eff.GetNbinsX()):
    for j in range(0, corrMatrix_statdata_pur_eff.GetNbinsY()):
        Multi = Covarinace_MC[1](i, i)*Covarinace_MC[1](j, j)
        if Multi != 0:
            corrMatrix_statMC_pur_eff.SetBinContent(i+1, j+1,    Covarinace_MC[1](i, j)    / (sqrt(Covarinace_MC[1](i, i))   * sqrt(Covarinace_MC[1](j, j))))
            print(i+1, j+1,    Covarinace_MC[1](i, j)    / (sqrt(Covarinace_MC[1](i, i))   * sqrt(Covarinace_MC[1](j, j))))
        if Multi == 0:
            corrMatrix_statMC_pur_eff.SetBinContent(i+1, j+1,   0)

unfolded_data_subtr_pur_eff.Write("unfolded_data_subtr_pur_eff")
StatData.Write("StatData")
StatMC.Write("StatMC")
corrMatrix_statdata_pur_eff.Write("corrMatrix_statdata_pur_eff")
corrMatrix_statMC_pur_eff.Write("corrMatrix_statMC_pur_eff")


print("end")
os._exit(0)
