#!/usr/bin/env python
# -*-coding:Latin-1 -*

import ROOT
import ROOT as root
from math import *
import os 
import numpy as np
import csv

from ROOT import TFile, TH1F, TH2F, TCanvas, TPad, TLegend, gStyle, gROOT, gPad, gDirectory, TVector2, TPaveStats, TStyle, TLatex
from ROOT import TColor, kBlack, kRed, kBlue, kMagenta, kYellow, kCyan, kGreen, kOrange, kTeal, kPink, kGray
from ROOT import TArrayD, TAxis, TMath, TVectorF, TMatrixF, TF1, TH2D, TH1D
from ROOT import kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal

from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
from root_numpy import *



# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def rebinTheMultiJetBackground(nominal, histogram_MJ):

    nominal.SetBinContent(1, 0)
    nominal.SetBinContent(2, histogram_MJ.GetBinContent(26) + histogram_MJ.GetBinContent(27) + histogram_MJ.GetBinContent(28) + histogram_MJ.GetBinContent(29) + histogram_MJ.GetBinContent(30))
    nominal.SetBinContent(3, histogram_MJ.GetBinContent(31) + histogram_MJ.GetBinContent(32) + histogram_MJ.GetBinContent(33) + histogram_MJ.GetBinContent(34) + histogram_MJ.GetBinContent(35))
    nominal.SetBinContent(4, histogram_MJ.GetBinContent(36) + histogram_MJ.GetBinContent(37) + histogram_MJ.GetBinContent(38) + histogram_MJ.GetBinContent(39) + histogram_MJ.GetBinContent(40))
    nominal.SetBinContent(5, histogram_MJ.GetBinContent(41) + histogram_MJ.GetBinContent(42) + histogram_MJ.GetBinContent(43) + histogram_MJ.GetBinContent(44) + histogram_MJ.GetBinContent(45))
    nominal.SetBinContent(6, histogram_MJ.GetBinContent(46) + histogram_MJ.GetBinContent(47) + histogram_MJ.GetBinContent(48) + histogram_MJ.GetBinContent(49) + histogram_MJ.GetBinContent(50))

    sum1 = 0
    for i in range(51, 61):
        sum1 = sum1 + histogram_MJ.GetBinContent(i)

    sum2 = 0
    for i in range(61, 81):
        sum2 = sum2 + histogram_MJ.GetBinContent(i)

    sum3 = 0
    for i in range(81, 101):
        sum3 = sum3 + histogram_MJ.GetBinContent(i)

    nominal.SetBinContent(7, sum1)
    nominal.SetBinContent(8, sum2)
    nominal.SetBinContent(9, sum3)

    nominal.SetBinError(1, 0)
    nominal.SetBinError(2, sqrt(pow(histogram_MJ.GetBinError(26), 2) + pow(histogram_MJ.GetBinError(27), 2) + pow(histogram_MJ.GetBinError(28), 2) + pow(histogram_MJ.GetBinError(29), 2) + pow(histogram_MJ.GetBinError(30), 2)))
    nominal.SetBinError(3, sqrt(pow(histogram_MJ.GetBinError(31), 2) + pow(histogram_MJ.GetBinError(32), 2) + pow(histogram_MJ.GetBinError(33), 2) + pow(histogram_MJ.GetBinError(34), 2) + pow(histogram_MJ.GetBinError(35), 2)))
    nominal.SetBinError(4, sqrt(pow(histogram_MJ.GetBinError(36), 2) + pow(histogram_MJ.GetBinError(37), 2) + pow(histogram_MJ.GetBinError(38), 2) + pow(histogram_MJ.GetBinError(39), 2) + pow(histogram_MJ.GetBinError(40), 2)))
    nominal.SetBinError(5, sqrt(pow(histogram_MJ.GetBinError(41), 2) + pow(histogram_MJ.GetBinError(42), 2) + pow(histogram_MJ.GetBinError(43), 2) + pow(histogram_MJ.GetBinError(44), 2) + pow(histogram_MJ.GetBinError(45), 2)))
    nominal.SetBinError(6, sqrt(pow(histogram_MJ.GetBinError(46), 2) + pow(histogram_MJ.GetBinError(47), 2) + pow(histogram_MJ.GetBinError(48), 2) + pow(histogram_MJ.GetBinError(49), 2) + pow(histogram_MJ.GetBinError(50), 2)))

    errsum1 = 0
    for i in range(51, 61):
        errsum1 = errsum1 + pow(histogram_MJ.GetBinError(i), 2)

    errsum2 = 0
    for i in range(61, 81):
        errsum2 = errsum2 + pow(histogram_MJ.GetBinError(i), 2)

    errsum3 = 0
    for i in range(81, 101):
        errsum3 = errsum3 + pow(histogram_MJ.GetBinError(i), 2)

    nominal.SetBinError(7, sqrt(errsum1))
    nominal.SetBinError(8, sqrt(errsum2))
    nominal.SetBinError(9, sqrt(errsum3))

    return nominal

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def GetTheRecoilVariation(channelName, Energy):

    RecoilSystVariation = []
    if(Energy == "5TeV"):
        #RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varSET_SYSbin1.root'))
        RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varRESPONSE_EXTSYS_DOWNbin1.root'))
        RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varRESOLUTION_EXTSYS_DOWNbin1.root'))
        RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varRESPONSE_SYS_DOWNbin1.root'))

        for i in range(3, 20):
            RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varRESPONSE_STAT0_DOWNbin'+str(i)+'.root'))
            RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varRESPONSE_STAT1_DOWNbin'+str(i)+'.root'))  # error

        for i in range(1, 13):
            RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varRESOLUTION_STAT0_DOWNbin'+str(i)+'.root'))
            RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varRESOLUTION_STAT1_DOWNbin'+str(i)+'.root'))

    if(Energy == "13TeV"):
        #RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varSET_SYSbin1.root'))
        RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varRESPONSE_EXTSYS_DOWNbin1.root'))
        RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varRESOLUTION_EXTSYS_DOWNbin1.root'))
        RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varRESPONSE_SYS_DOWNbin1.root'))

        for i in range(1, 16):
            RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varRESPONSE_STAT0_DOWNbin'+str(i)+'.root'))
            RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varRESPONSE_STAT1_DOWNbin'+str(i)+'.root'))

        for i in range(1, 15):
            RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varRESOLUTION_STAT0_DOWNbin'+str(i)+'.root'))
            RecoilSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+Energy+'/RecoilVar/merge/mc16_'+Energy+'.varRESOLUTION_STAT1_DOWNbin'+str(i)+'.root'))

    print(len(RecoilSystVariation))

    for i in range(0, len(RecoilSystVariation)):
        print(i, RecoilSystVariation[i])

    return RecoilSystVariation

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def GetTheCalibVariationMu(channelName, Energy):

    CalibSystVariation = []

    CalibSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channelName + '_MC_'+str(Energy)+'/MuCalibVar/merge/MC_varMUON_ID_1down.root'))
    CalibSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channelName + '_MC_'+str(Energy)+'/MuCalibVar/merge/MC_varMUON_MS_1down.root'))
    CalibSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channelName + '_MC_'+str(Energy)+'/MuCalibVar/merge/MC_varMUON_SCALE_1down.root'))
    CalibSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channelName + '_MC_'+str(Energy)+'/MuCalibVar/merge/MC_varSagittaBiasOffsetDown.root'))

    for i in range(1, 25):
        CalibSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channelName + '_MC_'+str(Energy)+'/MuCalibVar/merge/mc16_'+str(Energy)+'.varSagittaBiasstatDownbin'+str(i)+'.root'))

    return CalibSystVariation

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def GetTheCalibVariation(channelName, Energy):

    CalibSystVariation = []

    for i in range(1, 25):
        CalibSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+str(Energy)+'/ElCalibVar/merge/mc16_'+str(Energy)+'.varscaleDownbin'+str(i)+'.root'))
        print('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channelName +'_MC_'+str(Energy)+'/ElCalibVar/merge/mc16_'+str(Energy)+'.varscaleDownbin'+str(i)+'.root')
    for i in range(1, 25):
        CalibSystVariation.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' +channelName+'_MC_'+str(Energy)+'/ElCalibVar/merge/mc16_'+str(Energy)+'.varcDownbin'+str(i)+'.root'))

    return CalibSystVariation

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def GetDevSystematic(reco_hist, reco_Varied, mig_hist, responseM, itera):

    reco_hist_Eff = GetEfficieny(mig_hist, reco_hist)
    reco_Varied_Eff = GetEfficieny(mig_hist, reco_Varied)

    responseM = RooUnfoldResponse(0, 0, mig_hist, "UNFOLD", "UNFOLD")

    unfoldMC = RooUnfoldBayes(responseM, reco_hist_Eff,   itera)
    HistNominal = unfoldMC.Hreco()

    unfoldMC_Var = RooUnfoldBayes(responseM, reco_Varied_Eff, itera)
    HistVaried = unfoldMC_Var.Hreco()

    CovarianceMatrix = SysCovarianc(HistNominal, HistVaried)

    return CovarianceMatrix

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def SysCovarianc(HistoNominal, unfoldedSys_down):

    ntbins = HistoNominal.GetNbinsX()
    xaxis = HistoNominal.GetXaxis()

    Covariance = TH2D("CovarianceMatrix", "CovarianceMatrix", ntbins,
                      xaxis.GetXbins().GetArray(), ntbins, xaxis.GetXbins().GetArray())

    for j in range(1, 1 + HistoNominal.GetNbinsX()):
        for k in range(1, 1 + HistoNominal.GetNbinsX()):
            Covariance.SetBinContent(j, k, (unfoldedSys_down.GetBinContent(
                j)-HistoNominal.GetBinContent(j))*(unfoldedSys_down.GetBinContent(k)-HistoNominal.GetBinContent(k)))

    return Covariance

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def SumCovarianceMatrix(CalibCovarianceMatrix):

    CovarianceTotal = (CalibCovarianceMatrix[0]).Clone("CovarianceTotal")

    print("nomber of event", len(CalibCovarianceMatrix))
    for j in range(1, 1 + CovarianceTotal.GetNbinsX()):
        for k in range(1, 1 + CovarianceTotal.GetNbinsX()):
            CovarianceTotal.SetBinContent(j, k, 0)

    for i in range(1, 1 + CovarianceTotal.GetNbinsX()):
        for j in range(1, 1 + CovarianceTotal.GetNbinsX()):
            covsum = 0
            for k in range(0, len(CalibCovarianceMatrix)):
                covsum = covsum + CalibCovarianceMatrix[k].GetBinContent(i, j)
            CovarianceTotal.SetBinContent(i, j, covsum)

    return CovarianceTotal

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def GetPrinSystematic(CalibVariationFiles, reco_hist, mig_hist, Channel, Energy, Variable, Niter, VarMin, VarMax, Muon, Syst, CofM):

    print(len(CalibVariationFiles))

    # read the variations, and the migration matrix correspond to each variation:
    CalibVariation           = []
    CalibVariation_Eff       = []
    CalibVariation_Migration = []

    for i in range(0, len(CalibVariationFiles)):
        CalibVariation.append(CalibVariationFiles[i].Get(Channel + "Selection/"+Variable+"_Reco_cut7"))
        CalibVariation_Migration.append(CalibVariationFiles[i].Get(Channel + "Selection/"+Variable+"_Reco_v_Truth_cut7"))

    # applied the efficiency correction factor to all the distributions:
    for i in range(0, len(CalibVariation)):
        CalibVariation_Eff.append(ApplyEffecciency(CalibVariation[i], reco_hist, mig_hist))

    # correct the reco distribution:
    reco_hist_Eff = (mig_hist.ProjectionX()).Clone("reco_hist_Eff")

    # define the nominal unfolding response matrix for Syst unfo
    response_Nominal = RooUnfoldResponse(0, 0, mig_hist, "UNFOLD", "UNFOLD")
    unfoldMCNominal  = RooUnfoldBayes(response_Nominal, reco_hist_Eff, 2)
    UnfoldingNominal = unfoldMCNominal.Hreco()

    # define the varied unfolding, with varied migration matrix:
    Systematic_list      = []
    Systematic_Cov_list  = []
    Systmatic_Uncert     = UnfoldingNominal.Clone("Systmatic_Uncert")
    Systmatic_Covariance = mig_hist.Clone("Systmatic_Covariance")

    for i in range(0, len(CalibVariation_Eff)):
        response_varied  = RooUnfoldResponse(0, 0, CalibVariation_Migration[i], "UNFOLD", "UNFOLD")
        unfoldMCVaried   = RooUnfoldBayes(response_varied, CalibVariation_Eff[i], 2)
        UnfoldingVaried  = unfoldMCVaried.Hreco()

        for k in range(0, UnfoldingVaried.GetNbinsX()):
            Systmatic_Uncert.SetBinContent(k+1, UnfoldingVaried.GetBinContent(k+1) - UnfoldingNominal.GetBinContent(k+1))
            Systmatic_Uncert.SetBinError(k+1, 0)
        for k in range(0, Systmatic_Covariance.GetNbinsX()):
            for l in range(0, Systmatic_Covariance.GetNbinsY()):
                Systmatic_Covariance.SetBinContent(k+1, l+1, Systmatic_Uncert.GetBinContent(k+1)*Systmatic_Uncert.GetBinContent(l+1))
        Systematic_list.append(Systmatic_Uncert)
        Systematic_Cov_list.append(Systmatic_Covariance)
    
    OutputFile = TFile("Output/output_"+Channel+str(CofM) +"/"+Variable+"/"+Syst+"_Syst.root", 'RECREATE')
    for i in range(0, len(Systematic_list)):
        Systematic_list[i].Write( "unfoldedSys_Syst" + Syst+str(i))
        Systematic_Cov_list[i].Write( "unfoldedSys_Covariance_" + Syst+str(i))

    """
    CalibVariation = []
    CalibVariation_Eff = []

    for i in range(0, len(CalibVariationFiles)):
        CalibVariation.append(CalibVariationFiles[i].Get(Channel + "Selection/"+Variable+"_Reco_cut7"))
        print(Channel + "Selection/"+Variable+"_Reco_cut7")

    # Add variation to Nominal and apply efficiency corrections:
    for i in range(0, len(CalibVariation)):
        CalibVariation_Eff.append(ApplyEffecciency(CalibVariation[i], reco_hist, mig_hist))

    # correct the reco distribution:
    reco_hist_Eff = (mig_hist.ProjectionX()).Clone("reco_hist_Eff")

    # define the response matrix for Syst unfo
    response_Syst = RooUnfoldResponse(0, 0, mig_hist, "UNFOLD", "UNFOLD")

    # Calculate the Syst:
    Systematics = []
    CovarianceIter = []

    for i in range(1, 1 + int(Niter)):
        unfoldMCNominal  = RooUnfoldBayes(response_Syst, reco_hist_Eff, i)
        HistoNominal     = unfoldMCNominal.Hreco()

        unfoldedSys_down = GetUnfoldToys(response_Syst, CalibVariation_Eff, i)
        Systematic_down  = GetSystematics(HistoNominal, unfoldedSys_down)
        Covariance       = GetSystCovarianceMatrix(HistoNominal, unfoldedSys_down, mig_hist, Variable)

        Systematic_down.GetXaxis().SetRangeUser(VarMin, VarMax)
        Covariance.GetXaxis().SetRangeUser(VarMin, VarMax)
        Covariance.GetYaxis().SetRangeUser(VarMin, VarMax)
        CovarianceIter.append(Covariance)
        Systematics.append(Systematic_down)

    OutputFile = TFile("Output/output_"+Channel+str(CofM) +"/"+Variable+"/"+Syst+"_Syst.root", 'RECREATE')

    for i in range(1, int(Niter)):
        Systematics[i].Write(Syst + "_Systematics_Iter"+str(i))
        CovarianceIter[i].Write(Syst + "_Covariance_Iter"+str(i))

    return Systematics
    """
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def GetSFSystematic(InputSyst, reco_hist, mig_hist, Channel, Energy, Variable, Syst, Niter, VarMin, VarMax, CofM):

    # read the variation from the input files:
    Variation_1D_down = []
    Variation_1D_down_Eff = []
    Variation_2D_down = []

    # Get reco, migration, data with Syst File
    reco_Syst     = RecoDistribution( InputSyst, Channel, Variable)
    mig_hist_Syst = MigrationMatrix(  InputSyst, Channel, str(Energy), Variable, reco_hist)

    # 
    DireName   = Channel  + "Selection_WeightVariations"
    Hist1DName = Variable + "_Reco_v_Truth_" + str(Energy) + "TeV_cut7"
    Hist2DName = Variable + "_Reco_cut7"

    # read the directory keys:
    directory = InputSyst.GetDirectory(DireName)

    # read the variation from the file, and the migration matrix correspond to each variation
    for key in directory.GetListOfKeys():
        hist = key.ReadObj()
        if (hist.ClassName() == 'TH1D' and (hist.GetName()).find(Variable+"_") != -1 and (hist.GetName()).find("down") != -1 and  (hist.GetName()).find(Syst) != -1):
            Variation_1D_down.append(hist)
        if (hist.ClassName() == 'TH2D' and (hist.GetName()).find(Variable+"_") != -1 and (hist.GetName()).find("down") != -1 and  (hist.GetName()).find(Syst) != -1):
            Variation_2D_down.append(hist)

    # add the nominal to the variations, and apply the effeciency factor to reconstructed data
    for i in range(0, len(Variation_1D_down)):
        Variation_1D_down[i].Add(reco_hist)
        Variation_1D_down_Eff.append(ApplyEffecciency(Variation_1D_down[i], reco_hist, mig_hist))

    # add the variation to the nominal migration matrix:
    for i in range(0, len(Variation_2D_down)):
        Variation_2D_down[i].Add(mig_hist_Syst)

    # correct the reco distribution:
    reco_hist_Eff = (mig_hist.ProjectionX()).Clone("reco_hist_Eff")

    # define the nominal unfolding response matrix for Syst unfo
    response_Nominal = RooUnfoldResponse(0, 0, mig_hist, "UNFOLD", "UNFOLD")
    unfoldMCNominal  = RooUnfoldBayes(response_Nominal, reco_hist_Eff, 2)
    UnfoldingNominal = unfoldMCNominal.Hreco()

    # define the varied unfolding, with varied migration matrix:
    Systematic_list      = []
    Systematic_Cov_list  = []
    Systmatic_Uncert     = UnfoldingNominal.Clone("Systmatic_Uncert")
    Systmatic_Covariance = mig_hist_Syst.Clone("Systmatic_Covariance")

    for i in range(0, len(Variation_1D_down_Eff)):
        response_varied  = RooUnfoldResponse(0, 0, Variation_2D_down[i], "UNFOLD", "UNFOLD")
        unfoldMCVaried   = RooUnfoldBayes(response_varied, Variation_1D_down_Eff[i], 2)
        UnfoldingVaried  = unfoldMCVaried.Hreco()
        for k in range(0, UnfoldingVaried.GetNbinsX()):
            Systmatic_Uncert.SetBinContent(k+1, UnfoldingVaried.GetBinContent(k+1) - UnfoldingNominal.GetBinContent(k+1))
            Systmatic_Uncert.SetBinError(k+1, 0)
        for k in range(0, Systmatic_Covariance.GetNbinsX()):
            for l in range(0, Systmatic_Covariance.GetNbinsY()):
                Systmatic_Covariance.SetBinContent(k+1, l+1, Systmatic_Uncert.GetBinContent(k+1)*Systmatic_Uncert.GetBinContent(l+1))
        Systematic_list.append(Systmatic_Uncert)
        Systematic_Cov_list.append(Systmatic_Covariance)
    
    OutputFile = TFile("Output/output_"+Channel+str(CofM) +"/"+Variable+"/Syst_" + Syst + ".root", 'RECREATE')
    for i in range(0, len(Systematic_list)):
        Systematic_list[i].Write( "unfoldedSys_Syst" + Variation_1D_down[i].GetName())
        Systematic_Cov_list[i].Write( "unfoldedSys_Covariance_" + Variation_1D_down[i].GetName())

 
    """
    for key in directory.GetListOfKeys():
        hist = key.ReadObj()
        if ((hist.GetName()).find("Reco_cut7") != -1 and (hist.GetName()).find("sumEt") == -1 and (hist.GetName()).find("coarseHigh") == -1 and (hist.GetName()).find("eta") == -1 and (hist.GetName()).find('_down') != -1):
            if (hist.ClassName() == 'TH1D' and (hist.GetName()).find(Variable+"_") != -1 and (hist.GetName()).find(Syst) != -1):
                print(hist.GetName())
                Variation_1D_down.append(hist)
    """
    """
    for i in range(0, len(Variation_1D_down)):
        Variation_1D_down[i].Add(reco_hist)
        Variation_1D_down_Eff.append(ApplyEffecciency(Variation_1D_down[i], reco_hist, mig_hist))

    for i in range(0, len(Variation_1D_down_Eff)):
        print(i+1, Variation_1D_down_Eff[i], Variation_1D_down_Eff[i].GetNbinsX(), Variation_1D_down_Eff[i].GetMean())

    # Get reco, migration, data with Syst File
    reco_Syst     = RecoDistribution( InputSyst, Channel, Variable)
    mig_hist_Syst = MigrationMatrix(  InputSyst, Channel, str(Energy), Variable, reco_hist)

    # correct the reco distribution:
    reco_hist_Eff = (mig_hist.ProjectionX()).Clone("reco_hist_Eff")

    # define the response matrix for Syst unfo
    response_Syst = RooUnfoldResponse(0, 0, mig_hist, "UNFOLD", "UNFOLD")

    # Calculate the Syst:
    Systematics = []
    CovarianceIter = []

    for i in range(1, 1 + int(Niter)):
        # unfolded MC nominal
        unfoldMCNominal     = RooUnfoldBayes(response_Syst, reco_hist_Eff, i)
        HistoNominal        = unfoldMCNominal.Hreco()

        # Unfolded toys of variation
        unfoldedSys_down    = GetUnfoldToys(response_Syst, Variation_1D_down_Eff, i)
        Systematic_down     = GetSystematics(HistoNominal, unfoldedSys_down)
        Covariance          = GetSystCovarianceMatrix(HistoNominal, unfoldedSys_down, mig_hist, Variable)

        #Systematic_down.GetXaxis().SetRangeUser(VarMin, VarMax)
        #Covariance.GetXaxis().SetRangeUser(VarMin, VarMax)
        #Covariance.GetYaxis().SetRangeUser(VarMin, VarMax)
        #CovarianceIter.append(Covariance)
        #Systematics.append(Systematic_down)
    """
    """
    OutputFile = TFile("Output/output_"+Channel+str(CofM) +
                       "/"+Variable+"/Syst_" + Syst + ".root", 'RECREATE')
    for i in range(1, int(Niter)):
        print(Syst + "_Systematics_Iter"+str(i))
        Systematics[i].Write(Syst + "_Systematics_Iter"+str(i))
        CovarianceIter[i].Write(Syst + "_Covariance_Iter"+str(i))

    for i in range(1, int(Niter)):
        print(i, Systematics[i])
        Systematics[i].Write( Syst + "_Systematics_Iter"+str(i))
        CovarianceIter[i].Write( Syst + "_Covariance_Iter"+str(i))

    return Systematics
    """
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def GetSystCovarianceMatrix(HistoNominal, unfoldedSys_down, mig_hist, Variable):

    CovarianceTot = []

    ntbins = HistoNominal.GetNbinsX()
    xaxis = HistoNominal.GetXaxis()

    if(Variable != "WpT"):
        CovarianceMatrix = mig_hist.Clone("CovarianceMatrix")
        CovarianceSumme = mig_hist.Clone("CovarianceSumme")
        CovMatrix = mig_hist.Clone("CovMatrix")
    if(Variable == "WpT"):
        CovarianceMatrix = TH2D("CovarianceMatrix", "CovarianceMatrix", ntbins, xaxis.GetXbins(
        ).GetArray(), ntbins, xaxis.GetXbins().GetArray())
        CovarianceSumme = TH2D("CovarianceSumme", "CovarianceSumme", ntbins, xaxis.GetXbins(
        ).GetArray(), ntbins, xaxis.GetXbins().GetArray())

    for i in range(0, len(unfoldedSys_down)):
        Covariance = 0
        CovarianceMatIter = CovarianceMatrix.Clone("CovarianceMatIter")
        for j in range(1, 1 + HistoNominal.GetNbinsX()):
            for k in range(1, 1 + HistoNominal.GetNbinsX()):
                Covariance = (unfoldedSys_down[i].GetBinContent(j) - HistoNominal.GetBinContent(
                    j))*(unfoldedSys_down[i].GetBinContent(k) - HistoNominal.GetBinContent(k))
                CovarianceMatIter.SetBinContent(j, k, Covariance)
        CovarianceTot.append(CovarianceMatIter)

    for i in range(1, 1 + HistoNominal.GetNbinsX()):
        for j in range(1, 1 + HistoNominal.GetNbinsX()):
            Covariance = 0
            for k in range(0, len(CovarianceTot)):
                Covariance = Covariance + CovarianceTot[k].GetBinContent(i, j)
            CovarianceSumme.SetBinContent(i, j, Covariance)
            print("Covariance element %f" % Covariance)
    return CovarianceSumme

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def GetSystematics(HistoNominal, unfoldedSys_down):

    HistSyst = HistoNominal.Clone("HistSyst")

    for i in range(1, 1 + HistoNominal.GetNbinsX()):
        difference = 0
        for j in range(0, len(unfoldedSys_down)):
            if(HistoNominal.GetBinContent(i) != 0):
                difference = difference + pow(((unfoldedSys_down[j].GetBinContent(
                    i) - HistoNominal.GetBinContent(i)) / HistoNominal.GetBinContent(i)), 2)
        HistSyst.SetBinContent(i, 100*sqrt(difference))
        HistSyst.SetBinError(i, 0)
    return HistSyst

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def ApplyEffecciency(data_hist, reco_hist, mig_hist):

    data_hist_Corr = data_hist.Clone("data_hist_Corr")

    Efficiency_hist = (mig_hist.ProjectionX()).Clone("Efficiency_hist")
    Efficiency_hist.Divide(reco_hist)

    for i in range(1, 1 + reco_hist.GetNbinsX()):
        data_hist_Corr.SetBinContent(i, data_hist.GetBinContent(
            i)*Efficiency_hist.GetBinContent(i))
    return data_hist_Corr

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def GetCovarianceMatrix(unfoldDATA, UnfoldToys, VarMin, VarMax, Variable, mig_hist):

    ntbins = unfoldDATA.GetNbinsX()
    xaxis = unfoldDATA.GetXaxis()

    if(Variable != "WpT"):
        CovMatrix = mig_hist.Clone("CovMatrix")
    if(Variable == "WpT"):
        CovMatrix = TH2D("Covariance", "Covariance", ntbins,  xaxis.GetXbins(
        ).GetArray(), ntbins, xaxis.GetXbins().GetArray())

    for i in range(1, 1 + ntbins):
        for j in range(1, 1 + ntbins):
            MatrixValue = 0
            for k in range(0, len(UnfoldToys)):
                MatrixValue = MatrixValue + ((UnfoldToys[k].GetBinContent(i) - unfoldDATA.GetBinContent(
                    i)) * (UnfoldToys[k].GetBinContent(j) - unfoldDATA.GetBinContent(j)))
            CovMatrix.SetBinContent(i, j, MatrixValue/len(UnfoldToys))
            CovMatrix.GetXaxis().SetRangeUser(VarMin, VarMax)
            CovMatrix.GetYaxis().SetRangeUser(VarMin, VarMax)

    return CovMatrix

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def GetUnfoldToys(responseM, ToysOfData, Niter):
    UnfoldToys = []

    for i in range(0, len(ToysOfData)):
        unfoldToy = RooUnfoldBayes(responseM, ToysOfData[i], Niter)
        UnfoldToys.append((unfoldToy.Hreco()).Clone())

    return UnfoldToys

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def GetToysofData(fInput_Data_BS, Channel, Variable, Energy):
    ToysOfData = []
    director = fInput_Data_BS.GetDirectory(
        Channel+"Selection_WeightVariations")
    if(Variable == "MuEta"):
        Variable = "muEtaSF"
    #if(Variable == "elEta"): Variable = "elEtaSF"

    i = 0
    for key in director.GetListOfKeys():
        hist = key.ReadObj()
        if hist.ClassName() == 'TH1D':
            if not ((hist.GetName()).find(Variable + "_cut7_toy")):
                # print(hist.GetName())
                ToysOfData.append(hist)
                i = i+1
                if(i == 400):
                    break
    return ToysOfData

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def GetTheBias(responseM, dataCorrected, mig_hist, Niteration, VarMin, VarMax, Channel, Energy, Variable):

    # clone the nominal plots
    Migration           =  mig_hist.Clone("Response")
    truth_hist          = (mig_hist.ProjectionY()).Clone("truth_hist")
    truth_hist_Weighted = (mig_hist.ProjectionY()).Clone("truth_hist_Weighted")
    reco_hist           = (mig_hist.ProjectionX()).Clone("reco_hist")
    reco_hist_Weighted  = (mig_hist.ProjectionX()).Clone("reco_hist_Weighted")
    Response            =  mig_hist.Clone("Response")

    RatioData           = dataCorrected.Clone("RatioData")
    RatioMC             = (mig_hist.ProjectionX()).Clone("RatioMC")
    Bias                = (mig_hist.ProjectionY()).Clone("Bias")

    # fit the ratio data/MC
    RatioData.Divide(RatioMC)
    f1 = TF1("f1", "pol6", VarMin, VarMax)
    RatioData.Fit("f1")

    # reweight the truth distributions:
    for i in range(1, 1 + truth_hist.GetNbinsX()):
        binC = truth_hist.GetXaxis().GetBinCenter(i)
        print((f1.Eval(binC)))
        truth_hist_Weighted.SetBinContent(
            i, truth_hist.GetBinContent(i) * (f1.Eval(binC)))
        truth_hist_Weighted.SetBinError(i,   truth_hist.GetBinError(i))
        #print("Bin i: %d,  truth_Weighted: %f"%(i, truth_hist_Weighted.GetBinContent(i)))

    # calculate the response matrix:
    for i in range(1, 1 + reco_hist.GetNbinsX()):
        for j in range(1, 1 + truth_hist.GetNbinsX()):
            if(truth_hist.GetBinContent(j) != 0):
                Response.SetBinContent(i, j, mig_hist.GetBinContent(
                    i, j) / truth_hist.GetBinContent(j))
            if(truth_hist.GetBinContent(j) == 0):
                Response.SetBinContent(i, j, mig_hist.GetBinContent(i, j))

    # calculte the reco_weighted
    for i in range(1, 1 + reco_hist.GetNbinsX()):
        ElemV = 0
        for j in range(1, 1 + truth_hist.GetNbinsX()):
            ElemV = ElemV + \
                Response.GetBinContent(
                    i, j) * truth_hist_Weighted.GetBinContent(j)

        reco_hist_Weighted.SetBinContent(i, ElemV)
        reco_hist_Weighted.SetBinError(i,   reco_hist.GetBinError(i))

    # define the covariance matrix:
    ntbins = truth_hist_Weighted.GetNbinsX()
    xaxis = truth_hist_Weighted.GetXaxis()
    print(truth_hist.Integral(), ntbins)
    print(ntbins, xaxis, xaxis.GetXbins().GetArray())
    if(Variable != "WpT"):
        CovMatrix = mig_hist.Clone("CovMatrix")
    if(Variable == "WpT"):
        CovMatrix = TH2D("Covariance", "Covariance", ntbins,  xaxis.GetXbins(
        ).GetArray(), ntbins, xaxis.GetXbins().GetArray())
    print(" ligne 2 ")

    # define the response matrix:
    responseMN = RooUnfoldResponse(0, 0, Migration, "UNFOLD", "UNFOLD")

    # define the output:
    OutputFile = TFile("Output/output_"+Channel+str(Energy)+"/" +Variable+"/Bias_"+Channel+str(Energy)+".root", 'RECREATE')
    truth_hist.Write("truth_hist")
    truth_hist_Weighted.Write("truth_hist_Weighted")
    reco_hist.Write("reco_hist")
    reco_hist_Weighted.Write("reco_hist_Weighted")

    reco_Weighted = reco_hist_Weighted.Clone("reco_Weighted")

    # Calculate the Bias
    for i in range(1, 1 + Niteration):

        unfoldMC = RooUnfoldBayes(responseMN, reco_Weighted, i)
        UnfoldHisto = unfoldMC.Hreco()

        Bias.Add(UnfoldHisto, truth_hist_Weighted, 1, -1)

        for j in range(1, 1 + truth_hist.GetNbinsX()):
            for k in range(1, 1 + truth_hist.GetNbinsX()):
                CovMatrix.SetBinContent(
                    j, k, Bias.GetBinContent(j)*Bias.GetBinContent(k))

        for j in range(1, 1 + truth_hist.GetNbinsX()):
            if(truth_hist_Weighted.GetBinContent(j) != 0):
                Bias.SetBinContent(
                    j, 100*Bias.GetBinContent(j) / truth_hist_Weighted.GetBinContent(j))
                Bias.SetBinError(j, 0)

        Bias.GetXaxis().SetRangeUser(VarMin, VarMax)
        CovMatrix.GetXaxis().SetRangeUser(VarMin, VarMax)
        CovMatrix.GetYaxis().SetRangeUser(VarMin, VarMax)
        Bias.Write("Bias_Iter_"+str(i))
        CovMatrix.Write("CovMatrix_Iter_"+str(i))

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def BackgroundUncer(Input_Bkgd1, Input_Bkgd2, Input_Bkgd3, Input_Bkgd5, Input_Bkgd4, Channel, Var, Energy):

    Bkg_W_hist = Input_Bkgd1.Get(Channel + "Selection/" + Var + "_cut7")
    Bkg_Z_hist = Input_Bkgd2.Get(Channel + "Selection/" + Var + "_cut7")
    Bkg_D_hist = Input_Bkgd3.Get(Channel + "Selection/" + Var + "_cut7")
    Bkg_T_hist = Input_Bkgd5.Get(Channel + "Selection/" + Var + "_cut7")

    Total_Unce = Bkg_W_hist.Clone()

    for i in range(0, Bkg_W_hist.GetNbinsX()):
        Total_Unce.SetBinContent(i+1, sqrt( pow(0.05 * Bkg_W_hist.GetBinContent(i+1), 2) +
                                            pow(0.1  * Bkg_T_hist.GetBinContent(i+1), 2) +
                                            pow(0.05 * Bkg_Z_hist.GetBinContent(i+1), 2) +
                                            pow(0.1  * Bkg_D_hist.GetBinContent(i+1), 2)))
    return Total_Unce

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def LuminosityUncer(reco_hist, Input_Bkgd1, Input_Bkgd2, Input_Bkgd3, Input_Bkgd5, Input_Bkgd4, Channel, Var, Energy):

    Bkg_W_hist = Input_Bkgd1.Get(Channel + "Selection/" + Var + "_cut7")
    Bkg_Z_hist = Input_Bkgd2.Get(Channel + "Selection/" + Var + "_cut7")
    Bkg_D_hist = Input_Bkgd3.Get(Channel + "Selection/" + Var + "_cut7")
    Bkg_T_hist = Input_Bkgd5.Get(Channel + "Selection/" + Var + "_cut7")

    Total_Unce = Bkg_W_hist.Clone()

    if(Energy == "5TeV"):
        ErrForLumi = 0.016
    if(Energy == "13TeV"):
        ErrForLumi = 0.015

    for i in range(0, Bkg_W_hist.GetNbinsX()):
        Total_Unce.SetBinContent(i+1,  ErrForLumi*(reco_hist.GetBinContent(i+1) + Bkg_W_hist.GetBinContent(i+1) + Bkg_T_hist.GetBinContent(i+1) + Bkg_Z_hist.GetBinContent(i+1) + Bkg_D_hist.GetBinContent(i+1)) )

    return Total_Unce

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def MultiJetUncerta( Input_Bkgd5, Input_Bkgd4, Channel, Var):

    if "Eta"in Var:
        if "enu" in Channel:    varr = "elEta"
        if "munu" in Channel:   varr = "muEta"
    if "pT" in Var:
        if "enu" in Channel:    varr = "elPt"
        if "munu" in Channel:   varr = "muPt"

    Background_Top             = Input_Bkgd5.Get(Channel + "Selection/" + Var + "_cut7")
    Background_MjOrig          = Input_Bkgd4.Get("hist/" + varr)
    Background_shapeOrig       = Input_Bkgd4.Get("shape_uphist/" + varr)
    Background_ExtrapErrOrig   = Input_Bkgd4.Get("ExtrapErr_uphist/" + varr)
    Background_uTSliceErrOrig  = Input_Bkgd4.Get("uTSliceErr_uphist/" + varr)

    Background_Mj              = rebinTheMultiJetBackground(Background_Top.Clone(), Background_MjOrig)
    Background_shape           = rebinTheMultiJetBackground(Background_Top.Clone(), Background_shapeOrig)
    Background_ExtrapErr       = rebinTheMultiJetBackground(Background_Top.Clone(), Background_ExtrapErrOrig)
    Background_uTSliceErr      = rebinTheMultiJetBackground(Background_Top.Clone(), Background_uTSliceErrOrig)

    for i in range(0, Background_Mj.GetNbinsX()):
        Background_shape.SetBinContent(i+1,       Background_shape.GetBinContent(i)      - Background_Mj.GetBinContent(i))
        Background_ExtrapErr.SetBinContent(i+1,   Background_ExtrapErr.GetBinContent(i)  - Background_Mj.GetBinContent(i))
        Background_uTSliceErr.SetBinContent(i+1,  Background_uTSliceErr.GetBinContent(i) - Background_Mj.GetBinContent(i))

    for i in range(0, Background_shape.GetNbinsX()):
        Background_shape.SetBinContent(i+1, sqrt(pow(Background_shape.GetBinContent(i+1), 2) + pow(Background_ExtrapErr.GetBinContent(i+1), 2) + pow(Background_uTSliceErr.GetBinContent(i+1), 2)  ))

    return Background_shape

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def SumBackground(Input_Bkgd1, Input_Bkgd2, Input_Bkgd3, Input_Bkgd5, Input_Bkgd4, Channel, Var):

    if "Eta"in Var:
        if "enu" in Channel:
            varr = "hist/elEta"
        if "munu" in Channel:
            varr = "hist/muEta"
    if "pT" in Var:
        if "enu" in Channel:
            varr = "hist/elPt"
        if "munu" in Channel:
            varr = "hist/muPt"

    print("Path: ", Channel + "Selection/" + Var + "_cut7")
    Background_W        = Input_Bkgd1.Get(Channel + "Selection/" + Var + "_cut7")
    Background_Z        = Input_Bkgd2.Get(Channel + "Selection/" + Var + "_cut7")
    Background_Dilepton = Input_Bkgd3.Get(Channel + "Selection/" + Var + "_cut7")
    Background_Top      = Input_Bkgd5.Get(Channel + "Selection/" + Var + "_cut7")
    Background_MjOrig   = Input_Bkgd4.Get(varr)
    Background_Mj       = rebinTheMultiJetBackground(Background_W.Clone(), Background_MjOrig)

    print(Background_Mj.GetNbinsX())
    """
        print("Multi      bins: ",Background_Mj.GetNbinsX(),           "       Integral: ",Background_Mj.Integral())
        print("W          bins: ",Background_W.GetNbinsX(),            "       Integral: ",Background_W.Integral())
        print("Z          bins: ",Background_Z.GetNbinsX(),            "       Integral: ",Background_Z.Integral())
        print("Dilepton   bins: ",Background_Dilepton.GetNbinsX(),     "       Integral: ",Background_Dilepton.Integral())
        print("Top        bins: ",Background_Top.GetNbinsX(),          "       Integral: ",Background_Top.Integral())
	"""

    Background_Total = Background_W.Clone("Background_Total")
    for i in range(0, Background_Total.GetNbinsX()):
        Background_Total.SetBinContent(i+1, Background_W.GetBinContent(i+1) + Background_Z.GetBinContent(i+1) + Background_Dilepton.GetBinContent(i+1) + Background_Top.GetBinContent(i+1) + Background_Mj.GetBinContent(i+1))

    print("Sum the background ditributions: Done")
    return Background_Total

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def MigrationMatrix(fInput_MC, Channel, Energy, Var, truth_hist):

    migration = fInput_MC.Get(Channel + "Selection/"+Var+"_Reco_v_Truth_cut7")

    print("Get the migration Matrix: Done")
    return migration

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def TruthDistribution(fInput_MC, Var, Energy):
    truth_hist = fInput_MC.Get("TruthSelection/"+Var+"_Truth_cut4")
    print("Get the Truth Distribution: Done")
    return truth_hist

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def DataDistribution(fInput_MC, Channel, Var):
    reco_hist = fInput_MC.Get(Channel + "Selection/"+Var+"_Reco_cut7")
    print("Get the Reco Distribution: Done")
    return reco_hist

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def RecoDistribution(fInput_MC, Channel, Var):
    reco_hist = fInput_MC.Get(Channel + "Selection/"+Var+"_Reco_cut7")
    print("Get the Reco Distribution: Done")
    return reco_hist

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def GetAcceptance(mig_hist, truth_hist):
    Acceptance_hist = (mig_hist.ProjectionY()).Clone("Acceptance_hist")
    Acceptance_hist.Divide(truth_hist)
    print("Get the Acceptance Corrections: Done")
    return Acceptance_hist

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def GetEfficieny(mig_hist, reco_hist):
    Efficiency_hist = (mig_hist.ProjectionX()).Clone("Efficiency_hist")
    Efficiency_hist.Divide(reco_hist)
    print("Get the Efficiency Corrections: Done")
    return Efficiency_hist

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def CorrectData(data_hist, reco_hist, Background_Total, Efficiency_hist):

    dataCorrected = data_hist.Clone("dataCorrected")

    for i in range(1, 1 + reco_hist.GetNbinsX()):
        if(reco_hist.GetBinContent(i) != 0):
            rapportBGMC = (Background_Total.GetBinContent(
                i) / (reco_hist.GetBinContent(i)+Background_Total.GetBinContent(i)))
            dataCorrected.SetBinContent(
                i, data_hist.GetBinContent(i)*(1-rapportBGMC))

    for i in range(1, 1 + dataCorrected.GetNbinsX()):
        dataCorrected.SetBinContent(i, dataCorrected.GetBinContent(
            i)*Efficiency_hist.GetBinContent(i))

    print("Correct data and subtract background: Done")

    return dataCorrected
