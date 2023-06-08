#!/usr/bin/env python
# -*-coding:Latin-1 -*
import ROOT
import ROOT as root
import matplotlib.pyplot as plt
import atlasplots as aplt
import numpy as np
import pylab
import sys
import os
import time
import gc

from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TH1, TF1, TText
from math import *
from array import array
from ROOT import gStyle
import numpy as np

from ROOT import TCanvas, TGraph
from ROOT import gROOT
from math import sin

from array import array
from ROOT import gROOT, TCanvas, TFile, TLatex, THStack, TH1F, TPad, TLine, TAttFill, TMatrixD, THilbertMatrixD, TDecompSVD, TGraphErrors, gRandom, TText

def makeLegend(hists, xmin, ymin, xmax, ymax):
    legend = root.TLegend(xmin, ymin, xmax, ymax)
    legend.SetTextSize(0.05)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetBorderSize(0)
    for hist in hists:
        legend.AddEntry(hist, hist.GetName())
    return legend

# ----------------------------------------------------------------------------------------------------------------------------------

def StatStudy(inputFile, IterMin, IterMax, mTmin, mTmax, channel, Var):

    binContent  = []    # some bin content
    StatError   = []    # with correlation

    Range1      = []
    for ite in range(IterMin, IterMax):
        CovMatrix        = inputFile.Get("CovarianceMatrix_Data_Iter"+ str(ite))
        data_Unfolded    = inputFile.Get("unfolded_data"+ str(ite))

        Error = 0
        Value = 0

        jmin = data_Unfolded.FindBin(mTmin)
        jmax = data_Unfolded.FindBin(mTmax)
        for j in range(jmin, jmax):
            Value = Value + data_Unfolded.GetBinContent(j)
            for k in range(jmin, jmax):
                Error = Error + CovMatrix[j,k]

        print(Value, sqrt(Error), 100*sqrt(Error)/Value)
        Range1.append(Error)
        Range1.append(Value)

    return Range1

# ----------------------------------------------------------------------------------------------------------------------------------

def BiasStudy(Bias, inputFile, IterMin, IterMax, mTmin, mTmax, channel, Var):

    binContent  = []    # some bin content
    StatError   = []    # with correlation

    Range1      = []
    for ite in range(IterMin, IterMax):
        CovMatrix        = Bias.Get("CovMatrix_Iter_"+ str(ite))
        data_Unfolded    = inputFile.Get("unfolded_data"+ str(ite))

        Error = 0
        Value = 0

        jmin = data_Unfolded.FindBin(mTmin)
        jmax = data_Unfolded.FindBin(mTmax)
        for j in range(jmin, jmax):
            Value = Value + data_Unfolded.GetBinContent(j)
            for k in range(jmin, jmax):
                Error = Error + CovMatrix.GetBinContent(j,k)

        print(Value, sqrt(Error), 100*sqrt(Error)/Value)
        Range1.append(Error)
        Range1.append(Value)

    return Range1

# ----------------------------------------------------------------------------------------------------------------------------------

def SystStudy(IterMin, IterMax, mTmin, mTmax, Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum, Variable):

    binContent  = []    # some bin content
    StatError   = []    # with correlation

    Range1      = []
    for ite in range(IterMin, IterMax):

        data_Unfolded    = Summarize.Get("unfolded_data"+ str(ite))

        # define systematics:
        if  (Name.find("enu")  != -1):
            IDCovMatrix      = IdSF.Get(    "ElIDSys_Covariance_Iter"+ str(ite))
            TrigCovMatrix    = TrigSF.Get(  "TrigSys_Covariance_Iter"+ str(ite))
            RecoCovMatrix    = RecoSF.Get(  "RecoSys_Covariance_Iter"+ str(ite))
            IsoCovMatrix     = IsoSF.Get(   "IsoSys_Covariance_Iter"+  str(ite))
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter"+  str(ite))
            CalibCovMatrix   = Calib.Get(   "Calib_Covariance_Iter"+   str(ite))
        if  (Name.find("munu")  != -1):
            TrigCovMatrix    = TrigSF.Get(  "TrigSys_Covariance_Iter"+ str(ite))
            RecoCovMatrix    = RecoSF.Get(  "RecoSys_Covariance_Iter"+ str(ite))
            IsoCovMatrix     = IsoSF.Get(   "IsoSys_Covariance_Iter"+  str(ite))
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter"+  str(ite))

        CovMatrix   = RecoCovMatrix.Clone("CovMatrix")

        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                if  (Name.find("enu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) )
                if  (Name.find("munu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )


        Error = 0
        Value = 0

        jmin = data_Unfolded.FindBin(mTmin)
        jmax = data_Unfolded.FindBin(mTmax)
        for j in range(jmin, jmax):
            Value = Value + data_Unfolded.GetBinContent(j)
            for k in range(jmin, jmax):
                Error = Error + CovMatrix.GetBinContent(j,k)

        print(Value, sqrt(Error), 100*sqrt(Error)/Value)
        Range1.append(Error)
        Range1.append(Value)

    return Range1

# ----------------------------------------------------------------------------------------------------------------------------------

def PlotOptimisationStat(Statrange1, Statrange2, Statrange3, channel, Energy, Varibale):

    # Data Cooedinates
    x = np.arange(1, 6)

    relative1 = [] 
    relative2 = [] 
    relative3 = [] 

    for i in [0, 2, 4, 6, 8]:
        relative1.append(100.*sqrt(Statrange1[i])/Statrange1[i+1])
        relative2.append(100.*sqrt(Statrange2[i])/Statrange2[i+1])
        relative3.append(100.*sqrt(Statrange3[i])/Statrange3[i+1])

    # PLot
    if (Varibale == "Lepton_pT"):
        plt.plot(x,np.array(relative1), label="range [25, 60] GeV") 
        plt.plot(x,np.array(relative2), label="range [30, 50] GeV")  
        plt.plot(x,np.array(relative3), label="range [35, 45] GeV") 
    if (Varibale == "Lepton_Eta"):
        plt.plot(x,np.array(relative1), label=" range [-1.05, 1.05]") 
        plt.plot(x,np.array(relative2), label=" range [-1.36, 1.36]") 
        plt.plot(x,np.array(relative3), label=" range [-2.50, 2.50]")  

    plt.legend()


    # Add Axes Labels
    plt.xlabel("Number of Iterations") 
    plt.ylabel("Statistical uncertainty (%)") 

    # Display
    plt.savefig("Output/Optimisation/Optimisation_Stat_"+channel+"_"+Energy+"_"+Varibale+".pdf")
    plt.close()

# ----------------------------------------------------------------------------------------------------------------------------------

def PlotOptimisationBias(Statrange1, Statrange2, Statrange3, channel, Energy, Varibale):

    # Data Cooedinates
    x = np.arange(1, 6)

    relative1 = [] 
    relative2 = [] 
    relative3 = [] 

    for i in [0, 2, 4, 6, 8]:
        relative1.append(100.*sqrt(Statrange1[i])/Statrange1[i+1])
        relative2.append(100.*sqrt(Statrange2[i])/Statrange2[i+1])
        relative3.append(100.*sqrt(Statrange3[i])/Statrange3[i+1])

    # PLot
    if (Varibale == "Lepton_pT"):
        plt.plot(x,np.array(relative1), label="range [25, 60] GeV") 
        plt.plot(x,np.array(relative2), label="range [30, 50] GeV")  
        plt.plot(x,np.array(relative3), label="range [35, 45] GeV") 
    if (Varibale == "Lepton_Eta"):
        plt.plot(x,np.array(relative1), label=" range [-1.05, 1.05]") 
        plt.plot(x,np.array(relative2), label=" range [-1.36, 1.36]") 
        plt.plot(x,np.array(relative3), label=" range [-2.50, 2.50]")  
    plt.legend()


    # Add Axes Labels
    plt.xlabel("Number of Iterations") 
    plt.ylabel("Bias uncertainty (%)") 

    # Display
    plt.savefig("Output/Optimisation/Optimisation_Bias_"+channel+"_"+Energy+"_"+Varibale+".pdf")
    plt.close()

# ----------------------------------------------------------------------------------------------------------------------------------

def PlotOptimisationSyst(Statrange1, Statrange2, Statrange3, channel, Energy, Varibale):

    # Data Cooedinates
    x = np.arange(1, 6)

    relative1 = [] 
    relative2 = [] 
    relative3 = [] 

    for i in [0, 2, 4, 6, 8]:
        relative1.append(100.*sqrt(Statrange1[i])/Statrange1[i+1])
        relative2.append(100.*sqrt(Statrange2[i])/Statrange2[i+1])
        relative3.append(100.*sqrt(Statrange3[i])/Statrange3[i+1])

    # PLot
    if (Varibale == "Lepton_pT"):
        plt.plot(x,np.array(relative1), label="range [25, 60] GeV") 
        plt.plot(x,np.array(relative2), label="range [30, 50] GeV")  
        plt.plot(x,np.array(relative3), label="range [35, 45] GeV") 
    if (Varibale == "Lepton_Eta"):
        plt.plot(x,np.array(relative1), label=" range [-1.05, 1.05]") 
        plt.plot(x,np.array(relative2), label=" range [-1.36, 1.36]") 
        plt.plot(x,np.array(relative3), label=" range [-2.50, 2.50]")  
    plt.legend()


    # Add Axes Labels
    plt.xlabel("Number of Iterations") 
    plt.ylabel("Systematic uncertainty (%)") 

    # Display
    plt.savefig("Output/Optimisation/Optimisation_Syst_"+channel+"_"+Energy+"_"+Varibale+".pdf")
    plt.close()

# ----------------------------------------------------------------------------------------------------------------------------------

def PlotOptimisationtot(Statrange1, Statrange2, Statrange3, Biasrange1, Biasrange2, Biasrange3, Systerange1, Systerange2, Systerange3, channel, Energy, Varibale, Indice):

    # Data Cooedinates
    x = np.arange(1, 6)

    relative1 = [] 
    relative2 = [] 
    relative3 = [] 

    for i in [0, 2, 4, 6, 8]:
        relative1.append( 100. * sqrt(  Statrange1[i] + Biasrange1[i] + Systerange1[i]  )  /  Statrange1[i+1]  )
        relative2.append( 100. * sqrt(  Statrange2[i] + Biasrange2[i] + Systerange2[i]  )  /  Statrange2[i+1]  )
        relative3.append( 100. * sqrt(  Statrange3[i] + Biasrange3[i] + Systerange3[i]  )  /  Statrange3[i+1]  )

    maxx = max(np.amax(np.array(relative1)), np.amax(np.array(relative1)), np.amax(np.array(relative1)))

    # PLot
    if (Varibale == "Lepton_pT"):
        plt.plot(x,np.array(relative1), label="range [25, 60] GeV") 
        plt.plot(x,np.array(relative2), label="range [30, 50] GeV")  
        plt.plot(x,np.array(relative3), label="range [35, 45] GeV") 
    if (Varibale == "Lepton_Eta"):
        plt.plot(x,np.array(relative1), label=" range [-1.05, 1.05]") 
        plt.plot(x,np.array(relative2), label=" range [-1.36, 1.36]") 
        plt.plot(x,np.array(relative3), label=" range [-2.50, 2.50]")  
    plt.legend()


    # Add Axes Labels
    plt.xlabel("Number of Iterations") 
    plt.ylabel("Total uncertainty (%)") 

    """
    font = {'family': 'serif','weight': 'normal','size': 12,}
    if( Varibale == "Lepton_pT"):
        plt.text(1, maxx, r'ATLAS Internal', fontdict=font)
        plt.text(1, maxx-0.02, r'$ '+Indice+'$', fontdict=font)
    if( Varibale == "Lepton_Eta"):
        plt.text(1, maxx, r'ATLAS Internal', fontdict=font)
        plt.text(1, maxx-0.02, r'$ '+Indice+'$', fontdict=font)
    """
    print(maxx)
    # Display
    plt.savefig("Output/Optimisation/Optimisation_total_"+channel+"_"+Energy+"_"+Varibale+".pdf")
    plt.close()

# ----------------------------------------------------------------------------------------------------------------------------------

def GetDiffernetialXsMu(Summarize, Bias, TrigSF, RecoSF, IsoSF, Recoil, Energy, Indice, Name, Lum, Variable):

    data = Summarize.Get("data_hist")
    reco = Summarize.Get("reco_hist")
    truth = Summarize.Get("truth_hist")
    mig_mat = Summarize.Get("mig_hist")
    bg = Summarize.Get("Background_Total")

    Hist1 = truth.Clone("Hist1")
    Hist2 = truth.Clone("Hist2")
    Hist3 = truth.Clone("Hist3")

    Acceptance = Summarize.Get("Acceptance_hist")

    # define Unfolded:
    HUnfolded = Summarize.Get("unfolded_data4")

    # Stat:
    Covariance = Summarize.Get("CovarianceMatrix_Data_Iter4")
    #print(Covariance, Covariance[26,26])

    # Bias:
    CovBias = Bias.Get("CovMatrix_Iter_1")

    # define systematics:
    if (Name.find("enu") != -1):

        IDCovMatrix = IdSF.Get("ElIDSys_Covariance_Iter4")
        TrigCovMatrix = TrigSF.Get("TrigSys_Covariance_Iter4")
        RecoCovMatrix = RecoSF.Get("RecoSys_Covariance_Iter4")
        IsoCovMatrix = IsoSF.Get("IsoSys_Covariance_Iter4")
        RecoilCovMatrix = Recoil.Get("Recoil_Covariance_Iter4")
        CalibCovMatrix = Calib.Get("Calib_Covariance_Iter4")

    if (Name.find("munu") != -1):

        TrigCovMatrix = TrigSF.Get("TrigSys_Covariance_Iter4")
        RecoCovMatrix = RecoSF.Get("RecoSys_Covariance_Iter4")
        IsoCovMatrix = IsoSF.Get("IsoSys_Covariance_Iter4")
        RecoilCovMatrix = Recoil.Get("Recoil_Covariance_Iter4")

    CovMatrix = RecoCovMatrix.Clone("CovMatrix")

    for i in range(1, 1 + CovMatrix.GetNbinsX()):
        for j in range(1, 1+CovMatrix.GetNbinsX()):
            if (Name.find("enu") != -1):
                CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i, j) + RecoCovMatrix.GetBinContent(i, j) + IsoCovMatrix.GetBinContent(
                    i, j) + IDCovMatrix.GetBinContent(i, j) + CalibCovMatrix.GetBinContent(i, j) + RecoilCovMatrix.GetBinContent(i, j))
            if (Name.find("munu") != -1):
                CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(
                    i, j) + RecoCovMatrix.GetBinContent(i, j) + RecoilCovMatrix.GetBinContent(i, j))

    for i in range(0, Hist1.GetNbinsX()):
        if (HUnfolded.GetBinContent(i+1) != 0):
            Hist1.SetBinContent(
                i+1, 100*sqrt(Covariance[i, i]) / HUnfolded.GetBinContent(i+1))
            Hist2.SetBinContent(
                i+1, 100*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1))
            Hist3.SetBinContent(
                i+1, 100*sqrt(CovBias.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1))
            Hist1.SetBinError(i+1, 0)
            Hist2.SetBinError(i+1, 0)
            Hist3.SetBinError(i+1, 0)
        if (HUnfolded.GetBinContent(i+1) == 0):
            Hist1.SetBinContent(i+1, 0)
            Hist2.SetBinContent(i+1, 0)
            Hist3.SetBinContent(i+1, 0)
            Hist1.SetBinError(i+1, 0)
            Hist2.SetBinError(i+1, 0)
            Hist3.SetBinError(i+1, 0)

    latexFile = open("Output/LatexTableau/Differential_" +
                     Variable+"Xs_"+Name+".tex", "w+")
    latexFile.write("\\documentclass[12pt]{article} \n")
    latexFile.write("\\usepackage{amsmath}  \n")
    latexFile.write("\\usepackage{graphicx} \n")
    latexFile.write("\\usepackage{hyperref} \n")
    latexFile.write("\\usepackage{hyperref} \n")
    latexFile.write("\\usepackage[latin1]{inputenc} \n")
    latexFile.write("\\begin{document} \n")

    latexFile.write("\\begin{table}[] \n")
    latexFile.write("\\begin{tabular}{c|c|c|c|c|}\n")
    latexFile.write("\\cline{2-5}\n")

    latexFile.write(
        "& \\multicolumn{4}{c|}{%s}  \\\ \\hline \\hline \n" % (Indice))

    if ("pT" in Variable):
        latexFile.write(
            "\\multicolumn{1}{|c|}{  range } & $d\sigma$/$dp^{T}_{l}$ [GeV]     & Stat uncertainty     & Unfolding bias     & Syst uncertainty        \\\ \\hline \\hline \n")
        i = 25
        max = 50
    if ("Eta" in Variable):
        latexFile.write(
            "\\multicolumn{1}{|c|}{  range } & $d\sigma$/$d\eta_{l}$     & Stat uncertainty     & Unfolding bias     & Syst uncertainty        \\\ \\hline \\hline \n")
        i = 0
        max = HUnfolded.GetNbinsX()

    while i < max:
        if (HUnfolded.GetBinWidth(i+1) != 0 and Acceptance.GetBinContent(i+1) != 0):
            diffXs = HUnfolded.GetBinContent(
                i+1)/(Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Stat = 100*sqrt(Covariance[i, i]) / HUnfolded.GetBinContent(i+1)
            Bias = 100*sqrt(CovBias.GetBinContent(i+1, i+1)) / \
                HUnfolded.GetBinContent(i+1)
            Syst = 100*sqrt(CovMatrix.GetBinContent(i+1, i+1)
                            ) / HUnfolded.GetBinContent(i+1)
        if (HUnfolded.GetBinWidth(i+1) == 0 or Acceptance.GetBinContent(i+1) == 0):
            diffXs = 0
            Stat = 0
            Bias = 0
            Syst = 0
        lowedge = HUnfolded.GetBinLowEdge(i+1)
        highEdge = HUnfolded.GetBinLowEdge(i+2)

        latexFile.write("\\multicolumn{1}{|c|}{{[}%3.2f,  %3.2f{]}}  & %5.3f & %5.3f & %5.3f & %5.3f \\\ \\hline \n" % (
            lowedge, highEdge, diffXs, Stat, Bias, Syst))

        i = i+1

    latexFile.write("\\end{tabular}\n")
    latexFile.write("\\end{table}\n")
    latexFile.write("\\end{document}\n")
    latexFile.close()

# ----------------------------------------------------------------------------------------------------------------------------------

def GetFiducialXsEl(Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum, Variable):
        print("Summarising Table")
        data       = Summarize.Get("data_hist")
        reco       = Summarize.Get("reco_hist")
        truth      = Summarize.Get("truth_hist")
        mig_mat    = Summarize.Get("mig_hist")
        bg         = Summarize.Get("Background_Total")
        
        # **************************************** Bin-by-Bin resultats ****************************************

        # Cross section
        Xsbyb = ( (data.Integral()-bg.Integral())/Lum ) * ( truth.Integral() / reco.Integral() )
    
        # Stat uncertainties
        Statbyb = 100/sqrt(data.Integral())
    
        # **************************************** Unfolding resultats ****************************************

        # define Unfolded:
        HUnfolded  = Summarize.Get("unfolded_data1")
    
        # Stat:        
        Covariance = Summarize.Get("CovarianceMatrix_Data_Iter1")

        Stat = 0
        for i in range(0,  HUnfolded.GetNbinsX()):
            for j in range(0, HUnfolded.GetNbinsX()):
                Stat = Stat + Covariance[i, j]
        
        # Bias:
        CovBias      = Bias.Get("CovMatrix_Iter_1")
            
        Bias = 0
        for i in range(1, 1 + CovBias.GetNbinsX()):
            for j in range(1, 1 + CovBias.GetNbinsX()):
                Bias = Bias + CovBias.GetBinContent(i, j)

        # define systematics:
        if  (Name.find("enu")  != -1):

            IDCovMatrix      = IdSF.Get(    "ElIDSys_Covariance_Iter1"  )
            TrigCovMatrix    = TrigSF.Get(  "TrigSys_Covariance_Iter1"  )
            RecoCovMatrix    = RecoSF.Get(  "RecoSys_Covariance_Iter1"  )
            IsoCovMatrix     = IsoSF.Get(   "IsoSys_Covariance_Iter1" )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter1"   )
            CalibCovMatrix   = Calib.Get(   "Calib_Covariance_Iter1"    )
        if  (Name.find("munu")  != -1):
            TrigCovMatrix    = TrigSF.Get(  "TrigSys_Covariance_Iter1" )
            RecoCovMatrix    = RecoSF.Get(  "RecoSys_Covariance_Iter1" )
            IsoCovMatrix     = IsoSF.Get(   "IsoSys_Covariance_Iter1"  )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter1"    )


        CovMatrix   = RecoCovMatrix.Clone("CovMatrix")

        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                if  (Name.find("enu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) )
                if  (Name.find("munu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )

        Syst = 0
        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                Syst = Syst + CovMatrix.GetBinContent(i, j)
    
        # define the cross section:
        truth     = Summarize.Get("truth_hist")
        truth_Acc = (Summarize.Get("mig_hist")).ProjectionY()

        Xs = ( HUnfolded.Integral()/Lum ) * ( truth.Integral() / truth_Acc.Integral() )
        Syst = ( truth.Integral() / truth_Acc.Integral() ) * sqrt(Syst) /  Lum
        Bias = ( truth.Integral() / truth_Acc.Integral() ) * sqrt(Bias) /  Lum
        Stat =  ( truth.Integral() / truth_Acc.Integral() ) * sqrt(Stat) /  Lum
        LumError =  1.6 * HUnfolded.Integral()/100
        LumXsErr = ( truth.Integral() / truth_Acc.Integral() ) * LumError /  Lum

        print("Cross sections : ", Xs)
        print("Nombre d'events : ", HUnfolded.Integral())
        print("Acceptance : ", truth.Integral() / truth_Acc.Integral() )
        print("Stat Uncert : ", Stat)
        print("Syst Uncert : ", Syst)
        print("Bias Uncert : ", Bias)
        print("luminosity error on Nunfold : ",   LumError)
        print("luminosity error on Xs      : ",   LumXsErr)
        print("Summarize table for differential Xs : ")

        
        #astyle.SetAtlasStyle()
        c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )
        c1.SetTickx()

        Syst1 = Stat
        Syst2 = sqrt(Stat*Stat + Syst*Syst)        
        Syst3 = sqrt(Stat*Stat + Syst*Syst + LumXsErr*LumXsErr)

        n = 1;

        print(Xs)
        x  = array( 'f', [      Xs] )
        ex = array( 'f', [   Syst1] )
        y  = array( 'f', [      Xs] )
        ey = array( 'f', [    1000] )
        gr = TGraphErrors( n, x, y, ex, ey )

        xx  = array( 'f', [     Xs] )
        exx = array( 'f', [  Syst2] )
        yy  = array( 'f', [     Xs] )
        eyy = array( 'f', [   1000] )
        grr = TGraphErrors( n, xx, yy, exx, eyy )

        xxx  = array( 'f', [    Xs] )
        exxx = array( 'f', [ Syst3] )
        yyy  = array( 'f', [    Xs] )
        eyyy = array( 'f', [  1000] )
        grrr = TGraphErrors( n, xxx, yyy, exxx, eyyy )

        grrr.SetTitle( 'TGraphErrors Example' )
        grrr.SetTitle("")
        axis = grrr.GetXaxis();
        axis.SetLimits(Xs-200,Xs+70)
        grrr.GetHistogram().SetMinimum(Xs-150);
        grrr.GetHistogram().SetMaximum(Xs+150);
        grrr.GetYaxis().SetBinLabel(1, "")
        grrr.GetXaxis().SetTitle("#sigma^{fid}")

        gr.SetMarkerColor( 5 )
        gr.SetMarkerStyle( 21 )
        gr.SetFillColor(5)
        gr.SetFillStyle(1001)

        grr.SetMarkerColor( 3 )
        grr.SetMarkerStyle( 21 )
        grr.SetFillColor(3)
        grr.SetFillStyle(1001)

        grrr.SetMarkerColor( 7 )
        grrr.SetMarkerStyle( 21 )
        grrr.SetFillColor(7)
        grrr.SetFillStyle(1001)

        grrr.Draw("ap2")
        grr.Draw("CP2")
        gr.Draw("CP2")
        #ax.plot(grrr, label="Random Hist", labelfmt="F")
        #ax.plot(grr, label="Random Hist", labelfmt="F")
        #ax.plot(gr, label="Random Hist", labelfmt="F")


        line = root.TLine(Xs, Xs-150, Xs, Xs+150)
        line.SetLineWidth(2)

        # define the legend
        histsN = []
        gr.SetLineWidth(0)
        grr.SetLineWidth(0)
        grrr.SetLineWidth(0)
        gr.SetName("data \pm stat")
        grr.SetName("data \pm stat \pm Syst")
        grrr.SetName("data \pm stat \pm Syst \pm Lum")

        latex = TLatex()
        latex.SetTextSize(0.065)
        latex.SetTextAlign(13)
        latex.DrawLatex(Xs-180,Xs+130,"ATLAS Internal")
        latex.DrawLatex(Xs-180,Xs+105,Indice)

        th1x  = array( 'f', [      Xs-5] )
        th1ex = array( 'f', [   5] )
        th1y  = array( 'f', [      Xs] )
        th1ey = array( 'f', [    0.3] )
        Therogr1 = TGraphErrors( n, th1x, th1y, th1ex, th1ey )
        Therogr1.SetName("CT18NNLO")
        Therogr1.SetMarkerStyle(107)
        Therogr1.SetMarkerSize(1)

        th2x  = array( 'f', [      Xs+53] )
        th2ex = array( 'f', [   6.3] )
        th2y  = array( 'f', [      Xs-40] )
        th2ey = array( 'f', [    0.3] )
        Therogr2 = TGraphErrors( n, th2x, th2y, th2ex, th2ey )
        Therogr2.SetName("HERAPDF20NNLO")
        Therogr2.SetMarkerStyle(73)
        Therogr2.SetMarkerSize(1)

        th3x  = array( 'f', [      Xs-73] )
        th3ex = array( 'f', [   9.3] )
        th3y  = array( 'f', [      Xs-80] )
        th3ey = array( 'f', [    0.3] )
        Therogr3 = TGraphErrors( n, th3x, th3y, th3ex, th3ey )
        Therogr3.SetName("NNPDF31_nnlo")
        Therogr3.SetMarkerStyle(69)
        Therogr3.SetMarkerSize(1)


        histsN.append(gr)
        histsN.append(grr)
        histsN.append(grrr)
        histsN.append(Therogr1)
        histsN.append(Therogr2)
        histsN.append(Therogr3)
        legendN = makeLegend(histsN,0.12, 0.23,0.45,0.68)
        legendN.Draw("same")
        Therogr1.Draw("same P")
        Therogr2.Draw("same P")
        Therogr3.Draw("same P")
        line.Draw("same")

        # TCanvas.Update() draws the frame, after which one can change it
        c1.Update()
        c1.Print("Output/CrossSection/Fiducial_CrossSection_"+Name+".pdf")
        
# ----------------------------------------------------------------------------------------------------------------------------------

def GetFiducialXsMu(Summarize, Bias, TrigSF, RecoSF, IsoSF, Recoil, Energy, Indice, Name, Lum, Variable):
              
        HBias      = Bias.Get("Bias_Iter_1")
        HUnfolded  = Summarize.Get("unfolded_data1")

        data       = Summarize.Get("data_hist")
        reco       = Summarize.Get("reco_hist")
        truth      = Summarize.Get("truth_hist")

        HStatError = HUnfolded.Clone("HStatError")
        HSystTotal = HUnfolded.Clone("HSystTotal")
        Covariance = Summarize.Get("CovarianceMatrix_Data_Iter1")


    
        # Stat:        
        Covariance = Summarize.Get("CovarianceMatrix_Data_Iter1")

        Stat = 0
        for i in range(0,  HUnfolded.GetNbinsX()):
            for j in range(0, HUnfolded.GetNbinsX()):
                Stat = Stat + Covariance[i, j]
        

        
        if  (Name.find("enu")  != -1):            
            IDCovMatrix      = IdSF.Get(    "ElIDSys_Covariance_Iter1"  )   
            TrigCovMatrix    = TrigSF.Get(  "TrigSys_Covariance_Iter1"  )  
            RecoCovMatrix    = RecoSF.Get(  "RecoSys_Covariance_Iter1"  ) 
            IsoCovMatrix     = IsoSF.Get(   "IsoSys_Covariance_Iter1" ) 
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter1"   ) 
            CalibCovMatrix   = Calib.Get(   "Calib_Covariance_Iter1"    ) 
        
        if  (Name.find("munu")  != -1):
            TrigCovMatrix    = TrigSF.Get(  "TrigSys_Covariance_Iter1" ) 
            RecoCovMatrix    = RecoSF.Get(  "RecoSys_Covariance_Iter1" ) 
            IsoCovMatrix     = IsoSF.Get(   "IsoSys_Covariance_Iter1"  )  
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter1"   ) 
                    
        CovMatrix   = TrigCovMatrix.Clone("CovMatrix")
                    
        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                if  (Name.find("enu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) )
                if  (Name.find("munu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )
        Syst = 0
        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                Syst = Syst + CovMatrix.GetBinContent(i, j)
    
        # Bias:
        CovBias      = Bias.Get("CovMatrix_Iter_1")
            
        Bias = 0
        for i in range(1, 1 + CovBias.GetNbinsX()):
            for j in range(1, 1 + CovBias.GetNbinsX()):
                Bias = Bias + CovBias.GetBinContent(i, j)

        # define the cross section:
        truth     = Summarize.Get("truth_hist")
        truth_Acc = (Summarize.Get("mig_hist")).ProjectionY()

        Xs = ( HUnfolded.Integral()/Lum ) * ( truth.Integral() / truth_Acc.Integral() )
        Syst = ( truth.Integral() / truth_Acc.Integral() ) * sqrt(Syst) /  Lum
        Bias = ( truth.Integral() / truth_Acc.Integral() ) * sqrt(Bias) /  Lum
        Stat =  ( truth.Integral() / truth_Acc.Integral() ) * sqrt(Stat) /  Lum
        LumError =  1.6 * HUnfolded.Integral()/100
        LumXsErr = ( truth.Integral() / truth_Acc.Integral() ) * LumError /  Lum

        print("Cross sections : ", Xs)
        print("Nombre d'events : ", HUnfolded.Integral())
        print("Acceptance : ", truth.Integral() / truth_Acc.Integral() )
        print("Stat Uncert : ", Stat)
        print("Syst Uncert : ", Syst)
        print("Bias Uncert : ", Bias)
        print("luminosity error on Nunfold : ",   LumError)
        print("luminosity error on Xs      : ",   LumXsErr)
        print("Summarize table for differential Xs : ")

        
        #astyle.SetAtlasStyle()
        c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )
        c1.SetTickx()

        Syst1 = Stat
        Syst2 = sqrt(Stat*Stat + Syst*Syst)        
        Syst3 = sqrt(Stat*Stat + Syst*Syst + LumXsErr*LumXsErr)

        n = 1;

        print(Xs)
        x  = array( 'f', [      Xs] )
        ex = array( 'f', [   Syst1] )
        y  = array( 'f', [      Xs] )
        ey = array( 'f', [    1000] )
        gr = TGraphErrors( n, x, y, ex, ey )

        xx  = array( 'f', [     Xs] )
        exx = array( 'f', [  Syst2] )
        yy  = array( 'f', [     Xs] )
        eyy = array( 'f', [   1000] )
        grr = TGraphErrors( n, xx, yy, exx, eyy )

        xxx  = array( 'f', [    Xs] )
        exxx = array( 'f', [ Syst3] )
        yyy  = array( 'f', [    Xs] )
        eyyy = array( 'f', [  1000] )
        grrr = TGraphErrors( n, xxx, yyy, exxx, eyyy )

        grrr.SetTitle( 'TGraphErrors Example' )
        grrr.SetTitle("")
        axis = grrr.GetXaxis();
        axis.SetLimits(Xs-200,Xs+70)
        grrr.GetHistogram().SetMinimum(Xs-150);
        grrr.GetHistogram().SetMaximum(Xs+150);
        grrr.GetYaxis().SetBinLabel(1, "")
        grrr.GetXaxis().SetTitle("#sigma^{fid}")

        gr.SetMarkerColor( 5 )
        gr.SetMarkerStyle( 21 )
        gr.SetFillColor(5)
        gr.SetFillStyle(1001)

        grr.SetMarkerColor( 3 )
        grr.SetMarkerStyle( 21 )
        grr.SetFillColor(3)
        grr.SetFillStyle(1001)

        grrr.SetMarkerColor( 7 )
        grrr.SetMarkerStyle( 21 )
        grrr.SetFillColor(7)
        grrr.SetFillStyle(1001)

        grrr.Draw("ap2")
        grr.Draw("CP2")
        gr.Draw("CP2")
        #ax.plot(grrr, label="Random Hist", labelfmt="F")
        #ax.plot(grr, label="Random Hist", labelfmt="F")
        #ax.plot(gr, label="Random Hist", labelfmt="F")


        line = root.TLine(Xs, Xs-150, Xs, Xs+150)
        line.SetLineWidth(2)

        # define the legend
        histsN = []
        gr.SetLineWidth(0)
        grr.SetLineWidth(0)
        grrr.SetLineWidth(0)
        gr.SetName("data \pm stat")
        grr.SetName("data \pm stat \pm Syst")
        grrr.SetName("data \pm stat \pm Syst \pm Lum")

        latex = TLatex()
        latex.SetTextSize(0.065)
        latex.SetTextAlign(13)
        latex.DrawLatex(Xs-180,Xs+130,"ATLAS Internal")
        latex.DrawLatex(Xs-180,Xs+105,Indice)

        th1x  = array( 'f', [      Xs-5] )
        th1ex = array( 'f', [   5] )
        th1y  = array( 'f', [      Xs] )
        th1ey = array( 'f', [    0.3] )
        Therogr1 = TGraphErrors( n, th1x, th1y, th1ex, th1ey )
        Therogr1.SetName("CT18NNLO")
        Therogr1.SetMarkerStyle(107)
        Therogr1.SetMarkerSize(1)

        th2x  = array( 'f', [      Xs+53] )
        th2ex = array( 'f', [   6.3] )
        th2y  = array( 'f', [      Xs-40] )
        th2ey = array( 'f', [    0.3] )
        Therogr2 = TGraphErrors( n, th2x, th2y, th2ex, th2ey )
        Therogr2.SetName("HERAPDF20NNLO")
        Therogr2.SetMarkerStyle(73)
        Therogr2.SetMarkerSize(1)

        th3x  = array( 'f', [      Xs-73] )
        th3ex = array( 'f', [   9.3] )
        th3y  = array( 'f', [      Xs-80] )
        th3ey = array( 'f', [    0.3] )
        Therogr3 = TGraphErrors( n, th3x, th3y, th3ex, th3ey )
        Therogr3.SetName("NNPDF31_nnlo")
        Therogr3.SetMarkerStyle(69)
        Therogr3.SetMarkerSize(1)


        histsN.append(gr)
        histsN.append(grr)
        histsN.append(grrr)
        histsN.append(Therogr1)
        histsN.append(Therogr2)
        histsN.append(Therogr3)
        legendN = makeLegend(histsN,0.12, 0.23,0.45,0.68)
        legendN.Draw("same")
        Therogr1.Draw("same P")
        Therogr2.Draw("same P")
        Therogr3.Draw("same P")
        line.Draw("same")

        # TCanvas.Update() draws the frame, after which one can change it
        c1.Update()
        c1.Print("Output/CrossSection/Fiducial_CrossSection_"+Name+".pdf")
        

# ----------------------------------------------------------------------------------------------------------------------------------

def GetSummaringTableEl(Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum, Variable):
        print("Summarising Table")
        data       = Summarize.Get("data_hist")
        reco       = Summarize.Get("reco_hist")
        truth      = Summarize.Get("truth_hist")
        mig_mat    = Summarize.Get("mig_hist")
        bg         = Summarize.Get("Background_Total")
        
        # **************************************** Bin-by-Bin resultats ****************************************

        # Cross section
        Xsbyb = ( (data.Integral()-bg.Integral())/Lum ) * ( truth.Integral() / reco.Integral() )
    
        # Stat uncertainties
        Statbyb = 100/sqrt(data.Integral())
    
        # **************************************** Unfolding resultats ****************************************

        # define Unfolded:
        HUnfolded  = Summarize.Get("unfolded_data1")
    
        # Stat:        
        Covariance = Summarize.Get("CovarianceMatrix_Data_Iter1")

        Stat = 0
        for i in range(0,  HUnfolded.GetNbinsX()):
            for j in range(0, HUnfolded.GetNbinsX()):
                Stat = Stat + Covariance[i, j]
        
        # Bias:
        CovBias      = Bias.Get("CovMatrix_Iter_1")
            
        Bias = 0
        for i in range(1, 1 + CovBias.GetNbinsX()):
            for j in range(1, 1 + CovBias.GetNbinsX()):
                Bias = Bias + CovBias.GetBinContent(i, j)

        # define systematics:
        if  (Name.find("enu")  != -1):

            IDCovMatrix      = IdSF.Get(    "ElIDSys_Covariance_Iter1"  )
            TrigCovMatrix    = TrigSF.Get(  "TrigSys_Covariance_Iter1"  )
            RecoCovMatrix    = RecoSF.Get(  "RecoSys_Covariance_Iter1"  )
            IsoCovMatrix     = IsoSF.Get(   "IsoSys_Covariance_Iter1" )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter1"   )
            CalibCovMatrix   = Calib.Get(   "Calib_Covariance_Iter1"    )
        if  (Name.find("munu")  != -1):
            TrigCovMatrix    = TrigSF.Get(  "TrigSys_Covariance_Iter1" )
            RecoCovMatrix    = RecoSF.Get(  "RecoSys_Covariance_Iter1" )
            IsoCovMatrix     = IsoSF.Get(   "IsoSys_Covariance_Iter1"  )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter1"    )


        CovMatrix   = RecoCovMatrix.Clone("CovMatrix")

        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                if  (Name.find("enu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) )
                if  (Name.find("munu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )

        Syst = 0
        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                Syst = Syst + CovMatrix.GetBinContent(i, j)
    
        # define the cross section:
        truth     = Summarize.Get("truth_hist")
        truth_Acc = (Summarize.Get("mig_hist")).ProjectionY()

        Xs = ( HUnfolded.Integral()/Lum ) * ( truth.Integral() / truth_Acc.Integral() )
        Syst = ( truth.Integral() / truth_Acc.Integral() ) * sqrt(Syst) /  Lum
        Bias = ( truth.Integral() / truth_Acc.Integral() ) * sqrt(Bias) /  Lum
        Stat =  ( truth.Integral() / truth_Acc.Integral() ) * sqrt(Stat) /  Lum
        LumError =  1.6 * HUnfolded.Integral()/100
        LumXsErr = ( truth.Integral() / truth_Acc.Integral() ) * LumError /  Lum

        print("Cross sections : ", Xs)
        print("Nombre d'events : ", HUnfolded.Integral())
        print("Acceptance : ", truth.Integral() / truth_Acc.Integral() )

        print("Stat Uncert : ", Stat)
        print("Syst Uncert : ", Syst)
        print("Bias Uncert : ", Bias)

        print("luminosity error on Nunfold : ",   LumError)
        print("luminosity error on Xs      : ",   LumXsErr)
        print("Summarize table for differential Xs : ")

        latexFile = open("Output/LatexTableau/FiducialCross_Section_"+Name+".tex","w+")
        latexFile.write("\\documentclass[12pt]{article} \n")
        latexFile.write("\\usepackage{amsmath}\n")
        latexFile.write("\\usepackage{graphicx}\n")
        latexFile.write("\\usepackage{hyperref}\n")
        latexFile.write("\\usepackage{hyperref}\n")
        latexFile.write("\\usepackage[latin1]{inputenc}\n")
        latexFile.write("\\begin{document}\n")

        latexFile.write("\\begin{table}[ht]\n")
        latexFile.write("\\begin{tabular}{c|c|}\n")
        latexFile.write("\\cline{2-2}\n")
        latexFile.write("                                                           &    %s  \\\ \\hline \n"%(Indice+",  \\text{(value $\\pm$ stat $\\pm$ syst $\\pm$ lum) [pb]}"))
        latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ }                     &    %5.3f   $\\pm$ %3.2f $\\pm$ %3.2f  $\\pm$ %3.2f     \\\ \\hline \n" %( Xs,    Stat,    Syst,      LumXsErr))
        #latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ (Bin-by-Bin)}        &    %5.3f   $\\pm$ %3.2f(Stat) $\\pm$ %5.3f(Syst)  $\\pm$ %5.3f(Lum)     \\\ \\hline \n" %( Xsbyb, Statbyb, Systbyb,   LumXsErr))
        latexFile.write("\\end{tabular}\n")
        latexFile.write("\\end{table}\n")
        latexFile.write("\\end{document}\n")
        latexFile.close()
# ----------------------------------------------------------------------------------------------------------------------------------

def GetSummaringTableMu(Summarize, Bias, TrigSF, RecoSF, IsoSF,  Recoil, Energy, Indice, Name, Lum, Variable):
        print("Summarising Table")
        data       = Summarize.Get("data_hist")
        reco       = Summarize.Get("reco_hist")
        truth      = Summarize.Get("truth_hist")
        mig_mat    = Summarize.Get("mig_hist")
        bg         = Summarize.Get("Background_Total")
        
        # **************************************** Bin-by-Bin resultats ****************************************

        # Cross section
        Xsbyb = ( (data.Integral()-bg.Integral())/Lum ) * ( truth.Integral() / reco.Integral() )
    
        # Stat uncertainties
        Statbyb = 100/sqrt(data.Integral())
    
        # **************************************** Unfolding resultats ****************************************

        # define Unfolded:
        HUnfolded  = Summarize.Get("unfolded_data1")
    
        # Stat:        
        Covariance = Summarize.Get("CovarianceMatrix_Data_Iter1")

        Stat = 0
        for i in range(0,  HUnfolded.GetNbinsX()):
            for j in range(0, HUnfolded.GetNbinsX()):
                Stat = Stat + Covariance[i, j]
            
        # Bias:
        CovBias      = Bias.Get("CovMatrix_Iter_1")
            
        Bias = 0
        for i in range(1, 1 + CovBias.GetNbinsX()):
            for j in range(1, 1 + CovBias.GetNbinsX()):
                Bias = Bias + CovBias.GetBinContent(i, j)

        # define systematics:
        if  (Name.find("enu")  != -1):

            IDCovMatrix      = IdSF.Get(    "ElIDSys_Covariance_Iter1"  )
            TrigCovMatrix    = TrigSF.Get(  "TrigSys_Covariance_Iter1"  )
            RecoCovMatrix    = RecoSF.Get(  "RecoSys_Covariance_Iter1"  )
            IsoCovMatrix     = IsoSF.Get(   "IsoSys_Covariance_Iter1" )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter1"   )
            CalibCovMatrix   = Calib.Get(   "Calib_Covariance_Iter1"    )
        if  (Name.find("munu")  != -1):
            TrigCovMatrix    = TrigSF.Get(  "TrigSys_Covariance_Iter1" )
            RecoCovMatrix    = RecoSF.Get(  "RecoSys_Covariance_Iter1" )
            IsoCovMatrix     = IsoSF.Get(   "IsoSys_Covariance_Iter1"  )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter1"    )


        CovMatrix   = RecoCovMatrix.Clone("CovMatrix")

        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                if  (Name.find("enu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) )
                if  (Name.find("munu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )

        Syst = 0
        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                Syst = Syst + CovMatrix.GetBinContent(i, j)
    

        # define the cross section:
        truth     = Summarize.Get("truth_hist")
        truth_Acc = (Summarize.Get("mig_hist")).ProjectionY()
        Xs = ( HUnfolded.Integral()/Lum ) * ( truth.Integral() / truth_Acc.Integral() )
        Syst = ( truth.Integral() / truth_Acc.Integral() ) * sqrt(Syst) /  Lum
        Bias = ( truth.Integral() / truth_Acc.Integral() ) * sqrt(Bias) /  Lum
        Stat =  ( truth.Integral() / truth_Acc.Integral() ) * sqrt(Stat) /  Lum
        LumError =  1.6 * HUnfolded.Integral()/100
        LumXsErr = ( truth.Integral() / truth_Acc.Integral() ) * LumError /  Lum

        print("Cross sections : ", Xs)
        print("Nombre d'events : ", HUnfolded.Integral())
        print("Acceptance : ", truth.Integral() / truth_Acc.Integral() )

        print("Stat Uncert : ", Stat)
        print("Syst Uncert : ", Syst)
        print("Bias Uncert : ", Bias)

        print("luminosity error on Nunfold : ",   LumError)
        print("luminosity error on Xs      : ",   LumXsErr)
        print("Summarize table for differential Xs : ")

        latexFile = open("Output/LatexTableau/FiducialCross_Section_"+Name+".tex","w+")
        latexFile.write("\\documentclass[12pt]{article} \n")
        latexFile.write("\\usepackage{amsmath}\n")
        latexFile.write("\\usepackage{graphicx}\n")
        latexFile.write("\\usepackage{hyperref}\n")
        latexFile.write("\\usepackage{hyperref}\n")
        latexFile.write("\\usepackage[latin1]{inputenc}\n")
        latexFile.write("\\begin{document}\n")

        latexFile.write("\\begin{table}[ht]\n")
        latexFile.write("\\begin{tabular}{c|c|}\n")
        latexFile.write("\\cline{2-2}\n")
        latexFile.write("                                                           &    %s  \\\ \\hline \n"%(Indice+",  \\text{(value $\\pm$ stat $\\pm$ syst $\\pm$ lum) [pb]}"))
        latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ }                     &    %5.3f   $\\pm$ %3.2f $\\pm$ %3.2f  $\\pm$ %3.2f     \\\ \\hline \n" %( Xs,    Stat,    Syst,      LumXsErr))
        #latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ (Bin-by-Bin)}        &    %5.3f   $\\pm$ %3.2f(Stat) $\\pm$ %5.3f(Syst)  $\\pm$ %5.3f(Lum)     \\\ \\hline \n" %( Xsbyb, Statbyb, Systbyb,   LumXsErr))
        latexFile.write("\\end{tabular}\n")
        latexFile.write("\\end{table}\n")
        latexFile.write("\\end{document}\n")
        latexFile.close()
