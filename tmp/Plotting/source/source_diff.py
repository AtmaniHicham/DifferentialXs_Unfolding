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
# ----------------------------------------------------------------------------------------------------------------------------------

def GetDiffernetialXsEl(Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum, Variable):

    print("********************* start *********************")

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
    Covariance = Summarize.Get("CovarianceMatrix_Data_Iter2")
    #print(Covariance, Covariance[26,26])

    # Bias:
    CovBias = Bias.Get("CovMatrix_Iter_2")

    # define systematics:
    if (Name.find("enu") != -1):

        IDCovMatrix = IdSF.Get("ElIDSys_Covariance_Iter2")
        TrigCovMatrix = TrigSF.Get("TrigSys_Covariance_Iter2")
        RecoCovMatrix = RecoSF.Get("RecoSys_Covariance_Iter2")
        IsoCovMatrix = IsoSF.Get("IsoSys_Covariance_Iter2")
        RecoilCovMatrix = Recoil.Get("Recoil_Covariance_Iter2")
        CalibCovMatrix = Calib.Get("Calib_Covariance_Iter2")

    if (Name.find("munu") != -1):

        TrigCovMatrix = TrigSF.Get("MuTrigSys_Covariance_Iter2")
        RecoCovMatrix = RecoSF.Get("MuRecoSys_Covariance_Iter2")
        IsoCovMatrix = IsoSF.Get("MuIsoSys_Covariance_Iter2")
        RecoilCovMatrix = Recoil.Get("Recoil_Covariance_Iter2s")

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
        i = 25
        max = 50
        latexFile.write(
            "\\multicolumn{1}{|c|}{  range } & $d\sigma$/$dp^{T}_{l}$ [GeV]     & Stat uncertainty     & Unfolding bias     & Syst uncertainty        \\\ \\hline \\hline \n")

    if ("Eta" in Variable):
        i = 0
        max = HUnfolded.GetNbinsX()
        latexFile.write(
            "\\multicolumn{1}{|c|}{  range } & $d\sigma$/$d\eta_{l}$      & Stat uncertainty     & Unfolding bias     & Syst uncertainty        \\\ \\hline \\hline \n")

    while i < max:
        if (HUnfolded.GetBinWidth(i+1) != 0 and Acceptance.GetBinContent(i+1) != 0):
            diffXs = HUnfolded.GetBinContent(i+1)/(Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Stat   = (100./diffXs)*sqrt(Covariance[i, i])                  / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Bias   = (100./diffXs)*sqrt(CovBias.GetBinContent(i+1, i+1))   / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Syst   = (100./diffXs)*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            #Stat = 100*sqrt(Covariance[i, i]) / HUnfolded.GetBinContent(i+1)
            #Bias = 100*sqrt(CovBias.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1)
            #Syst = 100*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1)
        if (HUnfolded.GetBinWidth(i+1) == 0 or Acceptance.GetBinContent(i+1) == 0):
            diffXs = 0
            Stat = 0
            Bias = 0
            Syst = 0
        lowedge = HUnfolded.GetBinLowEdge(i+1)
        highEdge = HUnfolded.GetBinLowEdge(i+2)
        print(lowedge, highEdge)
        latexFile.write("\\multicolumn{1}{|c|}{{[}%3.2f,  %3.2f{]}}  & %5.3f & %5.3f & %5.3f & %5.3f \\\ \\hline \n" % (
            lowedge, highEdge, diffXs, Stat, Bias, Syst))

        i = i+1

    latexFile.write("\\end{tabular}\n")
    latexFile.write("\\end{table}\n")
    latexFile.write("\\end{document}\n")
    latexFile.close()

# ----------------------------------------------------------------------------------------------------------------------------------
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
            diffXs = HUnfolded.GetBinContent(i+1)/(Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Stat   = (100./diffXs)*sqrt(Covariance[i, i])                  / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Bias   = (100./diffXs)*sqrt(CovBias.GetBinContent(i+1, i+1))   / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Syst   = (100./diffXs)*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
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
# ----------------------------------------------------------------------------------------------------------------------------------


def GetDiffernetialXsElPlot(Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum, Variable):

    print("********************* start *********************")

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
    Covariance = Summarize.Get("CovarianceMatrix_Data_Iter2")
    #print(Covariance, Covariance[26,26])

    # Bias:
    CovBias = Bias.Get("CovMatrix_Iter_2")

    # define systematics:
    if (Name.find("enu") != -1):

        IDCovMatrix = IdSF.Get("ElIDSys_Covariance_Iter2")
        TrigCovMatrix = TrigSF.Get("TrigSys_Covariance_Iter2")
        RecoCovMatrix = RecoSF.Get("RecoSys_Covariance_Iter2")
        IsoCovMatrix = IsoSF.Get("IsoSys_Covariance_Iter2")
        RecoilCovMatrix = Recoil.Get("Recoil_Covariance_Iter2")
        CalibCovMatrix = Calib.Get("Calib_Covariance_Iter2")

    if (Name.find("munu") != -1):

        TrigCovMatrix = TrigSF.Get("MuTrigSys_Covariance_Iter2")
        RecoCovMatrix = RecoSF.Get("MuRecoSys_Covariance_Iter2")
        IsoCovMatrix = IsoSF.Get("MuIsoSys_Covariance_Iter2")
        RecoilCovMatrix = Recoil.Get("Recoil_Covariance_Iter2s")

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
        i = 25
        max = 50
        latexFile.write(
            "\\multicolumn{1}{|c|}{  range } & $d\sigma$/$dp^{T}_{l}$ [GeV]     & Stat uncertainty     & Unfolding bias     & Syst uncertainty        \\\ \\hline \\hline \n")

    if ("Eta" in Variable):
        i = 0
        max = HUnfolded.GetNbinsX()
        latexFile.write(
            "\\multicolumn{1}{|c|}{  range } & $d\sigma$/$d\eta_{l}$      & Stat uncertainty     & Unfolding bias     & Syst uncertainty        \\\ \\hline \\hline \n")

    while i < max:
        if (HUnfolded.GetBinWidth(i+1) != 0 and Acceptance.GetBinContent(i+1) != 0):
            diffXs = HUnfolded.GetBinContent(i+1)/(Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Stat   = (100./diffXs)*sqrt(Covariance[i, i])                  / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Bias   = (100./diffXs)*sqrt(CovBias.GetBinContent(i+1, i+1))   / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Syst   = (100./diffXs)*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            print(Stat, diffXs)
            #Stat = 100*sqrt(Covariance[i, i]) / HUnfolded.GetBinContent(i+1)
            #Bias = 100*sqrt(CovBias.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1)
            #Syst = 100*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1)
        if (HUnfolded.GetBinWidth(i+1) == 0 or Acceptance.GetBinContent(i+1) == 0):
            diffXs = 0
            Stat = 0
            Bias = 0
            Syst = 0
        lowedge = HUnfolded.GetBinLowEdge(i+1)
        highEdge = HUnfolded.GetBinLowEdge(i+2)

        latexFile.write("\\multicolumn{1}{|c|}{{[}%2d,  %2d{]}}  & %5.3f & %5.3f & %5.3f & %5.3f \\\ \\hline \n" % (
            lowedge, highEdge, diffXs, Stat, Bias, Syst))

        i = i+1

    latexFile.write("\\end{tabular}\n")
    latexFile.write("\\end{table}\n")
    latexFile.write("\\end{document}\n")
    latexFile.close()

# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
