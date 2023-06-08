#!/usr/bin/env python
# -*-coding:Latin-1 -*
import ROOT
import ROOT as root
import matplotlib.pyplot as plt
import atlasplots as aplt
import seaborn as sns
import numpy as np
import pylab
import sys
import os
import time
import gc

from atlasify import atlasify
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size': 14})

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

def FunctionConparisonNormalised(CrossSection, TruthHisto, indice, Name, Energy, Variable, truth_sherpa, truth_PhH7EG):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, (ax1, ax2) = aplt.ratio_plot(name="fig1", figsize=(1200, 800), hspace=0.05)

    # define the variables:
    if "pT" in Variable:
        ax1.set_yscale("log")
    #if "Eta" in Variable:
        #ax1.set_yscale("log")

    # Randomly fill two histograms according to the above distribution
    hist1 = CrossSection.Clone("hist1")
    hist2 = TruthHisto.Clone("hist2")
    hist3 = truth_sherpa.Clone("hist3")
    hist4 = truth_PhH7EG.Clone("hist4")

    # The line width
    hist1.SetLineWidth(2)
    hist2.SetLineWidth(2)
    hist1.SetMarkerSize(2)

    hist3.SetLineStyle(2)
    hist4.SetLineStyle(2)

    hist3.SetLineWidth(2)
    hist4.SetLineWidth(2)
    hist3.SetMarkerSize(2)

    # Draw the histograms on these axes
    hist_graph = aplt.root_helpers.hist_to_graph(hist1)
    ax1.plot(hist2, linecolor=2  )
    ax1.plot(hist3, linecolor=8  )
    ax1.plot(hist4, linecolor=4  )
    ax1.plot(hist_graph,"P", linecolor=1 )

    # Calculate and draw the ratio
    ratio_hist = hist1.Clone("ratio_hist")
    ratio_hist.Divide(hist2)

    ratio_hist1 = hist1.Clone("ratio_hist1")
    ratio_hist1.Divide(hist3)

    ratio_hist2 = hist1.Clone("ratio_hist2")
    ratio_hist2.Divide(hist4)

    #for i in range(0, ratio_hist.GetNbinsX()):
    #    ratio_hist.SetBinError(i+1, 0.1)
    ratio_graph  = aplt.root_helpers.hist_to_graph(ratio_hist)
    ratio_graph1 = aplt.root_helpers.hist_to_graph(ratio_hist1)
    ratio_graph2 = aplt.root_helpers.hist_to_graph(ratio_hist2)

    ax2.plot(ratio_graph,  "P", linecolor=2)
    ax2.plot(ratio_graph1, "P", linecolor=8)
    ax2.plot(ratio_graph2, "P", linecolor=4)
    ratio_graph.SetMarkerColor(2)
    ratio_graph1.SetMarkerColor(8)
    ratio_graph2.SetMarkerColor(4)

    ratio_graph1.SetLineStyle(2)
    ratio_graph2.SetLineStyle(2)

    # Plot the relative error on the ratio axes
    err_band_ratio = aplt.root_helpers.hist_to_graph(
            ratio_hist,
            show_bin_width=True,
            norm=True)
    ax2.plot(err_band_ratio, "2", fillcolor=1, fillstyle=3254)


    # Add extra space at top of plot to make room for labels
    ax1.add_margins(top=0.16)
    ax2.add_margins(top=0.1, bottom=0.3)

    if "pT" in Variable:
        # Add extra space at top and bottom of ratio panel
        min = 0
        max = 1
        for i in range(0, hist1.GetNbinsX()):
            if i == 0: 
                min = hist1.GetBinContent(i+1)    
                max = hist1.GetBinContent(i+1) 
            if hist1.GetBinContent(i+1) < min: min = hist1.GetBinContent(i+1)
            if hist1.GetBinContent(i+1) > max: max = hist1.GetBinContent(i+1)
        ax1.set_ylim(0.001, max+4)
        ax1.set_xlim(25, 60)
        ax2.set_xlim(25, 60)
        ax2.set_ylim(0.82, 1.18)
        # Set axis titles
        ax2.set_ylabel("Pred./Data", loc="centre",  titleoffset=1.8)
        ax1.set_ylabel("1/#sigma  d#sigma / dp_{T}^{l} [GeV^{-1}]", titleoffset=1.8)
        ax2.set_xlabel(" p_{T}^{l} [GeV]", titleoffset=2.8)

    if "Eta" in Variable:
        # Add extra space at top and bottom of ratio panel
        min = hist1.GetBinContent(1)
        max = hist1.GetBinContent(1)
        for i in range(0, hist1.GetNbinsX()):
            if hist1.GetBinContent(i+1) < min: min = hist1.GetBinContent(i+1)
            if hist1.GetBinContent(i+1) > max: max = hist1.GetBinContent(i+1)
        ax1.set_ylim(min-0.01, max+0.014)
        ax1.set_xlim(-2.50, 2.50)
        ax2.set_xlim(-2.50, 2.50)
        ax2.set_ylim(0.82, 1.18)
        # Set axis titles
        ax2.set_ylabel("Pred./Data", loc="centre",  titleoffset=1.4)
        ax1.set_ylabel("1/#sigma  d#sigma / d#eta_{l} ",      titleoffset=1.4)
        ax2.set_xlabel(" #eta_{l}",                 titleoffset=2.9)

    # Draw line at y=1 in ratio panel
    line = root.TLine(ax1.get_xlim()[0], 1, ax1.get_xlim()[1], 1)
    ax2.plot(line)

    # Go back to top axes to add labels
    ax1.cd()

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax1.text(0.2, 0.84, indice, size=22, align=13)

    # Add legend
    #legend = ax1.legend(loc=(0.7, 0.78, 0.94, 0.90))
    legend = ax1.legend(loc=(0.62, 0.54, 1 - root.gPad.GetRightMargin() - 0.03, 1 - root.gPad.GetTopMargin() - 0.03),textsize=22)
    legend.AddEntry(hist_graph,     "Unfolded Data",    "EP")
    legend.AddEntry(hist2,          "PowhegPythia8",    "L")
    legend.AddEntry(hist3,          "Sherpa",           "L")
    legend.AddEntry(hist4,          "PhH7EG",           "L")
    legend.AddEntry(err_band_ratio, "Total Unc.",       "F")

    # Save the plot as a PDF
    fig.savefig("Output/DiffCrossSection/Normalised_Differential_" + Name+"_Xs_"+Energy+"_"+Variable+".pdf")


def FunctionConparison(CrossSection, TruthHisto, indice, Name, Energy, Variable, truth_sherpa, truth_PhH7EG):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, (ax1, ax2) = aplt.ratio_plot(name="fig1", figsize=(1200, 800), hspace=0.05)

    # define the variables:
    if "pT" in Variable:
        ax1.set_yscale("log")
    #if "Eta" in Variable:
        #ax1.set_yscale("log")

    # Randomly fill two histograms according to the above distribution
    hist1 = CrossSection.Clone("hist1")
    hist2 = TruthHisto.Clone("hist2")
    hist3 = truth_sherpa.Clone("hist3")
    hist4 = truth_PhH7EG.Clone("hist4")

    # The line width
    hist1.SetLineWidth(2)
    hist2.SetLineWidth(2)
    hist1.SetMarkerSize(2)

    hist3.SetLineStyle(2)
    hist4.SetLineStyle(2)

    hist3.SetLineWidth(2)
    hist4.SetLineWidth(2)
    hist3.SetMarkerSize(2)

    # Draw the histograms on these axes
    hist_graph = aplt.root_helpers.hist_to_graph(hist1)
    ax1.plot(hist2, linecolor=2  )
    ax1.plot(hist3, linecolor=8  )
    ax1.plot(hist4, linecolor=4  )
    ax1.plot(hist_graph,"P", linecolor=1 )

    # Calculate and draw the ratio
    ratio_hist = hist1.Clone("ratio_hist")
    ratio_hist.Divide(hist2)

    ratio_hist1 = hist1.Clone("ratio_hist1")
    ratio_hist1.Divide(hist3)

    ratio_hist2 = hist1.Clone("ratio_hist2")
    ratio_hist2.Divide(hist4)

    #for i in range(0, ratio_hist.GetNbinsX()):
    #    ratio_hist.SetBinError(i+1, 0.1)
    ratio_graph  = aplt.root_helpers.hist_to_graph(ratio_hist)
    ratio_graph1 = aplt.root_helpers.hist_to_graph(ratio_hist1)
    ratio_graph2 = aplt.root_helpers.hist_to_graph(ratio_hist2)

    ax2.plot(ratio_graph,  "P", linecolor=2)
    ax2.plot(ratio_graph1, "P", linecolor=8)
    ax2.plot(ratio_graph2, "P", linecolor=4)
    ratio_graph.SetMarkerColor(2)
    ratio_graph1.SetMarkerColor(8)
    ratio_graph2.SetMarkerColor(4)

    ratio_graph1.SetLineStyle(2)
    ratio_graph2.SetLineStyle(2)

    # Plot the relative error on the ratio axes
    err_band_ratio = aplt.root_helpers.hist_to_graph(
            ratio_hist,
            show_bin_width=True,
            norm=True)
    ax2.plot(err_band_ratio, "2", fillcolor=1, fillstyle=3254)


    # Add extra space at top of plot to make room for labels
    ax1.add_margins(top=0.16)
    ax2.add_margins(top=0.1, bottom=0.3)

    if "pT" in Variable:
        # Add extra space at top and bottom of ratio panel
        min = 0
        max = 1
        for i in range(0, hist1.GetNbinsX()):
            if i == 0: 
                min = hist1.GetBinContent(i+1)    
                max = hist1.GetBinContent(i+1) 
            if hist1.GetBinContent(i+1) < min: min = hist1.GetBinContent(i+1)
            if hist1.GetBinContent(i+1) > max: max = hist1.GetBinContent(i+1)
        ax1.set_ylim(1, max+3000)
        ax1.set_xlim(25, 60)
        ax2.set_xlim(25, 60)
        ax2.set_ylim(0.82, 1.18)
        # Set axis titles
        ax2.set_ylabel("Pred./Data", loc="centre",  titleoffset=1.8)
        ax1.set_ylabel("d#sigma / dp_{T}^{l} [GeV^{-1}]", titleoffset=1.8)
        ax2.set_xlabel(" p_{T}^{l} [GeV]", titleoffset=2.8)

    if "Eta" in Variable:
        # Add extra space at top and bottom of ratio panel
        min = hist1.GetBinContent(1)  
        max = hist1.GetBinContent(1)  
        for i in range(0, hist1.GetNbinsX()):
            if hist1.GetBinContent(i+1) < min: min = hist1.GetBinContent(i+1)
            if hist1.GetBinContent(i+1) > max: max = hist1.GetBinContent(i+1)
        ax1.set_ylim(min-200, max+300)
        ax1.set_xlim(-2.50, 2.50)
        ax2.set_xlim(-2.50, 2.50)
        ax2.set_ylim(0.82, 1.18)
        # Set axis titles
        ax2.set_ylabel("Pred./Data", loc="centre",  titleoffset=1.4)
        ax1.set_ylabel("d#sigma / d#eta_{l} ",      titleoffset=1.4)
        ax2.set_xlabel(" #eta_{l}",                 titleoffset=2.9)

    # Draw line at y=1 in ratio panel
    line = root.TLine(ax1.get_xlim()[0], 1, ax1.get_xlim()[1], 1)
    ax2.plot(line)

    # Go back to top axes to add labels
    ax1.cd()

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax1.text(0.2, 0.84, indice, size=22, align=13)

    # Add legend
    #legend = ax1.legend(loc=(0.7, 0.78, 0.94, 0.90))
    legend = ax1.legend(loc=(0.62, 0.54, 1 - root.gPad.GetRightMargin() - 0.03, 1 - root.gPad.GetTopMargin() - 0.03),textsize=22)
    legend.AddEntry(hist_graph,     "Unfolded Data",    "EP")
    legend.AddEntry(hist2,          "PowhegPythia8",    "L")
    legend.AddEntry(hist3,          "Sherpa",           "L")
    legend.AddEntry(hist4,          "PhH7EG",           "L")
    legend.AddEntry(err_band_ratio, "Total Unc.",       "F")

    # Save the plot as a PDF
    fig.savefig("Output/DiffCrossSection/Differential_" + Name+"_Xs_"+Energy+"_"+Variable+".pdf")

def plotingfunction(ArrayDiff, ArrayTruth, Uncert, arrayBinCenter, xname, yname, indice, Name):

    #bin_counts, bin_edges, patches = plt.hist(arrayBinCenter, bins=[0, 10, 20, 30, 40, 50, 100])
    #bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2

    fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [2, 1]})
    atlasify("Internal",    axes=axs[0])
    atlasify("",            axes=axs[1])

    #axs[0].hist([ArrayDiff, ArrayTruth], arrayBinCenter, histtype='bar', stacked=True)
    axs[0].plot(ArrayDiff, arrayBinCenter)

    ratio = ArrayDiff/ArrayTruth
    axs[1].plot(arrayBinCenter,ratio)
    axs[0].set_ylabel('Events',     fontsize=14)
    axs[1].set_xlabel('$m_{T}$',    fontsize=14)
    axs[1].set_ylabel('ratio',      fontsize=14)
    #axs[1].set_ylim([0, 2])

    #plt.xlabel('$m_{T}$', fontsize=16)
    #plt.ylabel('Events', fontsize=16)

    plt.fill_between(arrayBinCenter, ratio - 0.1, ratio + 0.1, color='b', alpha=0.2)
    plt.savefig("plot.pdf")

def plotingratio(hist1, hist2, xname, yname, indice, Name):

    for i in range(0, hist1.GetNbinsX()):
        print(hist1.GetBinContent(i+1), hist2.GetBinContent(i+1))
   # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, (ax1, ax2) = aplt.ratio_plot(name="fig1", figsize=(800, 800), hspace=0.05)

    # Draw the histograms on these axes
    ax1.plot(hist1, "EP X0",  linewidth=1,  linecolor=root.kRed+1,  label="Red",  labelfmt="EP")
    ax1.plot(hist2,                         linecolor=1, label="Blue", labelfmt="EP")

    # Calculate and draw the ratio
    ratio_hist = hist1.Clone("ratio_hist")
    ratio_hist.Divide(hist2)
    ax2.plot(ratio_hist, linecolor=root.kBlack)

    # Add extra space at top of plot to make room for labels
    ax1.add_margins(top=0.16)

    # Set axis titles
    ax2.set_xlabel(xname)
    ax1.set_ylabel(yname)
    ax1.set_yscale("log")
    ax1.set_ylim(0.8, 300)
    ax1.set_xlim(25, 70)
    ax2.set_xlim(25, 70)
    ax2.set_ylim(0.7, 1.3)
    ax2.set_ylabel("Ratio", loc="centre")

    # Draw line at y=1 in ratio panel
    line = root.TLine(ax1.get_xlim()[0], 1, ax1.get_xlim()[1], 1)
    ax2.plot(line)

    # Add extra space at top and bottom of ratio panel
    ax2.add_margins(top=0.1, bottom=0.1)

    # Go back to top axes to add labels
    ax1.cd()

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax1.text(0.7, 0.84, indice, size=22, align=13)

    # Add legend
    ax1.legend(loc=(0.58, 0.78, 73, 0.90))

    # Save the plot as a PDF
    fig.savefig(Name+".pdf")


def GetDiffernetialXs(Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum, Variable):

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

    Lumiunc = 0

    # Bias:
    CovBias = Bias.Get("CovMatrix_Iter_2")

    # define systematics:
    if (Name.find("enu") != -1):
        IDCovMatrix     = IdSF.Get("ElIDSys_Covariance_Iter2")
        TrigCovMatrix   = TrigSF.Get("TrigSys_Covariance_Iter2")
        RecoCovMatrix   = RecoSF.Get("RecoSys_Covariance_Iter2")
        IsoCovMatrix    = IsoSF.Get("IsoSys_Covariance_Iter2")
        RecoilCovMatrix = Recoil.Get("Recoil_Covariance_Iter2")
        CalibCovMatrix  = Calib.Get("Calib_Covariance_Iter2")

    if (Name.find("munu") != -1):

        TrigCovMatrix   = TrigSF.Get("TrigSys_Covariance_Iter2")
        RecoCovMatrix   = RecoSF.Get("RecoSys_Covariance_Iter2")
        IsoCovMatrix    = IsoSF.Get("IsoSys_Covariance_Iter2")
        RecoilCovMatrix = Recoil.Get("Recoil_Covariance_Iter2s")

    CovMatrix = RecoCovMatrix.Clone("CovMatrix")

    for i in range(1, 1 + CovMatrix.GetNbinsX()):
        for j in range(1, 1+CovMatrix.GetNbinsX()):
            if (Name.find("enu") != -1):
                CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i, j) + RecoCovMatrix.GetBinContent(i, j) + IsoCovMatrix.GetBinContent(i, j) + IDCovMatrix.GetBinContent(i, j) + CalibCovMatrix.GetBinContent(i, j) + RecoilCovMatrix.GetBinContent(i, j))
            if (Name.find("munu") != -1):
                CovMatrix.SetBinContent(i, j, IsoCovMatrix.GetBinContent(i, j) + TrigCovMatrix.GetBinContent(i, j) + RecoCovMatrix.GetBinContent(i, j) + RecoilCovMatrix.GetBinContent(i, j))

    for i in range(0, Hist1.GetNbinsX()):
        if (HUnfolded.GetBinContent(i+1) != 0):
            Hist1.SetBinContent(i+1, 100*sqrt(Covariance[i, i]) / HUnfolded.GetBinContent(i+1))
            Hist2.SetBinContent(i+1, 100*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1))
            Hist3.SetBinContent(i+1, 100*sqrt(CovBias.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1))
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
                     Variable+"_Xs_"+Name+".tex", "w+")
    latexFile.write("\\documentclass[12pt]{article} \n")
    latexFile.write("\\usepackage{amsmath}  \n")
    latexFile.write("\\usepackage{graphicx} \n")
    latexFile.write("\\usepackage{hyperref} \n")
    latexFile.write("\\usepackage{hyperref} \n")
    latexFile.write("\\usepackage[latin1]{inputenc} \n")
    latexFile.write("\\begin{document} \n")

    latexFile.write("\\begin{table}[] \n")
    latexFile.write("\\begin{tabular}{c|c|c|c|c|c|}\n")
    latexFile.write("\\cline{2-5}\n")

    latexFile.write(
        "& \\multicolumn{4}{c|}{%s}  \\\ \\hline \\hline \n" % (Indice))

    if ("pT" in Variable):
        i = 0
        max = HUnfolded.GetNbinsX()
        latexFile.write(
            "\\multicolumn{1}{|c|}{  range } & $d\sigma$/$dp^{T}_{l}$ [GeV]     & Stat Uncer.     & Unfolding bias     & Syst Uncer.    & Lumi Uncer.        \\\ \\hline \\hline \n")

    if ("Eta" in Variable):
        i = 0
        max = HUnfolded.GetNbinsX()
        latexFile.write("\\multicolumn{1}{|c|}{  range } & $d\sigma$/$d\eta_{l}$      & Stat Uncer.      & Unfolding bias     & Syst Uncer.   & Lumi Uncer.        \\\ \\hline \\hline \n")

    while i < max:
        if (HUnfolded.GetBinWidth(i+1) != 0 and Acceptance.GetBinContent(i+1) != 0):
            diffXs = HUnfolded.GetBinContent(i+1)/(Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Stat   = (100./diffXs)*sqrt(Covariance[i, i])                  / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Bias   = (100./diffXs)*sqrt(CovBias.GetBinContent(i+1, i+1))   / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Syst   = (100./diffXs)*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            if Energy == "5TeV":
                LumError = 0.016 * HUnfolded.GetBinContent(i+1)
            if Energy == "13TeV":
                LumError = 0.015 * HUnfolded.GetBinContent(i+1)
            Lumiunc   = (100./diffXs)*( LumError ) / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))

        if (HUnfolded.GetBinWidth(i+1) == 0 or Acceptance.GetBinContent(i+1) == 0):
            diffXs = 0
            Stat = 0
            Bias = 0
            Syst = 0
        lowedge = HUnfolded.GetBinLowEdge(i+1)
        highEdge = HUnfolded.GetBinLowEdge(i+2)
        latexFile.write("\\multicolumn{1}{|c|}{{[}%3.2f,  %3.2f{]}}  & %5.2f & %5.2f & %5.2f & %5.2f & %5.2f \\\ \\hline \n" % (lowedge, highEdge, diffXs, Stat, Bias, Syst, Lumiunc))
        i = i+1

    latexFile.write("\\end{tabular}\n")
    latexFile.write("\\end{table}\n")
    latexFile.write("\\end{document}\n")
    latexFile.close()

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

    Lumiunc = 0

    # Bias:
    CovBias = Bias.Get("CovMatrix_Iter_2")

    # define systematics:
    if (Name.find("enu") != -1):

        IDCovMatrix     = IdSF.Get("ElIDSys_Covariance_Iter2")
        TrigCovMatrix   = TrigSF.Get("TrigSys_Covariance_Iter2")
        RecoCovMatrix   = RecoSF.Get("RecoSys_Covariance_Iter2")
        IsoCovMatrix    = IsoSF.Get("IsoSys_Covariance_Iter2")
        RecoilCovMatrix = Recoil.Get("Recoil_Covariance_Iter2")
        CalibCovMatrix  = Calib.Get("Calib_Covariance_Iter2")

    if (Name.find("munu") != -1):

        TrigCovMatrix   = TrigSF.Get("MuTrigSys_Covariance_Iter2")
        RecoCovMatrix   = RecoSF.Get("MuRecoSys_Covariance_Iter2")
        IsoCovMatrix    = IsoSF.Get("MuIsoSys_Covariance_Iter2")
        RecoilCovMatrix = Recoil.Get("Recoil_Covariance_Iter2s")

    CovMatrix = RecoCovMatrix.Clone("CovMatrix")

    for i in range(1, 1 + CovMatrix.GetNbinsX()):
        for j in range(1, 1+CovMatrix.GetNbinsX()):
            if (Name.find("enu") != -1):
                CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i, j) + RecoCovMatrix.GetBinContent(i, j) + IsoCovMatrix.GetBinContent(i, j) + IDCovMatrix.GetBinContent(i, j) + CalibCovMatrix.GetBinContent(i, j) + RecoilCovMatrix.GetBinContent(i, j))
            if (Name.find("munu") != -1):
                CovMatrix.SetBinContent(i, j, IsoCovMatrix.GetBinContent(i, j) + TrigCovMatrix.GetBinContent(i, j) + RecoCovMatrix.GetBinContent(i, j) + RecoilCovMatrix.GetBinContent(i, j))

    for i in range(0, Hist1.GetNbinsX()):
        if (HUnfolded.GetBinContent(i+1) != 0):
            Hist1.SetBinContent(i+1, 100*sqrt(Covariance[i, i]) / HUnfolded.GetBinContent(i+1))
            Hist2.SetBinContent(i+1, 100*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1))
            Hist3.SetBinContent(i+1, 100*sqrt(CovBias.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1))
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
                     Variable+"_Xs_"+Name+".tex", "w+")
    latexFile.write("\\documentclass[12pt]{article} \n")
    latexFile.write("\\usepackage{amsmath}  \n")
    latexFile.write("\\usepackage{graphicx} \n")
    latexFile.write("\\usepackage{hyperref} \n")
    latexFile.write("\\usepackage{hyperref} \n")
    latexFile.write("\\usepackage[latin1]{inputenc} \n")
    latexFile.write("\\begin{document} \n")

    latexFile.write("\\begin{table}[] \n")
    latexFile.write("\\begin{tabular}{c|c|c|c|c|c|}\n")
    latexFile.write("\\cline{2-5}\n")

    latexFile.write(
        "& \\multicolumn{4}{c|}{%s}  \\\ \\hline \\hline \n" % (Indice))

    if ("pT" in Variable):
        i = 0
        max = HUnfolded.GetNbinsX()
        latexFile.write(
            "\\multicolumn{1}{|c|}{  range } & $d\sigma$/$dp^{T}_{l}$ [GeV]     & Stat Uncer.     & Unfolding bias     & Syst Uncer.    & Lumi Uncer.        \\\ \\hline \\hline \n")

    if ("Eta" in Variable):
        i = 0
        max = HUnfolded.GetNbinsX()
        latexFile.write("\\multicolumn{1}{|c|}{  range } & $d\sigma$/$d\eta_{l}$      & Stat Uncer.      & Unfolding bias     & Syst Uncer.   & Lumi Uncer.        \\\ \\hline \\hline \n")

    while i < max:
        if (HUnfolded.GetBinWidth(i+1) != 0 and Acceptance.GetBinContent(i+1) != 0):
            diffXs = HUnfolded.GetBinContent(i+1)/(Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Stat   = (100./diffXs)*sqrt(Covariance[i, i])                  / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Bias   = (100./diffXs)*sqrt(CovBias.GetBinContent(i+1, i+1))   / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Syst   = (100./diffXs)*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            if Energy == "5TeV":
                LumError = 0.016 * HUnfolded.GetBinContent(i+1)
            if Energy == "13TeV":
                LumError = 0.015 * HUnfolded.GetBinContent(i+1)
            Lumiunc   = (100./diffXs)*( LumError ) / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))

        if (HUnfolded.GetBinWidth(i+1) == 0 or Acceptance.GetBinContent(i+1) == 0):
            diffXs = 0
            Stat = 0
            Bias = 0
            Syst = 0
        lowedge = HUnfolded.GetBinLowEdge(i+1)
        highEdge = HUnfolded.GetBinLowEdge(i+2)
        latexFile.write("\\multicolumn{1}{|c|}{{[}%3.2f,  %3.2f{]}}  & %5.2f & %5.2f & %5.2f & %5.2f & %5.2f \\\ \\hline \n" % (lowedge, highEdge, diffXs, Stat, Bias, Syst, Lumiunc))
        i = i+1

    latexFile.write("\\end{tabular}\n")
    latexFile.write("\\end{table}\n")
    latexFile.write("\\end{document}\n")
    latexFile.close()

def GetDiffernetialXsMu(Summarize, Bias, TrigSF, RecoSF, IsoSF,IdSF, Calib, Recoil, Energy, Indice, Name, Lum, Variable):

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

    Lumiunc = 0

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
            Hist1.SetBinContent(i+1, 100*sqrt(Covariance[i, i]) / HUnfolded.GetBinContent(i+1))
            Hist2.SetBinContent(i+1, 100*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1))
            Hist3.SetBinContent(i+1, 100*sqrt(CovBias.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1))
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
    latexFile.write("\\begin{tabular}{c|c|c|c|c|c|}\n")
    latexFile.write("\\cline{2-5}\n")

    latexFile.write(
        "& \\multicolumn{4}{c|}{%s}  \\\ \\hline \\hline \n" % (Indice))

    if ("pT" in Variable):
        i = 0
        max = HUnfolded.GetNbinsX()
        latexFile.write("\\multicolumn{1}{|c|}{  range } & $d\sigma$/$dp^{T}_{l}$ [GeV]     & Stat uncertainty     & Unfolding bias     & Syst uncertainty   & Lumi uncertainty       \\\ \\hline \\hline \n")

    if ("Eta" in Variable):
        i = 0
        max = HUnfolded.GetNbinsX()
        latexFile.write("\\multicolumn{1}{|c|}{  range } & $d\sigma$/$d\eta_{l}$      & Stat uncertainty     & Unfolding bias     & Syst uncertainty  & Lumi uncertainty       \\\ \\hline \\hline \n")

    while i < max:
        if (HUnfolded.GetBinWidth(i+1) != 0 and Acceptance.GetBinContent(i+1) != 0):
            diffXs = HUnfolded.GetBinContent(i+1)/(Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Stat   = (100./diffXs)*sqrt(Covariance[i, i])                  / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Bias   = (100./diffXs)*sqrt(CovBias.GetBinContent(i+1, i+1))   / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            Syst   = (100./diffXs)*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))

            if Energy == "5TeV":
                LumError = 0.016 * HUnfolded.GetBinContent(i+1)
            if Energy == "13TeV":
                LumError = 0.015 * HUnfolded.GetBinContent(i+1)
            Lumiunc   = (100./diffXs)*( LumError ) / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))

        if (HUnfolded.GetBinWidth(i+1) == 0 or Acceptance.GetBinContent(i+1) == 0):
            diffXs = 0
            Stat = 0
            Bias = 0
            Syst = 0
        lowedge = HUnfolded.GetBinLowEdge(i+1)
        highEdge = HUnfolded.GetBinLowEdge(i+2)

        latexFile.write("\\multicolumn{1}{|c|}{{[}%3.2f,  %3.2f{]}}  & %5.2f & %5.2f & %5.2f & %5.2f & %5.2f \\\ \\hline \n" % (lowedge, highEdge, diffXs, Stat, Bias, Syst, Lumiunc))

        i = i+1

    latexFile.write("\\end{tabular}\n")
    latexFile.write("\\end{table}\n")
    latexFile.write("\\end{document}\n")
    latexFile.close()

def GetDiffernetialXsPlot(Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum, Variable, pT_Sherpa, Eta_Sherpa, pT_PhH7EG, Eta_PhH7EG):

    # define the variables:
    if "pT" in Variable:
        truth_sherpa = pT_Sherpa
        truth_PhH7EG = pT_PhH7EG
    if "Eta" in Variable:
        truth_sherpa = Eta_Sherpa
        truth_PhH7EG = Eta_PhH7EG
    
    print("********************* start *********************")

    data = Summarize.Get("data_hist")
    reco = Summarize.Get("reco_hist")
    truth = Summarize.Get("truth_hist")
    mig_mat = Summarize.Get("mig_hist")
    bg = Summarize.Get("Background_Total")

    TruthHisto = truth.Clone("TruthHisto")

    Hist1 = truth.Clone("Hist1")
    Hist2 = truth.Clone("Hist2")
    Hist3 = truth.Clone("Hist3")

    Acceptance = Summarize.Get("Acceptance_hist")

    # define Unfolded:
    HUnfolded = Summarize.Get("unfolded_data4")

    # Stat:
    Covariance = Summarize.Get("CovarianceMatrix_Data_Iter2")
    #print(Covariance, Covariance[26,26])

    Lumiunc = 0

    # Bias:
    CovBias = Bias.Get("CovMatrix_Iter_2")

    # define systematics:
    if (Name.find("enu") != -1):
        IDCovMatrix     = IdSF.Get("ElIDSys_Covariance_Iter2")
        TrigCovMatrix   = TrigSF.Get("TrigSys_Covariance_Iter2")
        RecoCovMatrix   = RecoSF.Get("RecoSys_Covariance_Iter2")
        IsoCovMatrix    = IsoSF.Get("IsoSys_Covariance_Iter2")
        RecoilCovMatrix = Recoil.Get("Recoil_Covariance_Iter2")
        CalibCovMatrix  = Calib.Get("Calib_Covariance_Iter2")

    if (Name.find("munu") != -1):
        CalibCovMatrix  = Calib.Get("Calib_Covariance_Iter2")
        IDCovMatrix     = IdSF.Get("MuTTVASys_Covariance_Iter2")
        TrigCovMatrix   = TrigSF.Get("TrigSys_Covariance_Iter2")
        RecoCovMatrix   = RecoSF.Get("RecoSys_Covariance_Iter2")
        IsoCovMatrix    = IsoSF.Get("IsoSys_Covariance_Iter2")
        RecoilCovMatrix = Recoil.Get("Recoil_Covariance_Iter2")

    CovMatrix = RecoCovMatrix.Clone("CovMatrix")

    print(RecoilCovMatrix)
    
    for i in range(1, 1 + CovMatrix.GetNbinsX()):
        for j in range(1, 1+CovMatrix.GetNbinsX()):
            CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i, j) + RecoCovMatrix.GetBinContent(i, j) + IsoCovMatrix.GetBinContent(i, j) + IDCovMatrix.GetBinContent(i, j) + CalibCovMatrix.GetBinContent(i, j) + RecoilCovMatrix.GetBinContent(i, j))

    for i in range(0, Hist1.GetNbinsX()):
        if (HUnfolded.GetBinContent(i+1) != 0):
            Hist1.SetBinContent(i+1, 100*sqrt(Covariance[i, i]) / HUnfolded.GetBinContent(i+1))
            Hist2.SetBinContent(i+1, 100*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1))
            Hist3.SetBinContent(i+1, 100*sqrt(CovBias.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1))
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

    listBinCenter = []
    listdiffXs = []
    listStat   = []
    listBias   = []
    listSyst   = []
    listLumi   = []

    # define the cross section:
    truth     = Summarize.Get("truth_hist")
    truth_Acc = (Summarize.Get("mig_hist")).ProjectionY()
    Xs       = ( HUnfolded.Integral()/Lum ) * ( truth.Integral() / truth_Acc.Integral() )

    max = HUnfolded.GetNbinsX()
    for i in range(0, max):
        if (HUnfolded.GetBinWidth(i+1) != 0 and Acceptance.GetBinContent(i+1) != 0):

            if Energy == "5TeV":
                LumError = 0.016 * HUnfolded.GetBinContent(i+1)
            if Energy == "13TeV":
                LumError = 0.015 * HUnfolded.GetBinContent(i+1)

            listBinCenter.append(  HUnfolded.GetBinCenter(i+1) )
            diffXs = (HUnfolded.GetBinContent(i+1)/(Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1)))
            listdiffXs.append(  (HUnfolded.GetBinContent(i+1)/(Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1)))                               )
            listStat.append(    ((1)*sqrt(Covariance[i, i])                  / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1)))    )
            listBias.append(    ((1)*sqrt(CovBias.GetBinContent(i+1, i+1))   / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1)))    )
            listSyst.append(    ((1)*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1)))    )
            listLumi.append(    ((1)*( LumError ) / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1)))                                )

        if (HUnfolded.GetBinWidth(i+1) == 0 or Acceptance.GetBinContent(i+1) == 0):
            listBinCenter.append(  HUnfolded.GetBinCenter(i+1) )
            listdiffXs.append( 0 )
            listStat.append(   0 )
            listBias.append(   0 )
            listSyst.append(   0 )
            listLumi.append(   0 )

    arrayBinCenter = np.array(listBinCenter)
    arraydiffXs    = np.array(listdiffXs)
    arrayStat      = np.array(listStat)
    arrayBias      = np.array(listBias)
    arraySyst      = np.array(listSyst)
    arrayLumi      = np.array(listLumi) 
    TotalStatSyst  = np.sqrt(arrayStat*arrayStat + arrayBias*arrayBias + arraySyst*arraySyst)
    TotalStatSystLum = np.sqrt(arrayStat*arrayStat + arrayBias*arrayBias + arraySyst*arraySyst + arrayLumi*arrayLumi)

    CrossSection = HUnfolded.Clone("HUnfolded")
    arrayCrossSection = []
    for i in range(0, CrossSection.GetNbinsX()):
        CrossSection.SetBinContent( i+1, arraydiffXs[i]      )
        CrossSection.SetBinError(   i+1, TotalStatSystLum[i] )
        arrayCrossSection.append( arraydiffXs[i] )

    for i in range(0, TruthHisto.GetNbinsX()):
            TruthHisto.SetBinContent(i+1,   (TruthHisto.GetBinContent(i+1)   / (Lum*TruthHisto.GetBinWidth(i+1)        ))  )
            TruthHisto.SetBinError(i+1,     (TruthHisto.GetBinError(i+1)     / (Lum*TruthHisto.GetBinWidth(i+1)        ))  )
      
            truth_sherpa.SetBinContent(i+1, (truth_sherpa.GetBinContent(i+1) / (Lum*truth_sherpa.GetBinWidth(i+1)      ))  )
            truth_sherpa.SetBinError(i+1,   (truth_sherpa.GetBinError(i+1)   / (Lum*truth_sherpa.GetBinWidth(i+1)      ))  )

            truth_PhH7EG.SetBinContent(i+1, (truth_PhH7EG.GetBinContent(i+1) / (Lum*truth_PhH7EG.GetBinWidth(i+1)      ))  )
            truth_PhH7EG.SetBinError(i+1,   (truth_PhH7EG.GetBinError(i+1)   / (Lum*truth_PhH7EG.GetBinWidth(i+1)      ))  )

    print("********************* verify *********************")    
    for i in range(0, CrossSection.GetNbinsX()):
        if Acceptance.GetBinContent(i+1) != 0:
            print(i+1,  HUnfolded.GetBinContent(i+1)/(Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1)), TruthHisto.GetBinContent(i+1), truth_PhH7EG.GetBinContent(i+1), truth_sherpa.GetBinContent(i+1))

    FunctionConparison(CrossSection, TruthHisto,  Indice, Name, Energy, Variable, truth_sherpa, truth_PhH7EG)
    



def GetDiffernetialXsPlotNormalised(Summarize, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum, Variable, pT_Sherpa, Eta_Sherpa, pT_PhH7EG, Eta_PhH7EG):

    # define the variables:
    if "pT" in Variable:
        truth_sherpa = pT_Sherpa
        truth_PhH7EG = pT_PhH7EG
    if "Eta" in Variable:
        truth_sherpa = Eta_Sherpa
        truth_PhH7EG = Eta_PhH7EG
    
    print("********************* start *********************")

    data = Summarize.Get("data_hist")
    reco = Summarize.Get("reco_hist")
    truth = Summarize.Get("truth_hist")
    mig_mat = Summarize.Get("mig_hist")
    bg = Summarize.Get("Background_Total")

    TruthHisto = truth.Clone("TruthHisto")

    Hist1 = truth.Clone("Hist1")
    Hist2 = truth.Clone("Hist2")
    Hist3 = truth.Clone("Hist3")

    Acceptance = Summarize.Get("Acceptance_hist")

    # define Unfolded:
    HUnfolded = Summarize.Get("unfolded_data4")


    # Cross section
    Xsfid = (HUnfolded.Integral()/Lum ) * ( Acceptance.Integral() )
    
    # Stat:
    Covariance = Summarize.Get("CovarianceMatrix_Data_Iter2")
    #print(Covariance, Covariance[26,26])

    Lumiunc = 0

    # Bias:
    CovBias = Bias.Get("CovMatrix_Iter_2")

    # define systematics:
    if (Name.find("enu") != -1):
        IDCovMatrix     = IdSF.Get("ElIDSys_Covariance_Iter2")
        TrigCovMatrix   = TrigSF.Get("TrigSys_Covariance_Iter2")
        RecoCovMatrix   = RecoSF.Get("RecoSys_Covariance_Iter2")
        IsoCovMatrix    = IsoSF.Get("IsoSys_Covariance_Iter2")
        RecoilCovMatrix = Recoil.Get("Recoil_Covariance_Iter2")
        CalibCovMatrix  = Calib.Get("Calib_Covariance_Iter2")

    if (Name.find("munu") != -1):
        CalibCovMatrix  = Calib.Get("Calib_Covariance_Iter2")
        IDCovMatrix     = IdSF.Get("MuTTVASys_Covariance_Iter2")
        TrigCovMatrix   = TrigSF.Get("TrigSys_Covariance_Iter2")
        RecoCovMatrix   = RecoSF.Get("RecoSys_Covariance_Iter2")
        IsoCovMatrix    = IsoSF.Get("IsoSys_Covariance_Iter2")
        RecoilCovMatrix = Recoil.Get("Recoil_Covariance_Iter2")

    CovMatrix = RecoCovMatrix.Clone("CovMatrix")

    print(RecoilCovMatrix)
    
    for i in range(1, 1 + CovMatrix.GetNbinsX()):
        for j in range(1, 1+CovMatrix.GetNbinsX()):
            CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i, j) + RecoCovMatrix.GetBinContent(i, j) + IsoCovMatrix.GetBinContent(i, j) + IDCovMatrix.GetBinContent(i, j) + CalibCovMatrix.GetBinContent(i, j) + RecoilCovMatrix.GetBinContent(i, j))

    for i in range(0, Hist1.GetNbinsX()):
        if (HUnfolded.GetBinContent(i+1) != 0):
            Hist1.SetBinContent(i+1, (1/Xsfid)*100*sqrt(Covariance[i, i]) / HUnfolded.GetBinContent(i+1))
            Hist2.SetBinContent(i+1, (1/Xsfid)*100*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1))
            Hist3.SetBinContent(i+1, (1/Xsfid)*100*sqrt(CovBias.GetBinContent(i+1, i+1)) / HUnfolded.GetBinContent(i+1))
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

    listBinCenter = []
    listdiffXs = []
    listStat   = []
    listBias   = []
    listSyst   = []
    listLumi   = []

    # define the cross section:
    truth     = Summarize.Get("truth_hist")
    truth_Acc = (Summarize.Get("mig_hist")).ProjectionY()
    Xs       = ( HUnfolded.Integral()/Lum ) * ( truth.Integral() / truth_Acc.Integral() )

    max = HUnfolded.GetNbinsX()
    for i in range(0, max):
        if (HUnfolded.GetBinWidth(i+1) != 0 and Acceptance.GetBinContent(i+1) != 0):

            if Energy == "5TeV":
                LumError = 0.016 * HUnfolded.GetBinContent(i+1)
            if Energy == "13TeV":
                LumError = 0.015 * HUnfolded.GetBinContent(i+1)

            listBinCenter.append(  HUnfolded.GetBinCenter(i+1) )
            diffXs = (1/Xsfid)*(HUnfolded.GetBinContent(i+1)/(Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1)))
            listdiffXs.append(  (1/Xsfid)*(HUnfolded.GetBinContent(i+1)/(Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1)))                               )
            listStat.append(    (1/Xsfid)*((1.)*sqrt(Covariance[i, i])                  / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1)))    )
            listBias.append(    (1/Xsfid)*((1.)*sqrt(CovBias.GetBinContent(i+1, i+1))   / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1)))    )
            listSyst.append(    (1/Xsfid)*((1.)*sqrt(CovMatrix.GetBinContent(i+1, i+1)) / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1)))    )
            listLumi.append(    (1/Xsfid)*((1.)*( LumError ) / (Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1)))                                )

        if (HUnfolded.GetBinWidth(i+1) == 0 or Acceptance.GetBinContent(i+1) == 0):
            listBinCenter.append(  HUnfolded.GetBinCenter(i+1) )
            listdiffXs.append( 0 )
            listStat.append(   0 )
            listBias.append(   0 )
            listSyst.append(   0 )
            listLumi.append(   0 )

    arrayBinCenter = np.array(listBinCenter)
    arraydiffXs    = np.array(listdiffXs)
    arrayStat      = np.array(listStat)
    arrayBias      = np.array(listBias)
    arraySyst      = np.array(listSyst)
    arrayLumi      = np.array(listLumi) 
    TotalStatSyst  = np.sqrt(arrayStat*arrayStat + arrayBias*arrayBias + arraySyst*arraySyst)
    TotalStatSystLum = np.sqrt(arrayStat*arrayStat + arrayBias*arrayBias + arraySyst*arraySyst + arrayLumi*arrayLumi)

    CrossSection = HUnfolded.Clone("HUnfolded")
    arrayCrossSection = []
    for i in range(0, CrossSection.GetNbinsX()):
        CrossSection.SetBinContent( i+1, arraydiffXs[i]      )
        CrossSection.SetBinError(   i+1, TotalStatSystLum[i] )
        arrayCrossSection.append( arraydiffXs[i] )

    for i in range(0, TruthHisto.GetNbinsX()):
            TruthHisto.SetBinContent(i+1,   (1/Xsfid)*(TruthHisto.GetBinContent(i+1)   / (Lum*TruthHisto.GetBinWidth(i+1)        ))  )
            TruthHisto.SetBinError(i+1,     (1/Xsfid)*(TruthHisto.GetBinError(i+1)     / (Lum*TruthHisto.GetBinWidth(i+1)        ))  )
      
            truth_sherpa.SetBinContent(i+1, (1/Xsfid)*(truth_sherpa.GetBinContent(i+1) / (Lum*truth_sherpa.GetBinWidth(i+1)      ))  )
            truth_sherpa.SetBinError(i+1,   (1/Xsfid)*(truth_sherpa.GetBinError(i+1)   / (Lum*truth_sherpa.GetBinWidth(i+1)      ))  )

            truth_PhH7EG.SetBinContent(i+1, (1/Xsfid)*(truth_PhH7EG.GetBinContent(i+1) / (Lum*truth_PhH7EG.GetBinWidth(i+1)      ))  )
            truth_PhH7EG.SetBinError(i+1,   (1/Xsfid)*(truth_PhH7EG.GetBinError(i+1)   / (Lum*truth_PhH7EG.GetBinWidth(i+1)      ))  )


    FunctionConparisonNormalised(CrossSection, TruthHisto,  Indice, Name, Energy, Variable, truth_sherpa, truth_PhH7EG)
    
