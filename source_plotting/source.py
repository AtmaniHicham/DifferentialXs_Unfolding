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

from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TH1, TF1
from math import *

def test():
    print("sdean")


def GetAcceptanceFactors(inputFile, channel, Variable, Indice):
    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, ax = aplt.subplots(1, 1, name="fig1", figsize=(800, 600))

    # Get the histogram:
    hist   = inputFile.Get("Acceptance_hist")
    hist1   = inputFile.Get("Efficiency_hist")

    # Draw the histogram on these axes
    ax.plot(hist,  label="Acceptance",     linecolor=root.kBlue+1,    labelfmt="L")
    ax.plot(hist1, label="Efficiency",     linecolor=root.kRed+1,     labelfmt="L")
    ax.set_ylim(0.4, 1.3)

    # Add extra space at top of plot to make room for labels
    ax.add_margins(top=0.16)

    # Set axis titles
    ax.set_ylabel("Correction factors")
    if (Variable == "Lepton_pT"):
        ax.set_xlabel("p_{T}^{lepton}")
        ax.set_xlim(25, 100)
    if (Variable == "Lepton_Eta"):
        ax.set_xlabel("\eta^{lepton}")

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax.text(0.2, 0.86, Indice, size=22, align=13)

    # Add legend
    ax.legend(loc=(0.65, 0.8, 0.95, 0.92))

    # Save the plot as a PDF
    fig.savefig("Output/output_"+channel+"/"+Variable+"/"+channel+"_Acceptance.pdf")

# ----------------------------------------------------------------------------------------------------------------------------------

def ShowNominalDistribution(inputFile, channel, Variable, Indice):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, (ax1, ax2) = aplt.ratio_plot(name="fig1", figsize=(800, 800), hspace=0.05)

    # define the histograms:
    unfolded_Data1   = inputFile.Get("unfolded_data1")
    unfolded_Data2   = inputFile.Get("unfolded_data4")
    unfolded_MC      = inputFile.Get("unfolded_MC1")
    Truth            = inputFile.Get("truth_prjection")

    # Draw the histograms on these axes
    ax1.plot(Truth,             linecolor=root.kBlack,      label="Truth",          labelfmt="L")
    ax1.plot(unfolded_MC,       linecolor=root.kBlue+1,     label="unfolded MC",    labelfmt="L")
    ax1.plot(unfolded_Data1,    linecolor=root.kRed+1,      label="unfolded Data - Iter 1",  labelfmt="L")
    #ax1.plot(unfolded_Data2,    linecolor=8,                label="unfolded Data - Iter 4",  labelfmt="L")

    # legend position
    ax1.legend(loc=(0.55, 0.7, 0.8, 0.9))

    # Calculate and draw the ratio
    ratio_hist = unfolded_MC.Clone("ratio_hist")
    ratio_hist.Divide(Truth)
    ax2.plot(ratio_hist, linecolor=root.kBlue+1)

    ratio_hist1 = unfolded_Data1.Clone("ratio_hist1")
    ratio_hist1.Divide(Truth)
    ratio_hist2 = unfolded_Data2.Clone("ratio_hist2")
    ratio_hist2.Divide(Truth)

    ax2.plot(ratio_hist1, linecolor=root.kRed+1)
    #ax2.plot(ratio_hist2, linecolor=8)

    # Add extra space at top of plot to make room for labels
    ax1.add_margins(top=0.16)

    # Set axis titles
    if (Variable == "Lepton_pT"):
        ax2.set_xlabel("p_{T}^{lepton}")
        ax1.set_xlim(25, 100)
        ax2.set_xlim(25, 100)
        line = root.TLine(25, 1, 100, 1)
        ax2.plot(line)
    if (Variable == "Lepton_Eta"):
        ax2.set_xlabel("\eta^{lepton}")
        line = root.TLine(ax1.get_xlim()[0], 1, ax1.get_xlim()[1], 1)
        ax2.plot(line)

    ax1.set_ylabel("Events")
    ax2.set_ylabel("Unfolded / Truth", loc="centre")
    ax2.set_ylim(0.88, 1.12)


    # Add extra space at top and bottom of ratio panel
    ax2.add_margins(top=0.1, bottom=0.1)

    # Go back to top axes to add labels
    ax1.cd()

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax1.text(0.2, 0.84, Indice, size=22, align=13)

    # Save the plot as a PDF
    fig.savefig("Output/output_"+channel+"/"+Variable+"/"+channel+"_NominalPlots.pdf")

# ----------------------------------------------------------------------------------------------------------------------------------

def CompareBias(Summarize, inputFile, NumberOfIterationMstopinimal, NumberOfIterationMaximal, channel, Variable, Indice):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, ax     = aplt.subplots(1, 1, name="fig1", figsize=(800, 600))

    bias1       = Summarize.Get("unfolded_MC1")

    for k in range(NumberOfIterationMstopinimal, NumberOfIterationMaximal):
        hist        = inputFile.Get('Bias_Iter_' + str(k))
        for i in range(0, hist.GetNbinsX()):
            bias1.SetBinContent(i+1, hist.GetBinContent(i+1))
            bias1.SetBinError(i+1, 0)
        ax.plot(bias1.Clone(),  label="Iteration "+str(k),     linecolor=k,    labelfmt="L")
    
    #for i in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
    #    hist = inputFile.Get('Bias_Iter_' + str(i))
    #    ax.plot(hist,  label="Acceptance",     linecolor=root.kBlue+1,    labelfmt="L")

    # Add extra space at top of plot to make room for labels
    ax.add_margins(top=0.16)

    # Set axis titles
    ax.set_ylabel("Bias Uncertainty (%)")
    if (Variable == "Lepton_pT"):
        ax.set_xlabel("p_{T}^{lepton}")
        ax.set_xlim(25, 100)
    if (Variable == "Lepton_Eta"):
        ax.set_xlabel("\eta^{lepton}")

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax.text(0.2, 0.86, Indice, size=22, align=13)
    ax.set_ylim(-0.5, 1)

    # Add legend
    ax.legend(loc=(0.65, 0.6, 0.95, 0.92))

    # Save the plot as a PDF
    fig.savefig("Output/output_"+channel+"/"+Variable+"/"+channel+"_Bias.pdf")

# ----------------------------------------------------------------------------------------------------------------------------------

def CompareStatError(Summarize, NumberOfIterationMstopinimal, NumberOfIterationMaximal, channel, Variable, Indice):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, ax     = aplt.subplots(1, 1, name="fig1", figsize=(800, 600))

    bias1       = Summarize.Get("unfolded_MC1")
    Covariance  = Summarize.Get("CovarianceMatrix_Data_Iter1")
    
    for k in range(NumberOfIterationMstopinimal, NumberOfIterationMaximal):
        hist        = Summarize.Get('unfolded_data' + str(k))
        for i in range(0, hist.GetNbinsX()):
            if(hist.GetBinError(i+1) != 0):
                bias1.SetBinContent(i+1, 100.*sqrt(Covariance[i, i])/hist.GetBinContent(i+1))
                print(k, i+1, 100.*sqrt(Covariance[i, i])/hist.GetBinContent(i+1))
                bias1.SetBinError(i+1, 0)
            if(hist.GetBinError(i+1) == 0):
                bias1.SetBinContent(i+1, 0)
                bias1.SetBinError(i+1, 0)
        ax.plot(bias1.Clone(),  label="Iteration "+str(k),     linecolor=k,    labelfmt="L")
    
    #for i in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
    #    hist = inputFile.Get('Bias_Iter_' + str(i))
    #    ax.plot(hist,  label="Acceptance",     linecolor=root.kBlue+1,    labelfmt="L")

    # Add extra space at top of plot to make room for labels
    ax.add_margins(top=0.16)

    # Set axis titles
    ax.set_ylabel("Stat Uncertainty")
    if (Variable == "Lepton_pT"):
        ax.set_xlim(25, 100)
        ax.set_xlabel("p_{T}^{lepton}")
    if (Variable == "Lepton_Eta"):
        ax.set_xlabel("\eta^{lepton}")

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax.text(0.2, 0.86, Indice, size=22, align=13)

    # Add legend
    ax.legend(loc=(0.45, 0.6, 0.8, 0.92))

    
    # Save the plot as a PDF
    fig.savefig("Output/output_"+channel+"/"+Variable+"/"+channel+"_Stat.pdf")

# ----------------------------------------------------------------------------------------------------------------------------------

def BiasProcedure(Summarize, inputFile, channel, Variable, Indice):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, ax     = aplt.subplots(1, 1, name="fig1", figsize=(800, 600))

    reco    = inputFile.Get("reco_hist")
    reco_W  = inputFile.Get("reco_hist_Weighted")
    data    = Summarize.Get("dataCorrected")

    reco.Divide(data)
    reco_W.Divide(data)

    ax.plot(reco.Clone(),    label="Data/MC(reco) ",           linecolor=2,    labelfmt="L")
    ax.plot(reco_W.Clone(),  label="Data/MC(reco weighted) ",  linecolor=4,    labelfmt="L")

    # Add extra space at top of plot to make room for labels
    ax.add_margins(top=0.16)
    ax.set_ylim(0.8, 1.4)

    # Set axis titles
    ax.set_ylabel("Ratio Data/MC")
    if (Variable == "Lepton_pT"):
        ax.set_xlim(25, 100)
        ax.set_xlabel("p_{T}^{lepton}")
    if (Variable == "Lepton_Eta"):
        ax.set_xlabel("\eta^{lepton}")

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax.text(0.2, 0.86, Indice, size=22, align=13)

    # Add legend
    ax.legend(loc=(0.45, 0.6, 0.8, 0.92))

    
    # Save the plot as a PDF
    fig.savefig("Output/output_"+channel+"/"+Variable+"/"+channel+"_BiasProce.pdf")

# ----------------------------------------------------------------------------------------------------------------------------------

def BiasProcedureOld(Summarize, inputFile, channel, Variable, Indice):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, ax     = aplt.subplots(1, 1, name="fig1", figsize=(800, 600))

    reco    = inputFile.Get("reco_hist")
    reco_W  = inputFile.Get("reco_hist_Weighted")
    data    = Summarize.Get("dataCorrected")

    print(reco_W)
    for i in range(0, reco.GetNbinsX()):
        print(  reco.GetBinContent(i+1)  )

    reco.Divide(data)
    #eco_W.Divide(data)

    # define the bias for Unfolding
    if(sys.argv[4] == "Lepton_Eta"):
        VarMin = -2.4
        VarMax =  2.4
    if(sys.argv[4] == "Lepton_pT"):
        VarMin =  0
        VarMax =  100
    
    f1 = TF1("f1", "pol6", VarMin, VarMax)
    reco.Fit("f1")

    print(reco)
    print(reco_W)

    # plot 
    #ax.plot(reco,       label="Data/reconstructed level",               linecolor=1,    labelfmt="L")
    #ax.plot(reco_W,     label="Data/reconstructed level weighted",      linecolor=2,    labelfmt="L")
    ax.plot(reco_W.Clone(),  label="Data/reconstructed level ",     linecolor=1,    labelfmt="L")

    """
    # Add extra space at top of plot to make room for labels
    ax.add_margins(top=0.16)

    # Set axis titles
    ax.set_ylabel("Stat Uncertainty")
    if (Variable == "Lepton_pT"):
        ax.set_xlim(25, 100)
        ax.set_xlabel("p_{T}^{lepton}")
    if (Variable == "Lepton_Eta"):
        ax.set_xlabel("\eta{lepton}")

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax.text(0.2, 0.86, Indice, size=22, align=13)
    ax.set_xlim(25, 60)
    #ax.set_ylim(-0.3, 0.5)

    # Add legend
    ax.legend(loc=(0.65, 0.6, 0.95, 0.92))

    # Save the plot as a PDF
    fig.savefig("Output/output_"+channel+"/"+Variable+"/"+channel+"_Bias_Procedure.pdf")
    """
# ----------------------------------------------------------------------------------------------------------------------------------

def CompareSystId(inputFile, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Variable, Indice):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, ax     = aplt.subplots(1, 1, name="fig1", figsize=(800, 600))

    for k in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
        hist    = inputFile.Get('ElIDSys_Systematics_Iter' + str(k))
        ax.plot(hist.Clone(),  label="Iteration "+str(k),     linecolor=k,    labelfmt="L")
    
    #for i in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
    #    hist = inputFile.Get('Bias_Iter_' + str(i))
    #    ax.plot(hist,  label="Acceptance",     linecolor=root.kBlue+1,    labelfmt="L")

    # Add extra space at top of plot to make room for labels
    ax.add_margins(top=0.16)

    # Set axis titles
    ax.set_ylabel("ID SF Uncertainty")
    if (Variable == "Lepton_pT"):
        ax.set_xlim(25, 100)
        ax.set_xlabel("p_{T}^{lepton}")
    if (Variable == "Lepton_Eta"):
        ax.set_xlabel("\eta^{lepton}")

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax.text(0.2, 0.86, Indice, size=22, align=13)
    #ax.set_ylim(-0.3, 0.5)

    # Add legend
    ax.legend(loc=(0.45, 0.6, 0.8, 0.92))

    
    # Save the plot as a PDF
    fig.savefig("Output/output_"+channel+"/"+Variable+"/"+channel+"_ID.pdf")

# ----------------------------------------------------------------------------------------------------------------------------------

def CompareSystIso(inputFile, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Variable, Indice):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, ax     = aplt.subplots(1, 1, name="fig1", figsize=(800, 600))

    for k in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
        hist    = inputFile.Get('IsoSys_Systematics_Iter' + str(k))
        ax.plot(hist.Clone(),  label="Iteration "+str(k),     linecolor=k,    labelfmt="L")
    
    #for i in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
    #    hist = inputFile.Get('Bias_Iter_' + str(i))
    #    ax.plot(hist,  label="Acceptance",     linecolor=root.kBlue+1,    labelfmt="L")

    # Add extra space at top of plot to make room for labels
    ax.add_margins(top=0.16)

    # Set axis titles
    ax.set_ylabel("Iso SF Uncertainty")
    if (Variable == "Lepton_pT"):
        ax.set_xlim(25, 100)
        ax.set_xlabel("p_{T}^{lepton}")
    if (Variable == "Lepton_Eta"):
        ax.set_xlabel("\eta^{lepton}")

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax.text(0.2, 0.86, Indice, size=22, align=13)
    #ax.set_ylim(-0.3, 0.5)

    # Add legend
    ax.legend(loc=(0.45, 0.6, 0.8, 0.92))

    
    # Save the plot as a PDF
    fig.savefig("Output/output_"+channel+"/"+Variable+"/"+channel+"_Iso.pdf")

# ----------------------------------------------------------------------------------------------------------------------------------

def CompareSystReco(inputFile, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Variable, Indice):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, ax     = aplt.subplots(1, 1, name="fig1", figsize=(800, 600))

    for k in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
        hist    = inputFile.Get('RecoSys_Systematics_Iter' + str(k))
        ax.plot(hist.Clone(),  label="Iteration "+str(k),     linecolor=k,    labelfmt="L")
    
    #for i in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
    #    hist = inputFile.Get('Bias_Iter_' + str(i))
    #    ax.plot(hist,  label="Acceptance",     linecolor=root.kBlue+1,    labelfmt="L")

    # Add extra space at top of plot to make room for labels
    ax.add_margins(top=0.16)

    # Set axis titles
    ax.set_ylabel("Reco SF Uncertainty")
    if (Variable == "Lepton_pT"):
        ax.set_xlim(25, 100)
        ax.set_xlabel("p_{T}^{lepton}")
    if (Variable == "Lepton_Eta"):
        ax.set_xlabel("\eta^{lepton}")

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax.text(0.2, 0.86, Indice, size=22, align=13)
    #ax.set_ylim(-0.3, 0.5)

    # Add legend
    ax.legend(loc=(0.45, 0.6, 0.8, 0.92))

    
    # Save the plot as a PDF
    fig.savefig("Output/output_"+channel+"/"+Variable+"/"+channel+"_Reco.pdf")

# ----------------------------------------------------------------------------------------------------------------------------------

def CompareSystTrig(inputFile, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Variable, Indice):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, ax     = aplt.subplots(1, 1, name="fig1", figsize=(800, 600))

    for k in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
        hist    = inputFile.Get('TrigSys_Systematics_Iter' + str(k))
        ax.plot(hist.Clone(),  label="Iteration "+str(k),     linecolor=k,    labelfmt="L")
    
    #for i in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
    #    hist = inputFile.Get('Bias_Iter_' + str(i))
    #    ax.plot(hist,  label="Acceptance",     linecolor=root.kBlue+1,    labelfmt="L")

    # Add extra space at top of plot to make room for labels
    ax.add_margins(top=0.16)

    # Set axis titles
    ax.set_ylabel("Trig SF Uncertainty")
    if (Variable == "Lepton_pT"):
        ax.set_xlim(25, 100)
        ax.set_xlabel("p_{T}^{lepton}")
    if (Variable == "Lepton_Eta"):
        ax.set_xlabel("\eta^{lepton}")

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax.text(0.2, 0.86, Indice, size=22, align=13)
    #ax.set_ylim(-0.3, 0.5)

    # Add legend
    ax.legend(loc=(0.45, 0.6, 0.8, 0.92))

    # Save the plot as a PDF
    fig.savefig("Output/output_"+channel+"/"+Variable+"/"+channel+"_Trig.pdf")

# ----------------------------------------------------------------------------------------------------------------------------------

def CompareSystRecoil(inputFile, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Variable, Indice):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, ax     = aplt.subplots(1, 1, name="fig1", figsize=(800, 600))

    for k in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
        hist    = inputFile.Get('Recoil_Systematics_Iter' + str(k))
        ax.plot(hist.Clone(),  label="Iteration "+str(k),     linecolor=k,    labelfmt="L")
    
    #for i in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
    #    hist = inputFile.Get('Bias_Iter_' + str(i))
    #    ax.plot(hist,  label="Acceptance",     linecolor=root.kBlue+1,    labelfmt="L")

    # Add extra space at top of plot to make room for labels
    ax.add_margins(top=0.16)

    # Set axis titles
    ax.set_ylabel("Recoil Calibration Uncertainty [%]")
    if (Variable == "Lepton_pT"):
        ax.set_xlim(25, 100)
        ax.set_xlabel("p_{T}^{lepton}")
    if (Variable == "Lepton_Eta"):
        ax.set_xlabel("\eta^{lepton}")

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax.text(0.2, 0.86, Indice, size=22, align=13)
    #ax.set_ylim(-0.3, 0.5)

    # Add legend
    ax.legend(loc=(0.45, 0.6, 0.8, 0.92))

    
    # Save the plot as a PDF
    fig.savefig("Output/output_"+channel+"/"+Variable+"/"+channel+"_Recoil.pdf")

# ----------------------------------------------------------------------------------------------------------------------------------

def CompareSystCalib(inputFile, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Variable, Indice):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, ax     = aplt.subplots(1, 1, name="fig1", figsize=(800, 600))

    for k in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
        hist    = inputFile.Get('Calib_Systematics_Iter' + str(k))
        ax.plot(hist.Clone(),  label="Iteration "+str(k),     linecolor=k,    labelfmt="L")
    
    #for i in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
    #    hist = inputFile.Get('Bias_Iter_' + str(i))
    #    ax.plot(hist,  label="Acceptance",     linecolor=root.kBlue+1,    labelfmt="L")

    # Add extra space at top of plot to make room for labels
    ax.add_margins(top=0.16)

    # Set axis titles
    ax.set_ylabel("Lepton Calibration Uncertainty [%]")
    if (Variable == "Lepton_pT"):
        ax.set_xlim(25, 100)
        ax.set_xlabel("p_{T}^{lepton}")
    if (Variable == "Lepton_Eta"):
        ax.set_xlabel("\eta^{lepton}")

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax.text(0.2, 0.86, Indice, size=22, align=13)
    #ax.set_ylim(-0.3, 0.5)

    # Add legend
    ax.legend(loc=(0.45, 0.6, 0.8, 0.92))

    
    # Save the plot as a PDF
    fig.savefig("Output/output_"+channel+"/"+Variable+"/"+channel+"_Calib.pdf")

# ----------------------------------------------------------------------------------------------------------------------------------

def CompareSystEl(inputFile, Bias, IdSF, IsoSF, RecoSF, TrigSF, Recoil, Calib, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Variable, Indice):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, ax         = aplt.subplots(1, 1, name="fig1", figsize=(800, 600))

    TotalUncert     = inputFile.Get("unfolded_MC1")
    
    for k in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
        hist            = inputFile.Get('unfolded_data' + str(k))
        Covariance_stat = inputFile.Get("CovarianceMatrix_Data_Iter"+ str(k))
        Covariance_bias = Bias.Get("CovMatrix_Iter_"+ str(k))

        Covariance_Iso  = IsoSF.Get('IsoSys_Covariance_Iter' + str(k))
        Covariance_Reco = RecoSF.Get('RecoSys_Covariance_Iter' + str(k))
        Covariance_Trig = TrigSF.Get('TrigSys_Covariance_Iter' + str(k))
        Covariance_ID    = IdSF.Get('ElIDSys_Covariance_Iter' + str(k))

        Covariance_Calib = Calib.Get('Calib_Covariance_Iter' + str(k))
        Covariance_Recoil = Recoil.Get('Recoil_Covariance_Iter' + str(k))

        for i in range(0, hist.GetNbinsX()):
            if(hist.GetBinError(i+1) != 0):
                sum = Covariance_stat[i,i] + Covariance_bias.GetBinContent(i+1, i+1)+Covariance_Iso.GetBinContent(i+1, i+1)+Covariance_Reco.GetBinContent(i+1, i+1)+Covariance_Trig.GetBinContent(i+1, i+1)+Covariance_Recoil.GetBinContent(i+1, i+1)+Covariance_ID.GetBinContent(i+1, i+1)+Covariance_Calib.GetBinContent(i+1, i+1)
                TotalUncert.SetBinContent(i+1, 100.*sqrt(sum)/hist.GetBinContent(i+1))
                TotalUncert.SetBinError(i+1, 0)
            if(hist.GetBinError(i+1) == 0):
                TotalUncert.SetBinContent(i+1, 0)
                TotalUncert.SetBinError(i+1, 0)
        ax.plot(TotalUncert.Clone(),  label="Iteration "+str(k),     linecolor=k,    labelfmt="L")
    
    # Add extra space at top of plot to make room for labels
    ax.add_margins(top=0.16)

    # Set axis titles
    ax.set_ylim(0, 6)
    ax.set_ylabel("Total Uncertainty [%]")
    if (Variable == "Lepton_pT"):
        ax.set_xlim(25, 100)
        ax.set_xlabel("p_{T}^{lepton}")
    if (Variable == "Lepton_Eta"):
        ax.set_xlabel("\eta^{lepton}")

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax.text(0.2, 0.86, Indice, size=22, align=13)
    #ax.set_ylim(-0.3, 0.5)

    # Add legend
    ax.legend(loc=(0.45, 0.6, 0.8, 0.92))

    
    # Save the plot as a PDF
    fig.savefig("Output/output_"+channel+"/"+Variable+"/"+channel+"_TotalSystematics.pdf")


# ----------------------------------------------------------------------------------------------------------------------------------

def CompareSystMu(inputFile, Bias, IsoSF, RecoSF, TrigSF, Recoil, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Variable, Indice):

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, ax         = aplt.subplots(1, 1, name="fig1", figsize=(800, 600))

    TotalUncert     = inputFile.Get("unfolded_MC1")
    
    for k in range(NumberOfIterationMinimal, NumberOfIterationMaximal):
        hist            = inputFile.Get('unfolded_data' + str(k))
        Covariance_stat = inputFile.Get("CovarianceMatrix_Data_Iter"+ str(k))
        Covariance_bias = Bias.Get("CovMatrix_Iter_"+ str(k))

        Covariance_Iso  = IsoSF.Get('IsoSys_Covariance_Iter' + str(k))
        Covariance_Reco = RecoSF.Get('RecoSys_Covariance_Iter' + str(k))
        Covariance_Trig = TrigSF.Get('TrigSys_Covariance_Iter' + str(k))

        Covariance_Recoil = Recoil.Get('Recoil_Covariance_Iter' + str(k))

        for i in range(0, hist.GetNbinsX()):
            if(hist.GetBinError(i+1) != 0):
                sum = Covariance_stat[i,i] + Covariance_bias.GetBinContent(i+1, i+1)+Covariance_Iso.GetBinContent(i+1, i+1)+Covariance_Reco.GetBinContent(i+1, i+1)+Covariance_Trig.GetBinContent(i+1, i+1)+Covariance_Recoil.GetBinContent(i+1, i+1)
                TotalUncert.SetBinContent(i+1, 100.*sqrt(sum)/hist.GetBinContent(i+1))
                TotalUncert.SetBinError(i+1, 0)
            if(hist.GetBinError(i+1) == 0):
                TotalUncert.SetBinContent(i+1, 0)
                TotalUncert.SetBinError(i+1, 0)
        ax.plot(TotalUncert.Clone(),  label="Iteration "+str(k),     linecolor=k,    labelfmt="L")
    
    # Add extra space at top of plot to make room for labels
    ax.add_margins(top=0.16)

    # Set axis titles
    ax.set_ylim(0, 6)
    ax.set_ylabel("Total Uncertainty [%]")
    if (Variable == "Lepton_pT"):
        ax.set_xlim(25, 100)
        ax.set_xlabel("p_{T}^{lepton}")
    if (Variable == "Lepton_Eta"):
        ax.set_xlabel("\eta^{lepton}")

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax.text(0.2, 0.86, Indice, size=22, align=13)
    #ax.set_ylim(-0.3, 0.5)

    # Add legend
    ax.legend(loc=(0.45, 0.6, 0.8, 0.92))

    
    # Save the plot as a PDF
    fig.savefig("Output/output_"+channel+"/"+Variable+"/"+channel+"_TotalSystematics.pdf")




