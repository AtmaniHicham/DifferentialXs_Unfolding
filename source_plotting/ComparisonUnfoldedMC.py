#!/usr/bin/env python
# -*-coding:Latin-1 -*

import atlasplots
from   atlasplots import atlas_style as astyle
from   atlasplots import utils
from   atlasplots import config_reader as config

from math import *
import matplotlib.pyplot as plt

import ROOT
import ROOT as root
from ROOT import gStyle
import numpy as np
from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TAttFill, TMatrixD, THilbertMatrixD, TDecompSVD

def makeLegend(hists, xmin, ymin, xmax, ymax):
    legend = root.TLegend(xmin, ymin, xmax, ymax)
    legend.SetTextSize(0.03)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetBorderSize(0)
    for hist in hists:
        legend.AddEntry(hist, hist.GetName())
    return legend


class ComparisonUnfoldedMC:
    """Classe représentant une personne"""

    def __init__(self):
        """Constructeur de notre classe"""

    def ComparisonNominal(self, Summarize_minusenu5, MCsamples_minusenu5, MCsamples_minusenu5_DYturbo, MCsamples_minusenu5_Pythia8, MCsamples_minusenu5_Sherpa, Niter, Lum, Energy, Indice, Name):

        # define the binning
        Binning = [0, 5, 10, 15, 21, 29, 40, 52, 64, 77, 92, 115, 145, 175, 220, 600]

        # Define the Unfolded distribution
        UnfoldedDistribution = Summarize_minusenu5.Get("Unfolded_data_iteration_4")

        # read the Total Uncertainties
        CovarianceStatError  = Summarize_minusenu5.Get("Covariance_Matrix1")
        StatError            = UnfoldedDistribution.Clone("Stat Uncertainties")

        i     = 0
        while i < StatError.GetNbinsX():
            StatError.SetBinContent(i+1,  1)
            StatError.SetBinError(i+1,    sqrt(CovarianceStatError.GetBinContent(i+1, i+1))/UnfoldedDistribution.GetBinContent(i+1) )
            i = i+1

        i     = 0
        while i < StatError.GetNbinsX():
            StatError.SetBinContent(i+1,  1)
            #StatError.SetBinError(i+1,    sqrt(StatError.GetBinError(i+1)*StatError.GetBinError(i+1) + 0.05*0.05) )
            i = i+1

        # Correct the unfolded distribution
        Truth_Raw            = Summarize_minusenu5.Get("Truth_MC_Raw")
        Truth_Corrected      = Summarize_minusenu5.Get("htruth_noMiss")
        Truth_Corrected.Divide(Truth_Raw)

        i     = 0
        while i < UnfoldedDistribution.GetNbinsX():
            UnfoldedDistribution.SetBinContent(i+1, UnfoldedDistribution.GetBinContent(i+1)*(1/Truth_Corrected.GetBinContent(i+1)) )
            UnfoldedDistribution.SetBinError(i+1,   0 )
            i=i+1

        # define the Initial Truth distribution
        TruthPowhegPythia8   = MCsamples_minusenu5.Get("TruthSelection/WpT_Truth_5TeV_cut4")
        TruthDYturbo         = MCsamples_minusenu5_DYturbo.Get("TruthSelection/WpT_Truth_5TeV_cut4")
        TruthPythia8         = MCsamples_minusenu5_Pythia8.Get("TruthSelection/WpT_Truth_5TeV_cut4")
        TruthSherpa8         = MCsamples_minusenu5_Sherpa.Get("TruthSelection/WpT_Truth_5TeV_cut4")

        i     = 0
        while i < TruthPowhegPythia8.GetNbinsX():
            TruthPowhegPythia8.SetBinError(i+1,   0 )
            TruthDYturbo.SetBinError(i+1,   0 )
            TruthPythia8.SetBinError(i+1,   0 )
            TruthSherpa8.SetBinError(i+1,   0 )

            i=i+1

        # Define the ratio for each histos
        TruthPowhegPythia8.Divide(UnfoldedDistribution)
        TruthDYturbo.Divide(UnfoldedDistribution)
        TruthPythia8.Divide(UnfoldedDistribution)
        TruthSherpa8.Divide(UnfoldedDistribution)


        # Define the style and color for each hists
        TruthPowhegPythia8.SetLineWidth(2)
        TruthPowhegPythia8.SetLineColor(1)
        TruthPowhegPythia8.SetName("PowhegPythia8")

        TruthDYturbo.SetLineWidth(2)
        TruthDYturbo.SetLineColor(1)
        TruthDYturbo.SetLineStyle(2)
        TruthDYturbo.SetName("DYRES")

        TruthPythia8.SetLineWidth(2)
        TruthPythia8.SetLineColor(2)
        TruthPythia8.SetName("Pythia8")

        TruthSherpa8.SetLineWidth(2)
        TruthSherpa8.SetLineColor(2)
        TruthSherpa8.SetLineStyle(2)
        TruthSherpa8.SetName("Sherpa8")

        UnfoldedDistribution.SetLineWidth(2)
        UnfoldedDistribution.SetLineColor(4)
        UnfoldedDistribution.SetLineStyle(2)
        UnfoldedDistribution.SetName("Unfolded distribution")

        StatError.SetFillColor(4)
        StatError.SetFillStyle(3002)


        # define the legend
        histsN = []
        histsN.append(TruthPowhegPythia8)
        histsN.append(TruthDYturbo)
        histsN.append(TruthPythia8)
        histsN.append(TruthSherpa8)
        histsN.append(StatError)
        legendN = makeLegend(histsN,0.6, 0.67,0.84,0.86)

        # define the Canvas output
        c1N = root.TCanvas("c1N", "The FillRandom example", 0, 0, 800, 600)
        c1N.Draw()
        #c1N.SetLogy()
        c1N.cd()
        TruthPowhegPythia8.GetXaxis().SetRangeUser(0,100)
        TruthPowhegPythia8.GetYaxis().SetRangeUser(0.9,1.2)
        TruthPowhegPythia8.SetStats(0)
        TruthPowhegPythia8.SetTitle("")
        TruthPowhegPythia8.GetYaxis().SetTitle(" Theory / Data")
        TruthPowhegPythia8.GetXaxis().SetTitle("p_{w}^{T} [GeV]")
        TruthPowhegPythia8.Draw("same")
        TruthDYturbo.Draw("same")
        TruthPythia8.Draw("same")
        TruthSherpa8.Draw("same")
        StatError.Draw("same E1 E2")
        legendN.Draw("same")
        astyle.ATLASLabel(0.15, 0.8, "Internal")
        utils.DrawText(0.15, 0.73, Indice)
        c1N.Print("Output/CrossSection/normalized_Differential_CrossSection_"+Name+".pdf")
