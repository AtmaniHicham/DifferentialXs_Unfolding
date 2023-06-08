#!/usr/bin/env python
# -*-coding:Latin-1 -*
import numpy as np

import atlasplots
from   atlasplots import atlas_style as astyle
from   atlasplots import utils
from   atlasplots import config_reader as config

from math import *

import ROOT
import ROOT as root
from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TH2D, TMatrixD, TMatrixT, TMatrix
ROOT.gSystem.Load("libMatrix.so")

def makeLegend(hists, xmin, ymin, xmax, ymax):
    legend = root.TLegend(xmin, ymin, xmax, ymax)
    legend.SetTextSize(0.03)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetBorderSize(0)
    for hist in hists:
        legend.AddEntry(hist, hist.GetName())
    return legend

class RebenningClass:
    """Classe représentant une personne"""

    def __init__(self, SummarizLargeBins, SummarizSmallBins):
        """Initialisation """

    def RebenningOfSmallBins(self, SummarizLargeBins, SummarizSmallBins, NumberOfIterations):
        print(NumberOfIterations)
        # Get the covariance matrix and unfolded data with 'Large and small' bins
        CovStatSmall   = SummarizSmallBins.Get("Covariance_Matrix"+str(NumberOfIterations))
        CovStatLarge   = SummarizLargeBins.Get("Covariance_Matrix"+str(NumberOfIterations))
        HFinner        = SummarizSmallBins.Get("Unfolded_data_iteration_"+str(NumberOfIterations));
        HLarger        = SummarizLargeBins.Get("Unfolded_data_iteration_"+str(NumberOfIterations));
        # define the jacobian and covariance matrix
        JacobienMatrix = CovStatSmall.Clone("JacobienMatrix")
        NewCovarMatrix = CovStatLarge.Clone("NewCovarMatrix")

        #Define the binning
        LargeBins = 15
        SmallBins = 160
        Jaco            = np.full((LargeBins, SmallBins), 0)
        JacobianT       = np.full((SmallBins,LargeBins),0)
        CovarianceData  = np.full((SmallBins,SmallBins),0)
        ProduitA        = np.full((LargeBins,SmallBins),0)
        ProduitB        = np.full((LargeBins,LargeBins),0)
        i=0
        j=0
        for i in range(0,SmallBins):
            for j in range(0,SmallBins):
                CovarianceData[i,j] = CovStatSmall.GetBinContent(i,j)
        for k in range(0,LargeBins):
            for i in range(0,SmallBins):
                if(HLarger.GetXaxis().GetBinLowEdge(k) < HFinner.GetXaxis().GetBinUpEdge(i) and  HFinner.GetXaxis().GetBinUpEdge(i) <= HLarger.GetXaxis().GetBinUpEdge(k)):
                    Jaco[k,i]=1
        JacobianT = Jaco.transpose()
        ProduitA  = Jaco.dot(CovarianceData)
        ProduitB  = ProduitA.dot(JacobianT)

        for k in range(0,LargeBins):
            for i in range(0,LargeBins):
                NewCovarMatrix.SetBinContent(k, i, ProduitB[k,i])
        return NewCovarMatrix


    def ComparisonOfrebinnedresults(self, SummarizLargeBins, SummarizSmallBins, RebinnedSmallBinsMatrix, NumberOfIterations):
        CovStatSmall   = SummarizSmallBins.Get("Covariance_Matrix"+str(NumberOfIterations))
        CovStatLarge   = SummarizLargeBins.Get("Covariance_Matrix"+str(NumberOfIterations))

        Unfolded_data_Large = SummarizLargeBins.Get("Unfolded_data_iteration_"+str(NumberOfIterations))
        Unfolded_data_Small = SummarizSmallBins.Get("Unfolded_data_iteration_"+str(NumberOfIterations))

        HStatLarge          = SummarizLargeBins.Get("Unfolded_MC")
        HStatSmall          = SummarizSmallBins.Get("Unfolded_MC")
        HStatSmallrebinned  = HStatLarge.Clone("HStatSmallrebinned")

        HStatLarge.SetName("Large Bins")
        HStatSmall.SetName("Small Bins")
        HStatSmallrebinned.SetName("rebinned Small Bins")

        LargeBins = 15
        SmallBins = 160
        for i in range(0,SmallBins):
            if( Unfolded_data_Small.GetBinContent(i) != 0):
                HStatSmall.SetBinContent(i,100*sqrt(CovStatSmall.GetBinContent(i,i))/Unfolded_data_Small.GetBinContent(i) )
            HStatSmall.SetBinError(i,0)
        for i in range(0,LargeBins):
            if( Unfolded_data_Large.GetBinContent(i) != 0):
                    HStatLarge.SetBinContent(i,100*sqrt(CovStatLarge.GetBinContent(i,i))/Unfolded_data_Large.GetBinContent(i) )
                    HStatSmallrebinned.SetBinContent(i,100*sqrt(RebinnedSmallBinsMatrix.GetBinContent(i,i))/Unfolded_data_Large.GetBinContent(i) )
            HStatLarge.SetBinError(i,0)
            HStatSmallrebinned.SetBinError(i,0)

        hists = []
        hists.append(HStatSmall)
        hists.append(HStatLarge)
        hists.append(HStatSmallrebinned)
        Color = 0
        for hist in hists:
            Color +=1
            hist.SetMarkerColor(Color)
            hist.SetStats(0)
            hist.SetTitle("")
            hist.GetXaxis().SetRangeUser(0,60)
            hist.GetXaxis().SetTitle("p^{W}_{T} [GeV]")
            hist.GetYaxis().SetTitle("Stat [%]")
            utils.SetHistogramLine(hist, Color, 2, Color, 1) # def SetHistogramLine(hist, color=1, width=1, style=1, alpha=1):
        # Draw
        c1 = root.TCanvas("c1", "The FillRandom example", 0, 0, 800, 600)
        for hist in hists:
            hist.Draw("same")
        astyle.ATLASLabel(0.2, 0.84, "Internal")
        utils.DrawText(0.2, 0.78, "W^{-}#rightarrow e^{-}#nu, 5TeV, Iteration "+ str(NumberOfIterations),1)
        legend = makeLegend(hists,0.2, 0.45,0.4,0.6)
        legend.Draw("same")
        c1.Update()
        c1.Print("LargeAndSmallBins.pdf")

        #ProduitA=Jaco*CovarianceData
