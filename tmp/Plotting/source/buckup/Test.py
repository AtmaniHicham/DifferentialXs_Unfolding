#!/usr/bin/env python
# -*-coding:Latin-1 -*

import atlasplots
from   atlasplots import atlas_style as astyle
from   atlasplots import utils
from   atlasplots import config_reader as config

from math import *

import ROOT
import ROOT as root
from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TAttFill

def makeLegend(hists, xmin, ymin, xmax, ymax):
    legend = root.TLegend(xmin, ymin, xmax, ymax)
    legend.SetTextSize(0.045)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetBorderSize(0)
    for hist in hists:
        legend.AddEntry(hist, hist.GetName())
    return legend



    
def DiffernetialXsPlot(Summarize_minusenu5, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum):

        data       = Summarize_minusenu5.Get("data_hist")
        reco       = Summarize_minusenu5.Get("reco_hist")
        truth      = Summarize_minusenu5.Get("truth_hist")
        mig_mat    = Summarize_minusenu5.Get("mig_hist")
        bg         = Summarize_minusenu5.Get("Background_Total")

        # define Unfolded:
        HUnfolded  = Summarize_minusenu5.Get("unfolded_data2")

        # Stat:        
        Covariance = Summarize_minusenu5.Get("CovarianceMatrix_Iter2")

        # Bias:
        CovBias      = Bias.Get("CovMatrix_Iter_2")

        # define systematics:
        if  (Name.find("enu")  != -1):

            IDCovMatrix      = IdSF.Get(    "ElIDSys_Covariance_Iter2"    )
            TrigCovMatrix    = TrigSF.Get(  "elTrigSys_Covariance_Iter2"  )
            RecoCovMatrix    = RecoSF.Get(  "elRecoSys_Covariance_Iter2"  )
            IsoCovMatrix     = IsoSF.Get(   "elIsoSys_Covariance_Iter2"   )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter2"     )
            CalibCovMatrix   = Calib.Get(   "Calib_Covariance_Iter2"      )

        if  (Name.find("munu")  != -1):

            TrigCovMatrix    = TrigSF.Get(  "MuTrigSys_Covariance_Iter2"  )
            RecoCovMatrix    = RecoSF.Get(  "MuRecoSys_Covariance_Iter2"  )
            IsoCovMatrix     = IsoSF.Get(   "MuIsoSys_Covariance_Iter2"   )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter2"     )

        CovMatrix   = RecoCovMatrix.Clone("CovMatrix")

        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                if  (Name.find("enu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + 
						  CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) )
                if  (Name.find("munu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )

        # add STAT, Bias, Syst:
        Uncer = truth.Clone("Uncer")

        for i in range(1, 1 + Uncer.GetNbinsX()):
            if  (Name.find("enu")  != -1):
                Uncer.SetBinContent(i,  TrigCovMatrix.GetBinContent(i,i) + RecoCovMatrix.GetBinContent(i,i) + IsoCovMatrix.GetBinContent(i,i) + IDCovMatrix.GetBinContent(i,i) + 
					CalibCovMatrix.GetBinContent(i,i) + RecoilCovMatrix.GetBinContent(i,i) + Covariance.GetBinContent(i,i) + CovBias.GetBinContent(i,i) )
            if  (Name.find("munu")  != -1):
                Uncer.SetBinContent(i,  TrigCovMatrix.GetBinContent(i,i) + RecoCovMatrix.GetBinContent(i,i) + IsoCovMatrix.GetBinContent(i,i) + RecoilCovMatrix.GetBinContent(i,i) + 
					Covariance.GetBinContent(i,i) + CovBias.GetBinContent(i,i) )

        for i in range(1, 1 + Uncer.GetNbinsX()):
            if(HUnfolded.GetBinContent(i) != 0):
                Uncer.SetBinContent(i, 1 )
                Uncer.SetBinError(i, 100*sqrt(Uncer.GetBinContent(i)) / HUnfolded.GetBinContent(i))

        # define the acceptance correction:
        truth_noMiss = (mig_mat.ProjectionY()).Clone("truth_noMiss")
        Acctepatance = truth.Clone("Acctepatance")
        Acctepatance.Divide(truth_noMiss)

        DXs = truth.Clone("DXs")
        DMC = truth.Clone("DXs")

        for i in range(0, HUnfolded.GetNbinsX() ):
            if(HUnfolded.GetBinContent(i) != 0):
                DXs.SetBinContent(i+1, Acctepatance.GetBinContent(i+1)*HUnfolded.GetBinContent(i+1) / (Lum*HUnfolded.GetBinWidth(i+1)) )
                DMC.SetBinContent(i+1, truth.GetBinContent(i+1) / (Lum*truth.GetBinWidth(i+1)) )
                DMC.SetBinError(i+1, 0)

        # define the ratio:
        Ratio1 = DXs.Clone("Ratio1")
        Ratio1.Divide(DMC)

	        
	for i in range(0, Ratio1.GetNbinsX() ):
	    Ratio1.SetBinError(i+1, 0)


	for i in range(0, Uncer.GetNbinsX()):
	    print(i+1, Uncer.GetBinContent(i+1), Uncer.GetBinError(i+1))













        astyle.SetAtlasStyle()
        c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )

        #define the pad
        pad1 = root.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.Draw()
        pad1.SetLogy()
        pad1.SetFillStyle(4000)
        pad1.cd()

        DXs.SetStats(0)
        DXs.SetTitle("")
        DXs.SetLineColor(4)
        DXs.SetLineWidth(2)
        DXs.SetLineWidth(2)
        DXs.SetLineColor(2)
        DXs.GetXaxis().SetRangeUser(-2., 2.)
        DXs.GetXaxis().SetLabelSize(0.16)
        DXs.GetXaxis().SetTitleSize(0.15)
        #DXs.GetYaxis().SetLabelSize(0.2)
        #DXs.GetYaxis().SetTitleSize(0.07)
        #DXs.GetYaxis().SetTitle("d#sigma / d#eta_{l}")
        #DXs.GetYaxis().SetLabelOffset(-0.4)
        DXs.GetXaxis().SetTitle("#eta_{l}")
        DXs.GetXaxis().SetTitleOffset(1.5)
        DMC.SetLineWidth(2)
        DMC.SetMarkerColor(4)
        DMC.SetMarkerStyle(4)
        DMC.SetMarkerSize(1)
        DMC.SetLineColor(4)

	print(DXs.GetYaxis().GetLabelSize())
	
        DXs.Draw()
        DMC.Draw("same")
        astyle.ATLASLabel( 0.68, 0.84, "Internal")
        utils.DrawText(    0.68, 0.78,  Indice)
        DXs.SetName(" Unfolded data")
        DMC.SetName(" PowhegPythia8")
        hists = []
        hists.append(DXs)
        hists.append(DMC)
        legendN = makeLegend(hists,0.64, 0.5,0.92,0.75)
        legendN.Draw("same")

	'''
        c1.Update()
        c1.cd()
        pad2 = root.TPad("pad2", "pad2", 0, 0.06, 1, 0.32)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.4)
        pad2.Draw()
        pad2.SetFillStyle(4000)
        pad2.cd()
       
        Uncer.GetYaxis().SetNdivisions(5)
        Uncer.GetYaxis().SetTitle("Pred./Data")
        Uncer.GetXaxis().SetTitle("#eta_{l}")
        Uncer.GetXaxis().SetLabelSize(0.16)
        Uncer.GetXaxis().SetTitleSize(0.14)
        Uncer.GetXaxis().SetRangeUser(-2., 2.)
        Uncer.GetYaxis().SetLabelSize(0.16)
        Uncer.GetYaxis().SetTitleSize(0.14)
        Uncer.GetYaxis().SetTitleOffset(0.4)
        Uncer.GetXaxis().SetTitleOffset(1.5)
        Uncer.GetYaxis().SetRangeUser(0.92,1.08)
        Uncer.SetLineWidth(3)
        Uncer.SetLineColor(4)
	Uncer.SetFillStyle(4050)
	Uncer.SetFillColor(2)
	Uncer.Draw("same E1E3")
	Ratio1.SetLineWidth(2)
        Ratio1.Draw("same")
	
        line1 = ROOT.TLine(-2,0.95,2,0.95)
        line2 = ROOT.TLine(-2,1.05,2,1.05)
        line3 = ROOT.TLine(-2,1.00,2,1.00)
        line1.SetLineStyle(2)
        line2.SetLineStyle(2)
        line3.SetLineStyle(2)
        line1.Draw("same")
        line2.Draw("same")
        line3.Draw("same")
	'''
        c1.Print("Output/CrossSection/Differential_Xs_"+Name+".pdf")



