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

from ROOT import TCanvas, TGraph
from ROOT import gROOT
from math import sin

from array import array
from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TAttFill, TMatrixD, THilbertMatrixD, TDecompSVD, TGraphErrors, TH2D, TText

def makeLegend(hists, xmin, ymin, xmax, ymax):
    legend = root.TLegend(xmin, ymin, xmax, ymax)
    legend.SetTextSize(0.03)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetBorderSize(0)
    for hist in hists:
        legend.AddEntry(hist, hist.GetName())
    return legend


def TrasfertoTH2F(MCovarianceMatrix):
        Nbin = MCovarianceMatrix.GetNrows()
        Migration = TH2D("Migration", "Migration", Nbin, 0, Nbin, Nbin, 0, Nbin)
        for i in range(0, Nbin):
            for j in range(0, Nbin):
                    Migration.SetBinContent(i+1, j+1, MCovarianceMatrix[i][j])
        return Migration



class DiffCrossSection:
    """Classe représentant une personne"""

    def __init__(self):
        """Constructeur de notre classe"""


    def GetDiffernetialXsComp(self, Summarize_plusenu5, Summarize_plusmunu5, Summarize_minusenu5, Summarize_minusmunu5, list1, Energy, Indice, Name, Lum):

        # **************************************************************************************************************************************************
	# ******************************************************************* Wplusenu *********************************************************************
        # **************************************************************************************************************************************************

        data1       = Summarize_plusenu5.Get("data_hist")
        reco1       = Summarize_plusenu5.Get("reco_hist")
        truth1      = Summarize_plusenu5.Get("truth_hist")
        mig_mat1    = Summarize_plusenu5.Get("mig_hist")
        bg1         = Summarize_plusenu5.Get("Background_Total")

        # define Unfolded:
        HUnfolded1  = Summarize_plusenu5.Get("unfolded_data4")

        # define the acceptance correction:
        truth_noMiss1 = (mig_mat1.ProjectionY()).Clone("truth_noMiss")
        Acctepatance1 = truth1.Clone("Acctepatance")
        Acctepatance1.Divide(truth_noMiss1)

        DXs1 = truth1.Clone("DXs1")

	# Syst
        HStatError = HUnfolded1.Clone("HStatError")
        HSystTotal = HUnfolded1.Clone("HSystTotal")
        CovarianceMat = Summarize_plusenu5.Get('CovarianceMatrix_Iter3')
        Covariance    = TrasfertoTH2F(CovarianceMat)

	IDCovMatrix      = list1[0].Get( "ElIDSys_Covariance_Iter3"  )
        IsoCovMatrix     = list1[1].Get( "ElIsoSys_Covariance_Iter3" )
	RecoCovMatrix    = list1[2].Get( "ElRecoSys_Covariance_Iter3"  )
        TrigCovMatrix    = list1[3].Get( "ElTrigSys_Covariance_Iter3"  )
	RecoilCovMatrix  = list1[4].Get( "Recoil_Covariance_Iter3"   )
	CalibCovMatrix   = list1[5].Get( "Calib_Covariance_Iter3"    )
        CovBias      	 = list1[6].Get( "CovMatrix_Iter_3"	     )

        HistSyst 	 = truth1.Clone("Uncer")

        for i in range(1, 1 + RecoCovMatrix.GetNbinsX()):
                HistSyst.SetBinContent(i,  TrigCovMatrix.GetBinContent(i,i) + RecoCovMatrix.GetBinContent(i,i) + IsoCovMatrix.GetBinContent(i,i) + IDCovMatrix.GetBinContent(i,i) + CalibCovMatrix.GetBinContent(i,i) + RecoilCovMatrix.GetBinContent(i,i) + Covariance.GetBinContent(i,i) + CovBias.GetBinContent(i,i) )

	
        # **************************************************************************************************************************************************
        # ******************************************************************* Wplusmunu *********************************************************************
        # **************************************************************************************************************************************************

        data2       = Summarize_plusmunu5.Get("data_hist")
        reco2       = Summarize_plusmunu5.Get("reco_hist")
        truth2      = Summarize_plusmunu5.Get("truth_hist")
        mig_mat2    = Summarize_plusmunu5.Get("mig_hist")
        bg2         = Summarize_plusmunu5.Get("Background_Total")

        # define Unfolded:
        HUnfolded2  = Summarize_plusmunu5.Get("unfolded_data4")

        # define the acceptance correction:
        truth_noMiss2 = (mig_mat2.ProjectionY()).Clone("truth_noMiss")
        Acctepatance2 = truth2.Clone("Acctepatance")
        Acctepatance2.Divide(truth_noMiss2)

        DXs2 = truth2.Clone("DXs2")


        # **************************************************************************************************************************************************
        # ******************************************************************* Wminusenu *********************************************************************
        # **************************************************************************************************************************************************

        data3       = Summarize_minusenu5.Get("data_hist")
        reco3       = Summarize_minusenu5.Get("reco_hist")
        truth3      = Summarize_minusenu5.Get("truth_hist")
        mig_mat3    = Summarize_minusenu5.Get("mig_hist")
        bg3         = Summarize_minusenu5.Get("Background_Total")

        # define Unfolded:
        HUnfolded3  = Summarize_minusenu5.Get("unfolded_data4")

        # define the acceptance correction:
        truth_noMiss3 = (mig_mat3.ProjectionY()).Clone("truth_noMiss")
        Acctepatance3 = truth3.Clone("Acctepatance")
        Acctepatance3.Divide(truth_noMiss3)

        DXs3  = truth3.Clone("DXs2")

        # **************************************************************************************************************************************************
        # ******************************************************************* Wminusmunu *********************************************************************
        # **************************************************************************************************************************************************
                
        data4       = Summarize_minusmunu5.Get("data_hist")
        reco4       = Summarize_minusmunu5.Get("reco_hist")
        truth4      = Summarize_minusmunu5.Get("truth_hist")
        mig_mat4    = Summarize_minusmunu5.Get("mig_hist")
        bg4         = Summarize_minusmunu5.Get("Background_Total")

        # define Unfolded:
        HUnfolded4  = Summarize_minusmunu5.Get("unfolded_data4")
        
        # define the acceptance correction:
        truth_noMiss4 = (mig_mat4.ProjectionY()).Clone("truth_noMiss")
        Acctepatance4 = truth4.Clone("Acctepatance")
        Acctepatance4.Divide(truth_noMiss4)

        DXs4 = truth4.Clone("DXs2")


        # **************************************************************************************************************************************************
        # ******************************************************************* Cal the Xs *******************************************************************
        # **************************************************************************************************************************************************

        for i in range(0, HUnfolded1.GetNbinsX() ):
                DXs1.SetBinContent(i+1, Acctepatance1.GetBinContent(i+1)*HUnfolded1.GetBinContent(i+1) / (Lum*HUnfolded1.GetBinWidth(i+1)) )
                DXs2.SetBinContent(i+1, Acctepatance2.GetBinContent(i+1)*HUnfolded2.GetBinContent(i+1) / (Lum*HUnfolded2.GetBinWidth(i+1)) )
                DXs3.SetBinContent(i+1, Acctepatance3.GetBinContent(i+1)*HUnfolded3.GetBinContent(i+1) / (Lum*HUnfolded3.GetBinWidth(i+1)) )
                DXs4.SetBinContent(i+1, Acctepatance4.GetBinContent(i+1)*HUnfolded4.GetBinContent(i+1) / (Lum*HUnfolded4.GetBinWidth(i+1)) )
		DXs1.SetBinError(i+1,   Acctepatance1.GetBinContent(i+1)*sqrt(HistSyst.GetBinContent(i+1))   / (Lum*HUnfolded4.GetBinWidth(i+1)) )
                DXs2.SetBinError(i+1,   Acctepatance1.GetBinContent(i+1)*sqrt(HistSyst.GetBinContent(i+1))   / (Lum*HUnfolded4.GetBinWidth(i+1)) )
                DXs3.SetBinError(i+1,   Acctepatance1.GetBinContent(i+1)*sqrt(HistSyst.GetBinContent(i+1))   / (Lum*HUnfolded4.GetBinWidth(i+1)) )
                DXs4.SetBinError(i+1,   Acctepatance1.GetBinContent(i+1)*sqrt(HistSyst.GetBinContent(i+1))   / (Lum*HUnfolded4.GetBinWidth(i+1)) )

        #for i in range(0, HUnfolded1.GetNbinsX() ):
	#	DXs1.SetBinError(i+1, 100*DXs1.GetBinError(i+1)/ DXs1.GetBinContent(i+1))
        #        DXs2.SetBinError(i+1, 100*DXs2.GetBinError(i+1)/ DXs2.GetBinContent(i+1))
        #        DXs3.SetBinError(i+1, 100*DXs3.GetBinError(i+1)/ DXs3.GetBinContent(i+1))
        #        DXs4.SetBinError(i+1, 100*DXs4.GetBinError(i+1)/ DXs4.GetBinContent(i+1))

        for i in range(0, DXs1.GetNbinsX() ):
	     print("Bin ", i, DXs1.GetBinContent(i+1), HistSyst.GetBinContent(i+1))




        astyle.SetAtlasStyle()

        c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )
        c1.Draw()
        #c1.SetLogy()
        c1.SetFillStyle(4000)
        c1.cd()

        DXs1.GetYaxis().SetRangeUser(150, 600)
        DXs1.GetXaxis().SetRangeUser(-2.4,2.4)
        DXs1.SetStats(0)
        DXs1.SetTitle("")
        DXs1.GetYaxis().SetTitle("d#sigma / d#eta_{l} ")
        DXs1.GetXaxis().SetTitle("#eta_{l}")
        DXs1.GetYaxis().SetTitleOffset(0.8)
        DXs1.GetXaxis().SetTitleOffset(1.4)
        DXs1.GetYaxis().SetTitleSize(0.07)
        DXs1.GetYaxis().SetLabelSize(0.04)
        DXs1.GetYaxis().SetLabelOffset(0.01)

        DXs1.SetLineWidth(2)
        DXs2.SetLineWidth(2)

        DXs1.SetMarkerColor(2)
        DXs2.SetMarkerColor(4)

        #DXs1.SetMarkerStyle(4)
        #DXs2.SetMarkerStyle(4)

        DXs1.SetMarkerSize(4)
        DXs2.SetMarkerSize(4)

        DXs1.SetLineColor(2)
        DXs2.SetLineColor(4)

        DXs3.SetLineWidth(2)
        DXs3.SetMarkerColor(3)
        #DXs3.SetMarkerStyle(4)
        DXs3.SetMarkerSize(4)
        DXs3.SetLineColor(3)

        DXs4.SetLineWidth(2)
        DXs4.SetMarkerColor(1)
        #DXs4.SetMarkerStyle(4)
        DXs4.SetMarkerSize(4)
        DXs4.SetLineColor(1)

        DXs1.Draw("H")
        DXs2.Draw("same H")
        DXs3.Draw("same H")
        DXs4.Draw("same H")


        astyle.ATLASLabel( 0.5, 0.86, "Internal, 5TeV")

        DXs1.SetName(" W^{+} \\rightarrow e \\nu")
        DXs2.SetName(" W^{+} \\rightarrow \\mu \\nu")
        DXs3.SetName(" W^{-} \\rightarrow e \\nu")
        DXs4.SetName(" W^{-} \\rightarrow \\mu \\nu")

        hists = []
        hists.append(DXs1)
        hists.append(DXs2)

        hists1 = []
        hists1.append(DXs3)
        hists1.append(DXs4)

        legendN = makeLegend(hists,0.6, 0.6,0.82,0.5)
        legendNp = makeLegend(hists1,0.2, 0.6,0.42,0.5)
        legendN.Draw("same")
        legendNp.Draw("same")


        c1.Print("Output/CrossSection/Differential_Xs_"+Name+".pdf")

	'''	
        # define the ratio:
        Ratio1 = DXs1.Clone("Ratio1")
        Ratio1.Divide(DXs2)

        c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )
        astyle.SetAtlasStyle()

        #define the pad
        pad1 = root.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.Draw()
        pad1.SetLogy()
        pad1.SetFillStyle(4000)
        pad1.cd()

        DXs1.GetXaxis().SetRangeUser(-2.4,2.4)
        DXs1.SetStats(0)
        DXs1.SetTitle("")
        DXs1.GetYaxis().SetTitle("d#sigma / d#eta_{l} ")
        DXs1.GetXaxis().SetTitle("#eta_{l}")
        DXs1.GetYaxis().SetTitleOffset(0.9)
        DXs1.GetXaxis().SetTitleOffset(1.4)
        DXs1.GetYaxis().SetTitleSize(0.07)
        DXs1.GetYaxis().SetLabelSize(0.08)

        DXs1.SetLineWidth(2)
        DXs2.SetLineWidth(2)

        DXs1.SetMarkerColor(2)
        DXs2.SetMarkerColor(4)

        #DXs1.SetMarkerStyle(4)
        #DXs2.SetMarkerStyle(4)

        #DXs1.SetMarkerSize(1)
        #DXs2.SetMarkerSize(1)

        DXs1.SetLineColor(2)
        DXs2.SetLineColor(4)

	DXs2.SetFillColor(2)
        DXs1.Draw()
        DXs2.Draw("same")
        #DXs3.Draw("same")
        #DXs4.Draw("same")

        astyle.ATLASLabel( 0.6, 0.86, "Internal")
        utils.DrawText(    0.6, 0.81,  Indice)

        DXs1.SetName(" W plus enu")
        DXs2.SetName(" W plus munu")
        hists = []
        hists.append(DXs1)
        hists.append(DXs2)

        legendN = makeLegend(hists,0.6, 0.6,0.92,0.8)
        legendN.Draw("same")

        c1.Update()
        c1.cd()
        astyle.SetAtlasStyle()

        pad2 = root.TPad("pad2", "pad2", 0, 0.06, 1, 0.32)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.4)
        pad2.Draw()
        pad2.SetFillStyle(4000)
        pad2.cd()

        Ratio1.GetYaxis().SetRangeUser(0.91,1.09)
        Ratio1.GetYaxis().SetNdivisions(5)
        Ratio1.GetYaxis().SetTitle("Pred./Data")
        Ratio1.GetXaxis().SetTitle("#eta_{l}")
        Ratio1.GetXaxis().SetLabelSize(0.15)
        Ratio1.GetXaxis().SetTitleSize(0.16)
        Ratio1.GetXaxis().SetTitleOffset(1.1)

        Ratio1.GetYaxis().SetTitleOffset(0.5)
        Ratio1.GetYaxis().SetTitleSize(0.11)
        Ratio1.GetYaxis().SetLabelSize(0.15)
        Ratio1.GetYaxis().CenterTitle()

        Ratio1.GetXaxis().SetRangeUser(-1.8,1.8)
        Ratio1.SetLineWidth(2)
        Ratio1.SetLineColor(2)
        Ratio1.Draw("same")

        line1 = ROOT.TLine(-2,0.95,2,0.95)
        line2 = ROOT.TLine(-2,1.05,2,1.05)
        line1.SetLineStyle(2)
        line2.SetLineStyle(2)

        line1.Draw("same")
        line2.Draw("same")
	'''
