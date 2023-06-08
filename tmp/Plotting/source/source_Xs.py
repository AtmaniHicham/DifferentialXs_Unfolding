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
