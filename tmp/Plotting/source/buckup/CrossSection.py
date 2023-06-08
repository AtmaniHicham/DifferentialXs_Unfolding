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
from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TAttFill, TMatrixD, THilbertMatrixD, TDecompSVD, TGraphErrors, TText

def makeLegend(hists, xmin, ymin, xmax, ymax):
    legend = root.TLegend(xmin, ymin, xmax, ymax)
    legend.SetTextSize(0.03)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetBorderSize(0)
    for hist in hists:
        legend.AddEntry(hist, hist.GetName())
    return legend


class CrossSection:
    """Classe représentant une personne"""

    def __init__(self):
        """Constructeur de notre classe"""

    def GetFiducialXs(self, Summarize_minusenu5, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum):
              
        HBias      = Bias.Get("Bias_Iter_1")
        HUnfolded  = Summarize_minusenu5.Get("unfolded_data1")

        data       = Summarize_minusenu5.Get("data_hist")
        reco       = Summarize_minusenu5.Get("reco_hist")
        truth      = Summarize_minusenu5.Get("truth_hist")

        HStatError = HUnfolded.Clone("HStatError")
        HSystTotal = HUnfolded.Clone("HSystTotal")
        Covariance = Summarize_minusenu5.Get("CovarianceMatrix_Iter1")


	if  (Name.find("enu")  != -1):
                    
            IDCovMatrix      = IdSF.Get(    "ElIDSys_Covariance_Iter1"  )   
            TrigCovMatrix    = TrigSF.Get(  "TrigSys_Covariance_Iter1"  )  
	    RecoCovMatrix    = RecoSF.Get(  "RecoSys_Covariance_Iter1"  ) 
  	    IsoCovMatrix     = IsoSF.Get(   "ElIsoSys_Covariance_Iter1" ) 
	    RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter1"   ) 
	    CalibCovMatrix   = Calib.Get(   "Calib_Covariance_Iter1"    ) 
    
    	if  (Name.find("munu")  != -1):
	
	    TrigCovMatrix    = TrigSF.Get(  "MuTrigSys_Covariance_Iter1" ) 
	    RecoCovMatrix    = RecoSF.Get(  "MuRecoSys_Covariance_Iter1" ) 
	    IsoCovMatrix     = IsoSF.Get(   "MuIsoSys_Covariance_Iter1"  )  
	    RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter1"	) 
                
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

        # define the Canvas output
        astyle.SetAtlasStyle()
        c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )

	Xs = (data.Integral()/Lum) / (truth.Integral()/reco.Integral())
        DataStat = Xs/sqrt( data.Integral() )
	LumSyst	 = 0.016*Xs 
	
        Syst = Xs*sqrt(Syst)/HUnfolded.Integral() 

	Syst1 = DataStat
	Syst2 = sqrt(DataStat*DataStat + Syst*Syst)        
	Syst3 = sqrt(DataStat*DataStat + Syst*Syst + LumSyst*LumSyst)

	print(Xs, DataStat, LumSyst, Syst)
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

        xxx  = array( 'f', [ 	Xs] )
        exxx = array( 'f', [ Syst3] )
        yyy  = array( 'f', [ 	Xs] )
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

        gr.SetMarkerColor( 4 )
        gr.SetMarkerStyle( 21 )
        gr.SetFillColor(4)
        gr.SetFillStyle(1001)

        grr.SetMarkerColor( 5 )
        grr.SetMarkerStyle( 21 )
        grr.SetFillColor(5)
        grr.SetFillStyle(1001)

        grrr.SetMarkerColor( 2 )
        grrr.SetMarkerStyle( 21 )
        grrr.SetFillColor(2)
        grrr.SetFillStyle(1001)

        grrr.Draw("ap2")
        grr.Draw("CP2")
        gr.Draw("CP2")


        # define the legend
        histsN = []
        gr.SetLineWidth(0)
        grr.SetLineWidth(0)
        grrr.SetLineWidth(0)
        gr.SetName("data \pm stat")
        grr.SetName("data \pm stat \pm Syst")
        grrr.SetName("data \pm stat \pm Syst \pm Lum")

        histsN.append(gr)
        histsN.append(grr)
        histsN.append(grrr)
        legendN = makeLegend(histsN,0.2, 0.6,0.5,0.78)
        legendN.Draw("same")

        astyle.ATLASLabel(0.2, 0.86, "Internal")
        utils.DrawText(0.2, 0.81, Indice)
        # TCanvas.Update() draws the frame, after which one can change it
        c1.Update()
        c1.Print("Output/CrossSection/Fiducial_CrossSection_"+Name+".pdf")

    def GetSummaringTable(self, Summarize_minusenu5, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum, Systbyb):

        data       = Summarize_minusenu5.Get("data_hist")
        reco       = Summarize_minusenu5.Get("reco_hist")
        truth      = Summarize_minusenu5.Get("truth_hist")
        mig_mat    = Summarize_minusenu5.Get("mig_hist")
	bg	   = Summarize_minusenu5.Get("Background_Total")

        # ***************************************************************************************************************************************************************************************
        # ********************************************************************************* Bin-by-Bin resultats ********************************************************************************
        # ***************************************************************************************************************************************************************************************

	# Cross section
        Xsbyb = ( (data.Integral()-bg.Integral())/Lum ) * ( truth.Integral() / reco.Integral() )
	
	# Stat uncertainties
	Statbyb = 100/sqrt(data.Integral())
	
        # ***************************************************************************************************************************************************************************************
	# ********************************************************************************* Unfolding resultats *********************************************************************************
	# ***************************************************************************************************************************************************************************************

	# define Unfolded:
        HUnfolded  = Summarize_minusenu5.Get("unfolded_data1")
	
	# Stat:        
	Covariance = Summarize_minusenu5.Get("CovarianceMatrix_Iter1")
	        
	Stat = 0
        for i in range(1, 1 + Covariance.GetNbinsX()):
            for j in range(1, 1+Covariance.GetNbinsX()):
                Stat = Stat + Covariance.GetBinContent(i, j)
	
	Stat = 100*sqrt(Stat)/HUnfolded.Integral()

	# Bias:
        CovBias      = Bias.Get("CovMatrix_Iter_1")
	        
	Bias = 0
        for i in range(1, 1 + CovBias.GetNbinsX()):
            for j in range(1, 1 + CovBias.GetNbinsX()):
                Bias = Bias + CovBias.GetBinContent(i, j)
        
        Bias = 100*sqrt(Bias)/HUnfolded.Integral()


	# define systematics:
        if  (Name.find("enu")  != -1):

            IDCovMatrix      = IdSF.Get(    "ElIDSys_Covariance_Iter1"  )
            TrigCovMatrix    = TrigSF.Get(  "ElTrigSys_Covariance_Iter1"  )
            RecoCovMatrix    = RecoSF.Get(  "ElRecoSys_Covariance_Iter1"  )
            IsoCovMatrix     = IsoSF.Get(   "ElIsoSys_Covariance_Iter1" )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter1"   )
            CalibCovMatrix   = Calib.Get(   "Calib_Covariance_Iter1"    )

        if  (Name.find("munu")  != -1):

            TrigCovMatrix    = TrigSF.Get(  "MuTrigSys_Covariance_Iter1" )
            RecoCovMatrix    = RecoSF.Get(  "MuRecoSys_Covariance_Iter1" )
            IsoCovMatrix     = IsoSF.Get(   "MuIsoSys_Covariance_Iter1"  )
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
	
	Syst = 100*sqrt(Syst)/HUnfolded.Integral()

	# define the cross section:
	truth     = Summarize_minusenu5.Get("truth_hist")
	truth_Acc = (Summarize_minusenu5.Get("mig_hist")).ProjectionY()

	Xs = ( HUnfolded.Integral()/Lum ) * ( truth.Integral() / truth_Acc.Integral() )

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
        latexFile.write("                                                           &    %s  \\\ \\hline \n"%(Indice))
        latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ (Unfolding)}         &    %5.3f   $\\pm$ %5.3f(Stat) $\\pm$ %5.3f(Syst) $\\pm$ %5.3f(Unf) $\\pm$ %5.3f(Lum)     \\\ \\hline \n" %( Xs,    Stat,    Syst,     Bias, 1.6))
        latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ (Bin-by-Bin)}        &    %5.3f   $\\pm$ %5.3f(Stat) $\\pm$ %5.3f(Syst) $\\pm$ %5.3f(Unf) $\\pm$ %5.3f(Lum)     \\\ \\hline \n" %( Xsbyb, Statbyb, Systbyb,  Bias, 1.6))
        latexFile.write("\\end{tabular}\n")
        latexFile.write("\\end{table}\n")
        latexFile.write("\\end{document}\n")
        latexFile.close()
	

	'''
    def GetDiffernetialXs(self, Summarize_minusenu5, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum):

        data       = Summarize_minusenu5.Get("data_hist")
        reco       = Summarize_minusenu5.Get("reco_hist")
        truth      = Summarize_minusenu5.Get("truth_hist")
        mig_mat    = Summarize_minusenu5.Get("mig_hist")
        bg         = Summarize_minusenu5.Get("Background_Total")

        # define Unfolded:
        HUnfolded  = Summarize_minusenu5.Get("unfolded_data4")

        # Stat:        
        Covariance = Summarize_minusenu5.Get("CovarianceMatrix_Iter4")

        # Bias:
        CovBias      = Bias.Get("CovMatrix_Iter_1")

        # define systematics:
        if  (Name.find("enu")  != -1):

            IDCovMatrix      = IdSF.Get(    "ElIDSys_Covariance_Iter4"  )
            TrigCovMatrix    = TrigSF.Get(  "TrigSys_Covariance_Iter4"  )
            RecoCovMatrix    = RecoSF.Get(  "RecoSys_Covariance_Iter4"  )
            IsoCovMatrix     = IsoSF.Get(   "ElIsoSys_Covariance_Iter4" )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter4"   )
            CalibCovMatrix   = Calib.Get(   "Calib_Covariance_Iter4"    )

        if  (Name.find("munu")  != -1):

            TrigCovMatrix    = TrigSF.Get(  "MuTrigSys_Covariance_Iter4" )
            RecoCovMatrix    = RecoSF.Get(  "MuRecoSys_Covariance_Iter4" )
            IsoCovMatrix     = IsoSF.Get(   "MuIsoSys_Covariance_Iter4"  )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter4"    )

        CovMatrix   = RecoCovMatrix.Clone("CovMatrix")

        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                if  (Name.find("enu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) )
                if  (Name.find("munu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )

	# define the Diff Xs:
	
        latexFile = open("Output/LatexTableau/Differential_Xs_"+Name+".tex","w+")
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

        latexFile.write("& \\multicolumn{4}{c|}{%s}  \\\ \\hline \\hline \n" %(Indice))
        latexFile.write("\\multicolumn{1}{|c|}{  range } & $d\sigma$/$d#eta_{l}$ [GeV]     & Stat uncertainty     & Unfolding bias     & Syst uncertainty        \\\ \\hline \\hline \n")

        i = 0
        while i < HUnfolded.GetNbinsX():
	    diffXs   = HUnfolded.GetBinContent(i+1)/(Lum*HUnfolded.GetBinWidth(i+1))
	    lowedge  = HUnfolded.GetBinLowEdge(i+1)
	    highEdge = HUnfolded.GetBinLowEdge(i+2)
	    Stat     = 100*sqrt(Covariance.GetBinContent( i+1, i+1 )) / HUnfolded.GetBinContent(i+1)
	    Bias     = 100*sqrt(CovBias.GetBinContent(    i+1, i+1 )) / HUnfolded.GetBinContent(i+1)
	    Syst     = 100*sqrt(CovMatrix.GetBinContent(  i+1, i+1 )) / HUnfolded.GetBinContent(i+1)
            latexFile.write("\\multicolumn{1}{|c|}{{[}%d,  %d{]}}  & %5.3f & %5.3f & %5.3f & %5.3f \\\ \\hline \n"  %(lowedge, highEdge, diffXs, Stat, Bias, Syst))

            i=i+1

        latexFile.write("\\end{tabular}\n")
        latexFile.write("\\end{table}\n")
        latexFile.write("\\end{document}\n")
        latexFile.close()
	'''

    def GetDiffernetialXsTest(self, Summarize_plusenu5, Summarize_plusmunu5, Energy, Indice, Name, Lum):

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


        for i in range(0, HUnfolded1.GetNbinsX() ):
            if(HUnfolded1.GetBinContent(i) != 0):
                DXs1.SetBinContent(i+1, Acctepatance1.GetBinContent(i+1)*HUnfolded1.GetBinContent(i+1) / (Lum*HUnfolded1.GetBinWidth(i+1)) )
                DXs2.SetBinContent(i+1, Acctepatance2.GetBinContent(i+1)*HUnfolded2.GetBinContent(i+1) / (Lum*HUnfolded2.GetBinWidth(i+1)) )

        for i in range(0, DXs1.GetNbinsX() ):
	     print("Bin ", i, "     ",DXs1.GetBinContent(i+1))

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

        DXs1.SetMarkerStyle(4)
        DXs2.SetMarkerStyle(4)

        DXs1.SetMarkerSize(1)
        DXs2.SetMarkerSize(1)

        DXs1.SetLineColor(2)
        DXs2.SetLineColor(4)

        DXs1.Draw()
        DXs2.Draw("same")

        astyle.ATLASLabel( 0.6, 0.86, "Internal")
        #utils.DrawText(    0.6, 0.81,  Indice)

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

        c1.Print("Output/CrossSection/Differential_Xs_"+Name+".pdf")
	


	
    def GetDiffernetialXsPlot(self, Summarize_minusenu5, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum):

        data       = Summarize_minusenu5.Get("data_hist")
        reco       = Summarize_minusenu5.Get("reco_hist")
        truth      = Summarize_minusenu5.Get("truth_hist")
        mig_mat    = Summarize_minusenu5.Get("mig_hist")
        bg         = Summarize_minusenu5.Get("Background_Total")

        # define Unfolded:
        HUnfolded  = Summarize_minusenu5.Get("unfolded_data4")

        # Stat:        
        Covariance = Summarize_minusenu5.Get("CovarianceMatrix_Iter4")

        # Bias:
        CovBias      = Bias.Get("CovMatrix_Iter_1")

        # define systematics:
        if  (Name.find("enu")  != -1):

            IDCovMatrix      = IdSF.Get(    "ElIDSys_Covariance_Iter4"  )
            TrigCovMatrix    = TrigSF.Get(  "ElTrigSys_Covariance_Iter4"  )
            RecoCovMatrix    = RecoSF.Get(  "ElRecoSys_Covariance_Iter4"  )
            IsoCovMatrix     = IsoSF.Get(   "ElIsoSys_Covariance_Iter4" )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter4"   )
            CalibCovMatrix   = Calib.Get(   "Calib_Covariance_Iter4"    )

        if  (Name.find("munu")  != -1):

            TrigCovMatrix    = TrigSF.Get(  "MuTrigSys_Covariance_Iter4" )
            RecoCovMatrix    = RecoSF.Get(  "MuRecoSys_Covariance_Iter4" )
            IsoCovMatrix     = IsoSF.Get(   "MuIsoSys_Covariance_Iter4"  )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter4"    )

        CovMatrix   = RecoCovMatrix.Clone("CovMatrix")

        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                if  (Name.find("enu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) )
                if  (Name.find("munu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )

	# add STAT, Bias, Syst:
	Uncer = truth.Clone("Uncer")
	
	for i in range(1, 1 + Uncer.GetNbinsX()):
	    if  (Name.find("enu")  != -1):
		Uncer.SetBinContent(i,  TrigCovMatrix.GetBinContent(i,i) + RecoCovMatrix.GetBinContent(i,i) + IsoCovMatrix.GetBinContent(i,i) + IDCovMatrix.GetBinContent(i,i) + CalibCovMatrix.GetBinContent(i,i) + RecoilCovMatrix.GetBinContent(i,i) + Covariance.GetBinContent(i,i) + CovBias.GetBinContent(i,i) )
            if  (Name.find("munu")  != -1):
                Uncer.SetBinContent(i,  TrigCovMatrix.GetBinContent(i,i) + RecoCovMatrix.GetBinContent(i,i) + IsoCovMatrix.GetBinContent(i,i) + RecoilCovMatrix.GetBinContent(i,i) + Covariance.GetBinContent(i,i) + CovBias.GetBinContent(i,i) )
	
        for i in range(1, 1 + Uncer.GetNbinsX()):
	    if(HUnfolded.GetBinContent(i) != 0):
		Uncer.SetBinContent(i, 1 + 100*sqrt(Uncer.GetBinContent(i)) / HUnfolded.GetBinContent(i) )
	    	Uncer.SetBinError(i, 0)

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

        c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )
        astyle.SetAtlasStyle()

        #define the pad
        pad1 = root.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.Draw()
        pad1.SetLogy()
        pad1.SetFillStyle(4000)
        pad1.cd()
	
	DXs.GetXaxis().SetRangeUser(-2.,2.)
	DXs.SetStats(0)
	DXs.SetTitle("")
	DXs.GetYaxis().SetTitle("d#sigma / d#eta_{l} ")
	DXs.GetXaxis().SetTitle("#eta_{l}")
	DXs.GetYaxis().SetTitleOffset(0.9)
	DXs.GetXaxis().SetTitleOffset(1.4)
        DXs.GetYaxis().SetTitleSize(0.07)
        DXs.GetYaxis().SetLabelSize(0.08)

	DXs.SetLineWidth(2)
	DMC.SetLineWidth(2)
	
        DXs.SetMarkerColor(2)
        DMC.SetMarkerColor(4)

        DXs.SetMarkerStyle(4)
        DMC.SetMarkerStyle(4)

        DXs.SetMarkerSize(1)
        DMC.SetMarkerSize(1)

	DXs.SetLineColor(2)
        DMC.SetLineColor(4)

	DXs.Draw()
	DMC.Draw("same")

        #gr.Draw("ap2")
        #grb.Draw("CP2")
        #grs.Draw("CP2")
	astyle.ATLASLabel( 0.6, 0.86, "Internal")
        utils.DrawText(    0.6, 0.81,  Indice)

	DXs.SetName(" Unfolded data")
	DMC.SetName(" PowhegPythia8")
	hists = []
        hists.append(DXs)
        hists.append(DMC)
	
        legendN = makeLegend(hists,0.6, 0.5,0.92,0.75)
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

	Ratio1.GetXaxis().SetRangeUser(-2.,2.)
	Ratio1.SetLineWidth(2)
        Ratio1.SetLineColor(2)
        Ratio1.Draw("same")

        line1 = ROOT.TLine(-2,0.95,2,0.95)
        line2 = ROOT.TLine(-2,1.05,2,1.05)
        line1.SetLineStyle(2)
        line2.SetLineStyle(2)
	
	line1.Draw("same")
	line2.Draw("same")

        c1.Print("Output/CrossSection/Differential_Xs_"+Name+".pdf")



    def GetDiffernetialXs(self, Summarize_minusenu5, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum):

        data       = Summarize_minusenu5.Get("data_hist")
        reco       = Summarize_minusenu5.Get("reco_hist")
        truth      = Summarize_minusenu5.Get("truth_hist")
        mig_mat    = Summarize_minusenu5.Get("mig_hist")
        bg         = Summarize_minusenu5.Get("Background_Total")

        Hist1      = reco.Clone("Hist1")
        Hist2      = reco.Clone("Hist2")
        Hist3      = reco.Clone("Hist3")


        Acceptance = Summarize_minusenu5.Get("Acceptance_hist")

        # define Unfolded:
        Unfolded  = Summarize_minusenu5.Get("unfolded_data2")

        if  (Name.find("enu")  != -1):
        	binning = [-2.5, -2.18, -1.95, -1.74, -1.52, -1.37, -1.05, -0.84, -0.63, -0.42, -0.21, 0.0, 0.21, 0.42, 0.63, 0.84, 1.05, 1.37, 1.52, 1.74, 1.95, 2.18, 2.5]
        if  (Name.find("munu")  != -1):
		binning = [-2.4, -1.918, -1.348, -1.1479, -1.05,  -0.908, -0.476, 0, 0.476, 0.908, 1.05, 1.1479, 1.348, 1.918, 2.4]

        HUnfolded = ROOT.TH1D("h",";m [GeV];Events/ 50 GeV;;",len(binning)-1,array('d',binning))

        for i in range(1, 1 + Unfolded.GetNbinsX() ):
            HUnfolded.SetBinContent(i, Unfolded.GetBinContent(i))
            HUnfolded.SetBinError(i,   Unfolded.GetBinError(i))

        # Stat:        
        Covariance = Summarize_minusenu5.Get("CovarianceMatrix_Iter2")

        # Bias:
        CovBias      = Bias.Get("CovMatrix_Iter_2")

        # define systematics:
        if  (Name.find("enu")  != -1):

            IDCovMatrix      = IdSF.Get(    "ElIDSys_Covariance_Iter2"  )
            TrigCovMatrix    = TrigSF.Get(  "elTrigSys_Covariance_Iter2"  )
            RecoCovMatrix    = RecoSF.Get(  "elRecoSys_Covariance_Iter2"  )
            IsoCovMatrix     = IsoSF.Get(   "ElIsoSys_Covariance_Iter2" )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter2"   )
            CalibCovMatrix   = Calib.Get(   "Calib_Covariance_Iter2"    )

        if  (Name.find("munu")  != -1):

            TrigCovMatrix    = TrigSF.Get(  "MuTrigSys_Covariance_Iter2" )
            RecoCovMatrix    = RecoSF.Get(  "MuRecoSys_Covariance_Iter2" )
            #IsoCovMatrix     = IsoSF.Get(   "MuIsoSys_Covariance_Iter2"  )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter2"    )

        CovMatrix   = RecoCovMatrix.Clone("CovMatrix")

        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                if  (Name.find("enu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) )
                if  (Name.find("munu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )
                    #CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )

        for i in range(0, Hist1.GetNbinsX()):
            if(HUnfolded.GetBinContent(i+1) != 0):
                Hist1.SetBinContent(i+1, 100*sqrt(Covariance.GetBinContent( i+1, i+1 )) / HUnfolded.GetBinContent(i+1))
                Hist2.SetBinContent(i+1, 100*sqrt(CovMatrix.GetBinContent(  i+1, i+1 )) / HUnfolded.GetBinContent(i+1))
                Hist3.SetBinContent(i+1, 100*sqrt(CovBias.GetBinContent(    i+1, i+1 )) / HUnfolded.GetBinContent(i+1))
                Hist1.SetBinError(i+1, 0)
                Hist2.SetBinError(i+1, 0)
                Hist3.SetBinError(i+1, 0)
            if(HUnfolded.GetBinContent(i+1) == 0):
                Hist1.SetBinContent(i+1, 0)
                Hist2.SetBinContent(i+1, 0)
                Hist3.SetBinContent(i+1, 0)
                Hist1.SetBinError(i+1, 0)
                Hist2.SetBinError(i+1, 0)
                Hist3.SetBinError(i+1, 0)


        # define the Diff Xs:

        latexFile = open("Output/LatexTableau/Differential_Xs_"+Name+".tex","w+")
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

        latexFile.write("& \\multicolumn{4}{c|}{%s}  \\\ \\hline \\hline \n" %(Indice))
        latexFile.write("\\multicolumn{1}{|c|}{  range } & $d\sigma$/$d\\eta_{l}$   & Stat uncertainty     & Unfolding bias     & Syst uncertainty        \\\ \\hline \\hline \n")

        i = 0
        while i <= HUnfolded.GetNbinsX():
            if(HUnfolded.GetBinContent(i+1) != 0 and Acceptance.GetBinContent(i+1) != 0):
      		print(i+1, HUnfolded.GetNbinsX(), HUnfolded.GetBinLowEdge(i+1))
	      	diffXs   = HUnfolded.GetBinContent(i+1)/(Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
            	lowedge  = HUnfolded.GetBinLowEdge(i+1)
            	highEdge = HUnfolded.GetBinLowEdge(i+2)
            	Stat     = 100*sqrt(Covariance.GetBinContent( i+1, i+1 )) / (HUnfolded.GetBinContent(i+1))
            	Bias     = 100*sqrt(CovBias.GetBinContent(    i+1, i+1 )) / (HUnfolded.GetBinContent(i+1))
            	Syst     = 100*sqrt(CovMatrix.GetBinContent(  i+1, i+1 )) / (HUnfolded.GetBinContent(i+1))
            	latexFile.write("\\multicolumn{1}{|c|}{{[}%5.2f,  %5.2f{]}}  & %5.3f & %5.3f & %5.3f & %5.3f \\\ \\hline \n"  %(lowedge, highEdge, diffXs, Stat, Bias, Syst))
            i=i+1

        latexFile.write("\\end{tabular}\n")
        latexFile.write("\\end{table}\n")
        latexFile.write("\\end{document}\n")
        latexFile.close()

        astyle.SetAtlasStyle()
        c1 = TCanvas("c1", "C1", 1200, 800)
        Hist1.SetStats(0)
        Hist1.SetTitle("")
        Hist1.GetYaxis().SetTitle("d#sigma / d#eta_{l} ")
        Hist1.GetXaxis().SetTitle("#eta_{l}")
        Hist1.GetYaxis().SetTitleOffset(0.9)
        Hist1.GetXaxis().SetTitleOffset(1.4)
        Hist1.GetYaxis().SetTitleSize(0.06)
        Hist1.GetYaxis().SetLabelSize(0.05)
        Hist1.GetXaxis().SetTitleSize(0.05)
        Hist1.GetXaxis().SetLabelSize(0.05)

        Hist1.GetXaxis().SetRangeUser(25, 80)
        Hist1.GetYaxis().SetRangeUser(0, 1.5)
        Hist1.SetLineWidth(2)
        Hist2.SetLineWidth(2)
        Hist3.SetLineWidth(2)
        Hist1.SetLineColor(1)
        Hist2.SetLineColor(2)
        Hist3.SetLineColor(4)

        Hist1.SetName(" Stat uncertainties")
        Hist2.SetName(" Syst uncertainties")
        Hist3.SetName(" Unfolded bias")

        hists = []
        hists.append(Hist1)
        hists.append(Hist2)
        hists.append(Hist3)
        legendN = makeLegend(hists,0.2, 0.5,0.5,0.7)

        Hist1.Draw()
        Hist2.Draw("same")
        Hist3.Draw("same")
        legendN.Draw("same")
        astyle.ATLASLabel( 0.2, 0.86, "Internal")
        utils.DrawText(    0.2, 0.81,  Indice)

        c1.Print("Output/LatexTableau/Differential_Xs_"+Name+"_plot.pdf")

	print(Hist1.GetNbinsX())

    def GetDiffernetialXsPlotN(self, Summarize_minusenu5, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum):

        data       = Summarize_minusenu5.Get("data_hist")
        reco       = Summarize_minusenu5.Get("reco_hist")
        truth      = Summarize_minusenu5.Get("truth_hist")
        mig_mat    = Summarize_minusenu5.Get("mig_hist")
        bg         = Summarize_minusenu5.Get("Background_Total")

        # define Unfolded:
        HUnfolded  = Summarize_minusenu5.Get("unfolded_data2")

	# define the Xs:        
	Xs = ( HUnfolded.Integral()/Lum ) * ( truth.Integral() / (mig_mat.ProjectionY()).Integral() )

        # Stat:        
        Covariance = Summarize_minusenu5.Get("CovarianceMatrix_Iter2")

        # Bias:
        CovBias      = Bias.Get("CovMatrix_Iter_2")

        # define systematics:
        if  (Name.find("enu")  != -1):

            IDCovMatrix      = IdSF.Get(    "ElIDSys_Covariance_Iter2"  )
            TrigCovMatrix    = TrigSF.Get(  "elTrigSys_Covariance_Iter2"  )
            RecoCovMatrix    = RecoSF.Get(  "elRecoSys_Covariance_Iter2"  )
            IsoCovMatrix     = IsoSF.Get(   "elIsoSys_Covariance_Iter2" )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter2"   )
            CalibCovMatrix   = Calib.Get(   "Calib_Covariance_Iter2"    )

        if  (Name.find("munu")  != -1):

            TrigCovMatrix    = TrigSF.Get(  "MuTrigSys_Covariance_Iter2" )
            RecoCovMatrix    = RecoSF.Get(  "MuRecoSys_Covariance_Iter2" )
            IsoCovMatrix     = IsoSF.Get(   "MuIsoSys_Covariance_Iter2"  )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter2"    )

        CovMatrix   = RecoCovMatrix.Clone("CovMatrix")

        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                if  (Name.find("enu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) )
                if  (Name.find("munu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )

        # define the acceptance correction:
        truth_noMiss = (mig_mat.ProjectionY()).Clone("truth_noMiss")
        Acctepatance = truth.Clone("Acctepatance")
        Acctepatance.Divide(truth_noMiss)

        DXs = truth.Clone("DXs")
        DMC = truth.Clone("DXs")

        for i in range(0, HUnfolded.GetNbinsX() ):
            if(HUnfolded.GetBinContent(i) != 0):
	            DXs.SetBinContent(i+1, (Acctepatance.GetBinContent(i+1)*HUnfolded.GetBinContent(i+1) / (Lum*HUnfolded.GetBinWidth(i+1)) ) / Xs )
        	    DMC.SetBinContent(i+1, (truth.GetBinContent(i+1) / (Lum*truth.GetBinWidth(i+1)) ) / Xs )
        	    DMC.SetBinError(i+1, 0)

        # define the ratio:
        Ratio1 = DXs.Clone("Ratio1")
        Ratio1.Divide(DMC)

        c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )
        astyle.SetAtlasStyle()

        #define the pad
        pad1 = root.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.Draw()
        pad1.SetLogy()
        pad1.SetFillStyle(4000)
        pad1.cd()

        DXs.GetXaxis().SetRangeUser(-2.,2.)
        DXs.SetStats(0)
        DXs.SetTitle("")
        DXs.GetYaxis().SetTitle("1/d#sigma . d#sigma / d#eta_{l} ")
        DXs.GetXaxis().SetTitle("#eta_{l} [GeV]")
        DXs.GetYaxis().SetTitleOffset(0.9)
        DXs.GetXaxis().SetTitleOffset(1.4)
        DXs.GetYaxis().SetTitleSize(0.07)
        DXs.GetYaxis().SetLabelSize(0.08)

        DXs.SetLineWidth(2)
        DMC.SetLineWidth(2)

        DXs.SetMarkerColor(2)
        DMC.SetMarkerColor(4)

        DXs.SetMarkerStyle(4)
        DMC.SetMarkerStyle(4)

        DXs.SetMarkerSize(1)
        DMC.SetMarkerSize(1)

        DXs.SetLineColor(2)
        DMC.SetLineColor(4)

        DXs.Draw()
        DMC.Draw("same")

        astyle.ATLASLabel( 0.6, 0.86, "Internal")
        utils.DrawText(    0.6, 0.81,  Indice)

        DXs.SetName(" Unfolded data")
        DMC.SetName(" PowhegPythia8")
        hists = []
        hists.append(DXs)
        hists.append(DMC)

        legendN = makeLegend(hists,0.6, 0.5,0.92,0.75)
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
        Ratio1.GetXaxis().SetTitle("#eta_{l} [GeV]")
        Ratio1.GetXaxis().SetLabelSize(0.15)
        Ratio1.GetXaxis().SetTitleSize(0.16)
        Ratio1.GetXaxis().SetTitleOffset(1.1)

        Ratio1.GetYaxis().SetTitleOffset(0.5)
        Ratio1.GetYaxis().SetTitleSize(0.11)
        Ratio1.GetYaxis().SetLabelSize(0.15)
        Ratio1.GetYaxis().CenterTitle()

        Ratio1.GetXaxis().SetRangeUser(-2.,2.)
        Ratio1.SetLineWidth(2)
        Ratio1.SetLineColor(2)
        Ratio1.Draw("same")

        line1 = ROOT.TLine(-2,0.95,2,0.95)
        line2 = ROOT.TLine(-2,1.05,2,1.05)
        line1.SetLineStyle(2)
        line2.SetLineStyle(2)

        line1.Draw("same")
        line2.Draw("same")

        c1.Print("Output/CrossSection/Differential_Xs_"+Name+"_Normalised.pdf")






































