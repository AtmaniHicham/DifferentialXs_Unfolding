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
from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TAttFill, TMatrixD, THilbertMatrixD, TDecompSVD, TGraphErrors, gRandom, TText

def makeLegend(hists, xmin, ymin, xmax, ymax):
    legend = root.TLegend(xmin, ymin, xmax, ymax)
    legend.SetTextSize(0.05)
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

    def GetFiducialXs(self, Summarize_minusenu6, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum):
              
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
	print("Summarising Table")
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

	print(type(RecoCovMatrix))
	'''
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

	Hist1 	   = truth.Clone("Hist1")
        Hist2      = truth.Clone("Hist2")
        Hist3      = truth.Clone("Hist3")

        Acceptance = Summarize_minusenu5.Get("Acceptance_hist")

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
            #IsoCovMatrix     = IsoSF.Get(   "MuIsoSys_Covariance_Iter4"  )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter4"    )

        CovMatrix   = RecoCovMatrix.Clone("CovMatrix")

        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                if  (Name.find("enu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) )
                if  (Name.find("munu")  != -1):
                    #CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )

	for i in range(0, Hist1.GetNbinsX()):
	    Hist1.SetBinContent(i+1, 100*sqrt(Covariance.GetBinContent( i+1, i+1 )) / HUnfolded.GetBinContent(i+1))
	    Hist2.SetBinContent(i+1, 100*sqrt(CovMatrix.GetBinContent(  i+1, i+1 )) / HUnfolded.GetBinContent(i+1))
	    Hist3.SetBinContent(i+1, 100*sqrt(CovBias.GetBinContent(    i+1, i+1 )) / HUnfolded.GetBinContent(i+1))
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
        latexFile.write("\\multicolumn{1}{|c|}{  range } & $d\sigma$/$dp^{T}_{l}$ [GeV]     & Stat uncertainty     & Unfolding bias     & Syst uncertainty        \\\ \\hline \\hline \n")

        i = 0
        while i < HUnfolded.GetNbinsX():
	    diffXs   = HUnfolded.GetBinContent(i+1)/(Lum*HUnfolded.GetBinWidth(i+1)*Acceptance.GetBinContent(i+1))
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

        astyle.SetAtlasStyle()
	c1 = Canvas("c1", "C1", 1200, 800)
	Hist1.Draw()
	Hist2.Draw("same")        
	Hist3.Draw("same")
	c1.Print("Output/LatexTableau/Differential_Xs_"+Name+"_plot.pdf")

	
    def GetDiffernetialXsPlot(self, Summarize_minusenu5, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum):

        data       = Summarize_minusenu5.Get("data_hist")
        reco       = Summarize_minusenu5.Get("reco_hist")
        truth      = Summarize_minusenu5.Get("truth_hist")
        mig_mat    = Summarize_minusenu5.Get("mig_hist")
        bg         = Summarize_minusenu5.Get("Background_Total")

        # define Unfolded:
        HUnfolded  = Summarize_minusenu5.Get("unfolded_data4")

        # Stat:        
        CovarianceM = Summarize_minusenu5.Get("CovarianceMatrix_Iter4")
        #Covariance  = Summarize_minusenu5.Get("CovarianceMatrix_Iter4")

        # Bias:
        CovBias      = Bias.Get("CovMatrix_Iter_1")

	# define the covariance
	Covariance = CovBias.Clone("Covariance")
	for i in range(0, truth.GetNbinsX()):
	    for j in range(0,  truth.GetNbinsX()):
		Covariance.SetBinContent(i+1, j+1, CovarianceM[i][j])

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
            IsoCovMatrix     = IsoSF.Get(   "MuRecoSys_Covariance_Iter4"  )
            #IsoCovMatrix     = IsoSF.Get(   "MuIsoSys_Covariance_Iter4"  )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter4"    )

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
		Uncer.SetBinContent(i, 1  )
	    	Uncer.SetBinError(i, 100*sqrt(Uncer.GetBinContent(i)) / HUnfolded.GetBinContent(i) )

        UncerUp   = Uncer.Clone("UncerUp")
        UncerDown = Uncer.Clone("UncerDown")

        for i in range(1, 1 + Uncer.GetNbinsX()):
            UncerUp.SetBinContent(i, 1+Uncer.GetBinError(i))
            UncerUp.SetBinError(i, 0)
            UncerDown.SetBinContent(i, 1-Uncer.GetBinError(i))
            UncerDown.SetBinError(i, 0)

	# define the acceptance correction:
	truth_noMiss = (mig_mat.ProjectionY()).Clone("truth_noMiss")
	Acctepatance = truth.Clone("Acctepatance")
	Acctepatance.Divide(truth_noMiss)

	DXs = truth.Clone("DXs")	
	DMC = truth.Clone("DXs")

	for i in range(0, HUnfolded.GetNbinsX() ):
            if(HUnfolded.GetBinContent(i+1) != 0):
		print(i+1, Acctepatance.GetBinContent(i+1)*HUnfolded.GetBinContent(i+1) / (Lum*HUnfolded.GetBinWidth(i+1)) )
	    	DXs.SetBinContent(i+1, Acctepatance.GetBinContent(i+1)*HUnfolded.GetBinContent(i+1) / (Lum*HUnfolded.GetBinWidth(i+1)) )
	    	DMC.SetBinContent(i+1, truth.GetBinContent(i+1) / (Lum*truth.GetBinWidth(i+1)) )
		DXs.SetBinError(i+1, Acctepatance.GetBinContent(i+1)*(HUnfolded.GetBinError(i+1)) / (Lum*HUnfolded.GetBinWidth(i+1)) )
	    	DMC.SetBinError(i+1, (truth.GetBinError(i+1))  / (Lum*truth.GetBinWidth(i+1)) )

        #for i in range(0, DXs.GetNbinsX() ):
        #    DXs.SetBinContent(i+1,   gRandom.Poisson(2)+DMC.GetBinContent(i+1))

	
	# define the ratio:
	Ratio1 = DXs.Clone("Ratio1")
        Ratio1.Divide(DMC)

	gROOT.SetStyle("ATLAS");
        astyle.SetAtlasStyle()
        c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )

        #define the pad
        pad1 = root.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.Draw()
        pad1.SetLogy()
        pad1.SetFillStyle(4000)
        pad1.cd()
	
        DMC.SetStats(0)
        DMC.SetTitle("")
        DXs.SetLineWidth(0)
        DXs.SetLineColor(1)
        DXs.SetMarkerStyle(8)
        DXs.SetMarkerSize(0.8)
        DXs.SetMarkerColor(1)
        DMC.SetLineWidth(2)
        DMC.SetLineColor(2)

        DMC.GetXaxis().SetLabelSize(0.05)
        DMC.GetXaxis().SetTitleSize(0.045)
        DMC.GetYaxis().SetLabelSize(0.045)
        DMC.GetYaxis().SetTitleSize(0.075)
        DMC.GetXaxis().SetRangeUser(30,60)
        DMC.GetYaxis().SetTitle("d#sigma / d#eta_{l}[pb]")
        DMC.GetXaxis().SetTitle("#eta_{l}")
        DMC.GetXaxis().SetTitleOffset(1.5)

	print(DMC.GetMinimum())
        print(DMC.GetMaximum())

	DMC.SetMinimum(70)
	DMC.SetMaximum(1500)
        DMC.GetYaxis().SetLabelSize(0.065)


	DMC.Draw("H")
        DXs.Draw("psame")

	for i in range(0, DMC.GetNbinsX()):
	   print(i+1, DMC.GetBinContent(i+1), DXs.GetBinContent(i+1))

	l = root.TLatex()
    	l.SetNDC()
   	l.SetTextSize(0.09)
    	l.DrawLatex(0.6, 0.77, r"#bf{#it{ATLAS}} " + "Internal")
        l.DrawLatex(0.6, 0.67, Indice)

	#astyle.ATLASLabel( 0.7, 0.77, "Internal")
        #utils.DrawText(    0.7, 0.70,  Indice)

	DXs.SetName(" Unfolded data")
	DMC.SetName(" PowhegPythia8")
	UncerUp.SetName(" total uncertainty")
	hists = []
        hists.append(DXs)
        hists.append(DMC)
	hists.append(UncerUp)

        legendN = makeLegend(hists,0.4, 0.2,0.6,0.45)
        legendN.Draw("same")
        c1.Update()        
	c1.cd()

        pad2 = root.TPad("pad2", "pad2", 0, 0.06, 1, 0.32)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.4)
        pad2.Draw()
        pad2.SetFillStyle(4000)
        pad2.cd()

	UncerUp.SetStats(0)
        UncerUp.GetYaxis().SetRangeUser(0.91,1.09)
        UncerUp.GetYaxis().SetNdivisions(5)
        UncerUp.GetYaxis().SetTitle("Pred./Data")
        UncerUp.GetXaxis().SetTitle("#eta_{l}")

        UncerUp.GetXaxis().SetLabelSize(0.15)
        UncerUp.GetXaxis().SetTitleSize(0.15)
        UncerUp.GetYaxis().SetLabelSize(0.15)
        UncerUp.GetYaxis().SetTitleSize(0.15)
        UncerUp.GetXaxis().SetTitleOffset(1.2)
        UncerUp.GetYaxis().SetTitleOffset(0.4)
	UncerUp.GetXaxis().SetRangeUser(30,60)

        Ratio1.SetMarkerStyle(8)
        Ratio1.SetMarkerSize(0.8)
        Ratio1.SetMarkerColor(1)
        Ratio1.SetLineWidth(2)
        Ratio1.SetLineStyle(1)
        Ratio1.SetLineColor(1)

        UncerUp.SetLineColorAlpha(1,1)
        UncerDown.SetLineColorAlpha(1,1)
        UncerUp.SetFillColorAlpha(2,1)
        UncerDown.SetFillColorAlpha(0,1)
	UncerUp.SetLineWidth(0)
	UncerUp.SetMarkerSize(0)
        UncerUp.SetMarkerColor(2)
        UncerUp.Draw("same")
        UncerDown.Draw("same")        
	Ratio1.Draw("Psame")
	
        line1 = ROOT.TLine(30,0.95,60,0.95)
        line2 = ROOT.TLine(30,1.05,60,1.05)
        line3 = ROOT.TLine(30,1.00,60,1.00)
        line1.SetLineStyle(2)
        line2.SetLineStyle(2)        
	line3.SetLineStyle(2)

	line1.Draw("same")
	line2.Draw("same")
        line3.Draw("same")
        c1.Print("Output/CrossSection/Differential_Xs_"+Name+".pdf")


    def GetDiffernetialXs(self, Summarize_minusenu5, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum, Niter):

        data       = Summarize_minusenu5.Get("data_hist")
        reco       = Summarize_minusenu5.Get("reco_hist")
        truth      = Summarize_minusenu5.Get("truth_hist")
        mig_mat    = Summarize_minusenu5.Get("mig_hist")
        bg         = Summarize_minusenu5.Get("Background_Total")

        Hist1      = truth.Clone("Hist1")
        Hist2      = truth.Clone("Hist2")
        Hist3      = truth.Clone("Hist3")

        # define Unfolded:
        HUnfolded  = Summarize_minusenu5.Get("unfolded_data"+str(Niter))

        # Stat:        
        Covariance = Summarize_minusenu5.Get("CovarianceMatrix_Iter"+str(Niter))

        # Bias:
        CovBias      = Bias.Get("CovMatrix_Iter_"+str(Niter))


        # define systematics:
        if  (Name.find("enu")  != -1):

            IDCovMatrix      = IdSF.Get(    "ElIDSys_Covariance_Iter"  +str(Niter) )
            TrigCovMatrix    = TrigSF.Get(  "ElTrigSys_Covariance_Iter"+str(Niter) )
            RecoCovMatrix    = RecoSF.Get(  "ElRecoSys_Covariance_Iter"+str(Niter) )
            IsoCovMatrix     = IsoSF.Get(   "ElIsoSys_Covariance_Iter" +str(Niter) )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter"   +str(Niter) )
            CalibCovMatrix   = Calib.Get(   "Calib_Covariance_Iter"    +str(Niter) )

        if  (Name.find("munu")  != -1):

            TrigCovMatrix    = TrigSF.Get(  "MuTrigSys_Covariance_Iter"+str(Niter) )
            RecoCovMatrix    = RecoSF.Get(  "MuRecoSys_Covariance_Iter"+str(Niter) )
            #IsoCovMatrix     = IsoSF.Get(   "MuIsoSys_Covariance_Iter" +str(Niter) )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter"   +str(Niter) )

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
        latexFile.write("\\multicolumn{1}{|c|}{  range } & $d\sigma$/$dp^{T}_{l}$ [GeV]     & Stat uncertainty     & Unfolding bias     & Syst uncertainty        \\\ \\hline \\hline \n")

        i = 25
        while i < 61:
            if(HUnfolded.GetBinContent(i) != 0):
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


        astyle.SetAtlasStyle()
        c1 = TCanvas("c1", "C1", 1200, 800)
        Hist1.SetStats(0)
        Hist1.SetTitle("")
        Hist1.GetYaxis().SetTitle("d#sigma / dp_{l}^{T} ")
        Hist1.GetXaxis().SetTitle("p_{l}^{T} [GeV]")
        Hist1.GetYaxis().SetTitleOffset(0.9)
        Hist1.GetXaxis().SetTitleOffset(1.4)
        Hist1.GetYaxis().SetTitleSize(0.06)
        Hist1.GetYaxis().SetLabelSize(0.05)
        Hist1.GetXaxis().SetTitleSize(0.05)
        Hist1.GetXaxis().SetLabelSize(0.05)

	Hist1.GetXaxis().SetRangeUser(25, 80)
	Hist1.GetYaxis().SetRangeUser(0, 2.5)
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




    def GetDiffernetialXsPlotN(self, Summarize_minusenu5, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, Energy, Indice, Name, Lum):


        data       = Summarize_minusenu5.Get("data_hist")
        reco       = Summarize_minusenu5.Get("reco_hist")
        truth      = Summarize_minusenu5.Get("truth_hist")
        mig_mat    = Summarize_minusenu5.Get("mig_hist")
        bg         = Summarize_minusenu5.Get("Background_Total")

        # define Unfolded:
        HUnfolded  = Summarize_minusenu5.Get("unfolded_data4")

        # define the Xs:        
        Xs = ( HUnfolded.Integral()/Lum ) * ( truth.Integral() / (mig_mat.ProjectionY()).Integral() )

        XsTruth = ( truth.Integral()/Lum ) #* ( truth.Integral() / (mig_mat.ProjectionY()).Integral() )


        # Stat:        
        CovarianceM = Summarize_minusenu5.Get("CovarianceMatrix_Iter4")

        # Bias:
        CovBias      = Bias.Get("CovMatrix_Iter_1")

        # define the covariance
        Covariance = CovBias.Clone("Covariance")
        for i in range(0, 14):
            for j in range(0, 14):
                Covariance.SetBinContent(i+1, j+1, CovarianceM[i][j])



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
            IsoCovMatrix     = IsoSF.Get(   "MuRecoSys_Covariance_Iter4"  )
            #IsoCovMatrix     = IsoSF.Get(   "MuIsoSys_Covariance_Iter4"  )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter4"    )

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
                Uncer.SetBinContent(i, 1  )
                Uncer.SetBinError(i, 100*sqrt(Uncer.GetBinContent(i)) / HUnfolded.GetBinContent(i) )

        UncerUp   = Uncer.Clone("UncerUp")
        UncerDown = Uncer.Clone("UncerDown")

        for i in range(1, 1 + Uncer.GetNbinsX()):
            UncerUp.SetBinContent(i, 1+Uncer.GetBinError(i))
            UncerUp.SetBinError(i, 0)
            UncerDown.SetBinContent(i, 1-Uncer.GetBinError(i))
            UncerDown.SetBinError(i, 0)

        # define the acceptance correction:
        truth_noMiss = (mig_mat.ProjectionY()).Clone("truth_noMiss")
        Acctepatance = truth.Clone("Acctepatance")
        Acctepatance.Divide(truth_noMiss)

        DXs = truth.Clone("DXs")
        DMC = truth.Clone("DXs")

        for i in range(0, HUnfolded.GetNbinsX() ):
            if(HUnfolded.GetBinContent(i+1) != 0):
                DXs.SetBinContent(i+1, Acctepatance.GetBinContent(i+1)*HUnfolded.GetBinContent(i+1) / (Lum*HUnfolded.GetBinWidth(i+1)) / Xs )
                DMC.SetBinContent(i+1, truth.GetBinContent(i+1) / (Lum*truth.GetBinWidth(i+1)) / XsTruth )
                DXs.SetBinError(i+1, Acctepatance.GetBinContent(i+1)*HUnfolded.GetBinError(i+1) / (Lum*HUnfolded.GetBinWidth(i+1)) / Xs  )
                DMC.SetBinError(i+1, truth.GetBinError(i+1)  / (Lum*truth.GetBinWidth(i+1)) / XsTruth )

        # define the ratio:
        Ratio1 = DXs.Clone("Ratio1")
        Ratio1.Divide(DMC)

        gROOT.SetStyle("ATLAS");
        astyle.SetAtlasStyle()
        c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )

        #define the pad
        pad1 = root.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.Draw()
        pad1.SetLogy()
        pad1.SetFillStyle(4000)
        pad1.cd()

        DMC.SetStats(0)
        DMC.SetTitle("")
        DXs.SetLineWidth(0)
        DXs.SetLineColor(1)
        DXs.SetMarkerStyle(8)
        DXs.SetMarkerSize(0.8)
        DXs.SetMarkerColor(1)
        DMC.SetLineWidth(2)
        DMC.SetLineColor(2)

        DMC.GetXaxis().SetLabelSize(0.05)
        DMC.GetXaxis().SetTitleSize(0.045)
        DMC.GetYaxis().SetLabelSize(0.045)
        DMC.GetYaxis().SetTitleSize(0.075)
        DMC.GetXaxis().SetRangeUser(30,60)
        DMC.GetYaxis().SetTitle("1 / #sigma d#sigma / d#eta_{l}[pb]")
        DMC.GetXaxis().SetTitle("#eta_{l}^")
        DMC.GetXaxis().SetTitleOffset(1.5)

        DMC.SetMinimum(0.02)
        DMC.SetMaximum(0.8)
        DMC.GetYaxis().SetLabelSize(0.065)



        DMC.Draw("H")
        DXs.Draw("psame")

        l = root.TLatex()
        l.SetNDC()
        l.SetTextSize(0.09)
        l.DrawLatex(0.6, 0.77, r"#bf{#it{ATLAS}} " + "Internal")
        l.DrawLatex(0.6, 0.67, Indice)

        #astyle.ATLASLabel( 0.7, 0.77, "Internal")
        #utils.DrawText(    0.7, 0.70,  Indice)

        DXs.SetName(" Unfolded data")
        DMC.SetName(" PowhegPythia8")
        UncerUp.SetName(" total uncertainty")
        hists = []
        hists.append(DXs)
        hists.append(DMC)
        hists.append(UncerUp)

        legendN = makeLegend(hists,0.4, 0.2,0.6,0.45)
        legendN.Draw("same")
        c1.Update()
        c1.cd()

        pad2 = root.TPad("pad2", "pad2", 0, 0.06, 1, 0.32)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.4)
        pad2.Draw()
        pad2.SetFillStyle(4000)
        pad2.cd()

        UncerUp.GetYaxis().SetRangeUser(0.91,1.09)
        UncerUp.GetYaxis().SetNdivisions(5)
        UncerUp.GetYaxis().SetTitle("Pred./Data")
        UncerUp.GetXaxis().SetTitle("#eta_{l}")

        UncerUp.GetXaxis().SetLabelSize(0.15)
        UncerUp.GetXaxis().SetTitleSize(0.15)
        UncerUp.GetYaxis().SetLabelSize(0.15)
        UncerUp.GetYaxis().SetTitleSize(0.15)
        UncerUp.GetXaxis().SetTitleOffset(1.2)
        UncerUp.GetYaxis().SetTitleOffset(0.4)
        UncerUp.GetXaxis().SetRangeUser(30,60)

        Ratio1.SetMarkerStyle(8)
        Ratio1.SetMarkerSize(0.8)
        Ratio1.SetMarkerColor(1)
        Ratio1.SetLineWidth(2)
        Ratio1.SetLineStyle(1)
        Ratio1.SetLineColor(1)

        UncerUp.SetLineColorAlpha(1,1)
        UncerDown.SetLineColorAlpha(1,1)
        UncerUp.SetFillColorAlpha(2,1)
        UncerDown.SetFillColorAlpha(0,1)
        UncerUp.SetLineWidth(0)
        UncerUp.SetMarkerSize(0)
        UncerUp.SetMarkerColor(2)
        UncerUp.Draw("same")
        UncerDown.Draw("same")
        Ratio1.Draw("Psame")

        line1 = ROOT.TLine(30,0.95,60,0.95)
        line2 = ROOT.TLine(30,1.05,60,1.05)
        line3 = ROOT.TLine(30,1.00,60,1.00)
        line1.SetLineStyle(2)
        line2.SetLineStyle(2)
        line3.SetLineStyle(2)

        line1.Draw("same")
        line2.Draw("same")
        line3.Draw("same")
        c1.Print("Output/CrossSection/Differential_Xs_"+Name+"_Normalised.pdf")




	'''
        data       = Summarize_minusenu5.Get("data_hist")
        reco       = Summarize_minusenu5.Get("reco_hist")
        truth      = Summarize_minusenu5.Get("truth_hist")
        mig_mat    = Summarize_minusenu5.Get("mig_hist")
        bg         = Summarize_minusenu5.Get("Background_Total")

        # define Unfolded:
        HUnfolded  = Summarize_minusenu5.Get("unfolded_data4")

	# define the Xs:        
	Xs = ( HUnfolded.Integral()/Lum ) * ( truth.Integral() / (mig_mat.ProjectionY()).Integral() )

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
            #IsoCovMatrix     = IsoSF.Get(   "MuIsoSys_Covariance_Iter4"  )
            RecoilCovMatrix  = Recoil.Get(  "Recoil_Covariance_Iter4"    )

        CovMatrix   = RecoCovMatrix.Clone("CovMatrix")

        for i in range(1, 1 + CovMatrix.GetNbinsX()):
            for j in range(1, 1+CovMatrix.GetNbinsX()):
                if  (Name.find("enu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) )
                if  (Name.find("munu")  != -1):
                    CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )
                    #CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j)  )

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

        DXs.GetXaxis().SetRangeUser(30,60)
        DXs.SetStats(0)
        DXs.SetTitle("")
        DXs.GetYaxis().SetTitle("1/d#sigma . d#sigma / dp_{l}^{T} ")
        DXs.GetXaxis().SetTitle("p_{l}^{T} [GeV]")
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
        Ratio1.GetXaxis().SetTitle("p_{l}^{T} [GeV]")
        Ratio1.GetXaxis().SetLabelSize(0.15)
        Ratio1.GetXaxis().SetTitleSize(0.16)
        Ratio1.GetXaxis().SetTitleOffset(1.1)

        Ratio1.GetYaxis().SetTitleOffset(0.5)
        Ratio1.GetYaxis().SetTitleSize(0.11)
        Ratio1.GetYaxis().SetLabelSize(0.15)
        Ratio1.GetYaxis().CenterTitle()

        Ratio1.GetXaxis().SetRangeUser(30,60)
        Ratio1.SetLineWidth(3)
        Ratio1.SetLineColor(1)
	
        Ratio1.Draw("same")

        line1 = ROOT.TLine(30,0.95,60,0.95)
        line2 = ROOT.TLine(30,1.05,60,1.05)
        line1.SetLineStyle(2)
        line2.SetLineStyle(2)

        line1.Draw("same")
        line2.Draw("same")

        c1.Print("Output/CrossSection/Differential_Xs_"+Name+"_Normalised.pdf")

	'''




































