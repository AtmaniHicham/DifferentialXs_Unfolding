#!/usr/bin/env python
# -*-coding:Latin-1 -*
import numpy as np
import matplotlib.pyplot as plt
import pylab

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

class OpitmisationStudy1B:
        """Classe représentant une personne"""

        def __init__(self):
            """Initialisation """

        def StatStudy(self, inputFile, IterMin, IterMax, Nbin, channel, Var):
            binContent  = []    # some bin content
            StatError   = []    # with correlation

            compteur  = IterMin
            while compteur <= IterMax:
                value   = 0
                Error   = 0

                CovMatrix        = inputFile.Get("CovarianceMatrix_Iter"+ str(compteur))
                data_Unfolded    = inputFile.Get("unfolded_data"+ str(compteur))
		
		value = value + data_Unfolded.GetBinContent(Nbin)
		Error = Error + CovMatrix.GetBinContent(Nbin,Nbin)
	
                binContent.append(value)
                StatError.append(Error)
                compteur += 1

            i = 0
            while i < IterMax:
		if( binContent[i] != 0):
                    StatError[i]   = 100*sqrt(StatError[i])/binContent[i]
                    i=i+1

            x=np.linspace(1,IterMax,IterMax)
            pylab.plot(x, StatError,   marker='o', label =  '[' + str( (data_Unfolded.GetBinLowEdge(Nbin)) ) + ',' + str( (data_Unfolded.GetBinLowEdge(Nbin+1))) + ']')
            pylab.xlabel('Iterations');
            pylab.ylabel('Stat [%]');
            pylab.legend(loc='center right')
            pylab.savefig("Output/"+channel+"/Stat_1B_Optimisation.pdf")

        def BiasStudy(self, inputFile, Bias, IterMin, IterMax, Nbin, channel, Var):
            binContent  = []    # some bin content
            StatError   = []    # with correlation

            compteur  = IterMin
            while compteur <= IterMax:
                value   = 0
                Error   = 0

                CovMatrix        = Bias.Get("CovMatrix_Iter_"+ str(compteur))
                data_Unfolded    = inputFile.Get("unfolded_data"+ str(compteur))

		print(type(data_Unfolded))

		value = value + data_Unfolded.GetBinContent(Nbin)
		Error = Error + CovMatrix.GetBinContent(Nbin,Nbin)

                binContent.append(value)
                StatError.append(Error)
                compteur += 1

            i = 0
            while i < IterMax:
                StatError[i]   = 100*sqrt(StatError[i])/binContent[i]
                i=i+1

            x=np.linspace(1,IterMax,IterMax)
	    pylab.plot(x, StatError,   marker='o', label =  '[' + str( (data_Unfolded.GetBinLowEdge(Nbin)) ) + ',' + str( (data_Unfolded.GetBinLowEdge(Nbin+1))) + ']')
            pylab.xlabel('Iterations');
            pylab.ylabel('Bias [%]');
            pylab.legend(loc='center right')
            pylab.savefig("Output/"+channel+"/Bias_1B_Optimisation.pdf")


        def EffSystematicStudy(self, inputFile, Syst_TrigSys, Syst_RecoSys, Syst_ElIsoSys, Syst_ElIDSys, IterMin, IterMax, Nbin, channel, Var):

            binContent  = []    # some bin content
            StatError   = []    # with correlation

            compteur  = IterMin
            while compteur <= IterMax:
                value   = 0
                Error   = 0

		if  (channel.find("enu")  != -1):        
		    IDCovMatrix      = Syst_ElIDSys.Get(  "ElIDSys_Covariance_Iter"   + str(compteur))
		    TrigCovMatrix    = Syst_TrigSys.Get(  "ElTrigSys_Covariance_Iter"   + str(compteur))
                    RecoCovMatrix    = Syst_RecoSys.Get(  "ElRecoSys_Covariance_Iter"   + str(compteur))
                    IsoCovMatrix     = Syst_ElIsoSys.Get( "ElIsoSys_Covariance_Iter"  + str(compteur))

                if  (channel.find("munu")  != -1):
                    TrigCovMatrix    = Syst_TrigSys.Get(  "MuTrigSys_Covariance_Iter" + str(compteur))
                    RecoCovMatrix    = Syst_RecoSys.Get(  "MuRecoSys_Covariance_Iter" + str(compteur))
                    IsoCovMatrix     = Syst_ElIsoSys.Get( "MuIsoSys_Covariance_Iter"  + str(compteur))

		CovMatrix = RecoCovMatrix.Clone("CovMatrix")

		for i in range(1, 1 + CovMatrix.GetNbinsX()):
		    for j in range(1, 1+CovMatrix.GetNbinsX()):
	                if  (channel.find("enu")  != -1):
                         	CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j)+RecoCovMatrix.GetBinContent(i,j)+IsoCovMatrix.GetBinContent(i,j)+IDCovMatrix.GetBinContent(i,j) ) 		
	                if  (channel.find("munu")  != -1):
                                CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j)+RecoCovMatrix.GetBinContent(i,j)+IsoCovMatrix.GetBinContent(i,j) )

                data_Unfolded    = inputFile.Get("unfolded_data"+ str(compteur))

		value = value + data_Unfolded.GetBinContent(Nbin)
		Error = Error + CovMatrix.GetBinContent(Nbin, Nbin)

                binContent.append(value)
                StatError.append(Error)
                compteur += 1

            i = 0
            while i < IterMax:
                StatError[i]   = 100*sqrt(StatError[i])/binContent[i]
                i=i+1

            x=np.linspace(1,IterMax,IterMax)
            pylab.plot(x, StatError,   marker='o', label =  '[' + str( (data_Unfolded.GetBinLowEdge(Nbin)) ) + ',' + str( (data_Unfolded.GetBinLowEdge(Nbin+1))) + ']')
            pylab.xlabel('Iterations');
            pylab.ylabel('Eff Syst [%]');
            pylab.legend(loc='center right')
            pylab.savefig("Output/"+channel+"/EffSyst_1B_Optimisation.pdf")


        def CalibSystematicStudy(self, inputFile, Calib, IterMin, IterMax, Nbin, channel, Var):
            binContent  = []    # some bin content
            StatError   = []    # with correlation

            compteur  = IterMin
            while compteur <= IterMax:
                value   = 0
                Error   = 0

                CovMatrix        = Calib.Get("Calib_Covariance_Iter"+ str(compteur))
                data_Unfolded    = inputFile.Get("unfolded_data"+ str(compteur))

		value = value + data_Unfolded.GetBinContent(Nbin)
		Error = Error + CovMatrix.GetBinContent(Nbin,Nbin)

                binContent.append(value)
                StatError.append(Error)
                compteur += 1

            i = 0
            while i < IterMax:
                StatError[i]   = 100*sqrt(StatError[i])/binContent[i]
                i=i+1

            x=np.linspace(1,IterMax,IterMax)
            pylab.plot(x, StatError,   marker='o', label =  '[' + str( (data_Unfolded.GetBinLowEdge(Nbin)) ) + ',' + str( (data_Unfolded.GetBinLowEdge(Nbin+1))) + ']')
            pylab.xlabel('Iterations');
            pylab.ylabel('Calib [%]');
            pylab.legend(loc='center right')
            pylab.savefig("Output/"+channel+"/Calib_1B_Optimisation.pdf")


        def RecoilSystematicStudy(self, inputFile, Recoil, IterMin, IterMax, Nbin, channel, Var):
            binContent  = []    # some bin content
            StatError   = []    # with correlation

            compteur  = IterMin
            while compteur <= IterMax:
                value   = 0
                Error   = 0

                CovMatrix        = Recoil.Get("Recoil_Covariance_Iter"+ str(compteur))
                data_Unfolded    = inputFile.Get("unfolded_data"+ str(compteur))

		value = value + data_Unfolded.GetBinContent(Nbin)
		Error = Error + CovMatrix.GetBinContent(Nbin,Nbin)

                binContent.append(value)
                StatError.append(Error)
                compteur += 1

            i = 0
            while i < IterMax:
                StatError[i]   = 100*sqrt(StatError[i])/binContent[i]
                i=i+1

            x=np.linspace(1,IterMax,IterMax)
            pylab.plot(x, StatError,   marker='o', label =  '[' + str( (data_Unfolded.GetBinLowEdge(Nbin)) ) + ',' + str( (data_Unfolded.GetBinLowEdge(Nbin+1))) + ']')
            pylab.xlabel('Iterations');
            pylab.ylabel('Recoil [%]');
            pylab.legend(loc='center right')
            pylab.savefig("Output/"+channel+"/Recoil_1B_Optimisation.pdf")


        def TotalSystematicStudy(self, inputFile, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, IterMin, IterMax, Nbin, channel, Var):
            binContent  = []    # some bin content
            StatError   = []    # with correlation

            compteur  = IterMin
            while compteur <= IterMax:
                value   = 0
                Error   = 0

                if  (channel.find("enu")  != -1):
                    IDCovMatrix      = IdSF.Get(  "ElIDSys_Covariance_Iter"   + str(compteur))
                    TrigCovMatrix    = TrigSF.Get( "ElTrigSys_Covariance_Iter"   + str(compteur))
                    RecoCovMatrix    = RecoSF.Get(  "ElRecoSys_Covariance_Iter"   + str(compteur))
                    IsoCovMatrix     = IsoSF.Get( "ElIsoSys_Covariance_Iter"  + str(compteur))
		    RecoilCovMatrix  = Recoil.Get("Recoil_Covariance_Iter"+ str(compteur))
		    CalibCovMatrix   = Calib.Get("Calib_Covariance_Iter"+ str(compteur))

                if  (channel.find("munu")  != -1):
                    TrigCovMatrix    = TrigSF.Get(  "MuTrigSys_Covariance_Iter" + str(compteur))
                    RecoCovMatrix    = RecoSF.Get(  "MuRecoSys_Covariance_Iter" + str(compteur))
                    IsoCovMatrix     = IsoSF.Get( "MuIsoSys_Covariance_Iter"  + str(compteur))
		    RecoilCovMatrix  = Recoil.Get("Recoil_Covariance_Iter"+ str(compteur))
   
		StatCovMatrix	 =  inputFile.Get("CovarianceMatrix_Iter"+ str(compteur))
		BiasCovMatrix    =  Bias.Get("CovMatrix_Iter_"+ str(compteur))		

		CovMatrix	 = BiasCovMatrix.Clone("CovMatrix")

                for i in range(1, 1 + CovMatrix.GetNbinsX()):
                    for j in range(1, 1+CovMatrix.GetNbinsX()):
			if  (channel.find("enu")  != -1):
                        	CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) + StatCovMatrix.GetBinContent(i,j) )
                	if  (channel.find("munu")  != -1):
                                CovMatrix.SetBinContent(i, j, TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) + StatCovMatrix.GetBinContent(i,j) )

                data_Unfolded    = inputFile.Get("unfolded_data"+ str(compteur))

		value = value + data_Unfolded.GetBinContent(Nbin)
		Error = Error + CovMatrix.GetBinContent(Nbin, Nbin)

                binContent.append(value)
                StatError.append(Error)
                compteur += 1

            i = 0
            while i < IterMax:
                StatError[i]   = 100*sqrt(StatError[i])/binContent[i]
                i=i+1

            x=np.linspace(1,IterMax,IterMax)
            pylab.plot(x, StatError,   marker='o', label =  '[' + str( (data_Unfolded.GetBinLowEdge(Nbin)) ) + ',' + str( (data_Unfolded.GetBinLowEdge(Nbin+1))) + ']')
            pylab.xlabel('Iterations');
            pylab.ylabel('Stat + Systematics [%]');
            pylab.legend(loc='center right')
            pylab.savefig("Output/"+channel+"/SysteTotal_1B_Optimisation.pdf")

        def TotalStudy(self, inputFile, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, IterMin, IterMax, Nbin, channel, Var):
            binContent  = []    # some bin content
            StatError   = []    # with correlation

            compteur  = IterMin
            while compteur <= IterMax:
                value   = 0
                Error   = 0

                if  (channel.find("enu")  != -1):
                    IDCovMatrix      = IdSF.Get(  "ElIDSys_Covariance_Iter"   + str(compteur))
                    TrigCovMatrix    = TrigSF.Get(  "ElTrigSys_Covariance_Iter"   + str(compteur))
                    RecoCovMatrix    = RecoSF.Get(  "ElRecoSys_Covariance_Iter"   + str(compteur))
                    IsoCovMatrix     = IsoSF.Get( "ElIsoSys_Covariance_Iter"  + str(compteur))
                    RecoilCovMatrix  = Recoil.Get("Recoil_Covariance_Iter"+ str(compteur))
                    CalibCovMatrix   = Calib.Get("Calib_Covariance_Iter"+ str(compteur))

                if  (channel.find("munu")  != -1):
                    TrigCovMatrix    = TrigSF.Get(  "MuTrigSys_Covariance_Iter" + str(compteur))
                    RecoCovMatrix    = RecoSF.Get(  "MuRecoSys_Covariance_Iter" + str(compteur))
                    IsoCovMatrix     = IsoSF.Get( "MuIsoSys_Covariance_Iter"  + str(compteur))
                    RecoilCovMatrix  = Recoil.Get("Recoil_Covariance_Iter"+ str(compteur))

                StatCovMatrix    =  inputFile.Get("CovarianceMatrix_Iter"+ str(compteur))
                BiasCovMatrix    =  Bias.Get("CovMatrix_Iter_"+ str(compteur))

                CovMatrix        = BiasCovMatrix.Clone("CovMatrix")

                for i in range(1, 1 + CovMatrix.GetNbinsX()):
                    for j in range(1, 1+CovMatrix.GetNbinsX()):
                        if  (channel.find("enu")  != -1):
                                CovMatrix.SetBinContent(i, j, BiasCovMatrix.GetBinContent(i,j) + TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) + StatCovMatrix.GetBinContent(i,j) )
                        if  (channel.find("munu")  != -1):
                                CovMatrix.SetBinContent(i, j, BiasCovMatrix.GetBinContent(i,j) + TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) + StatCovMatrix.GetBinContent(i,j) )

                data_Unfolded    = inputFile.Get("unfolded_data"+ str(compteur))

                value = value + data_Unfolded.GetBinContent(Nbin)
                Error = Error + CovMatrix.GetBinContent(Nbin, Nbin)

                binContent.append(value)
                StatError.append(Error)
                compteur += 1

            i = 0
            while i < IterMax:
                StatError[i]   = 100*sqrt(StatError[i])/binContent[i]
                i=i+1

            x=np.linspace(1,IterMax,IterMax)
            pylab.plot(x, StatError,   marker='o', label =  '[' + str( (data_Unfolded.GetBinLowEdge(Nbin)) ) + ',' + str( (data_Unfolded.GetBinLowEdge(Nbin+1))) + ']')
            pylab.xlabel('Iterations');
            pylab.ylabel('Bias + Stat + Systematics [%]');
            pylab.legend(loc='center right')
            pylab.savefig("Output/"+channel+"/Total_1B_Optimisation.pdf")


