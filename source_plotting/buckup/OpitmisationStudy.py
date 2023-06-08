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

class OpitmisationStudy:
        """Classe représentant une personne"""

        def __init__(self):
            """Initialisation """

        def StatStudy(self, inputFile, IterMin, IterMax, mTmin, mTmax, channel, Var):
            binContent  = []    # some bin content
            StatError   = []    # with correlation

            compteur  = IterMin
            while compteur <= IterMax:
                value   = 0
                Error   = 0

                CovMatrix        = inputFile.Get("CovarianceMatrix_Iter"+ str(compteur))
                data_Unfolded    = inputFile.Get("unfolded_data"+ str(compteur))
		
		jmin = data_Unfolded.FindBin(mTmin)
                jmax = data_Unfolded.FindBin(mTmax)
	        kmin = data_Unfolded.FindBin(mTmin)                
		kmax = data_Unfolded.FindBin(mTmax)
	
                #print("jmin: %d, jmax:%d, kmin:%d, kmax:%d"%(jmin, jmax, kmin, kmax))

		for j in range(jmin, jmax):
                    value   = value   + data_Unfolded.GetBinContent(j)
		    for k in range(kmin, kmax):
                        Error = Error + CovMatrix.GetBinContent(j,k)
	
                binContent.append(value)
                StatError.append(Error)
                compteur += 1

            i = 0
            while i < IterMax:
                StatError[i]   = 100*sqrt(StatError[i])/binContent[i]
                i=i+1

            x=np.linspace(1,IterMax,IterMax)
            pylab.plot(x, StatError,   marker='o', label= str(mTmin) + ' < ' + Var + ' < ' + str(mTmax) )
            pylab.xlabel('Iterations');
            pylab.ylabel('Stat [%]');
            pylab.legend(loc='center right')
            pylab.savefig("Output/"+channel+"/Stat_Optimisation.pdf")

        def BiasStudy(self, inputFile, Bias, IterMin, IterMax, mTmin, mTmax, channel, Var):
            binContent  = []    # some bin content
            StatError   = []    # with correlation

            compteur  = IterMin
            while compteur <= IterMax:
                value   = 0
                Error   = 0

                CovMatrix        = Bias.Get("CovMatrix_Iter_"+ str(compteur))
                data_Unfolded    = inputFile.Get("unfolded_data"+ str(compteur))

                jmin = data_Unfolded.FindBin(mTmin)
                jmax = data_Unfolded.FindBin(mTmax)
                kmin = data_Unfolded.FindBin(mTmin)
                kmax = data_Unfolded.FindBin(mTmax)


                #print("jmin: %d, jmax:%d, kmin:%d, kmax:%d"%(jmin, jmax, kmin, kmax))
                for j in range(jmin, jmax):
                    value   = value   + data_Unfolded.GetBinContent(j)
                    for k in range(kmin, kmax):
                        print(" Iteration %d, bin j:%d, bin k:%d"%(compteur, j, k))
                        Error = Error + CovMatrix.GetBinContent(j,k)

                binContent.append(value)
                StatError.append(Error)
                compteur += 1

            i = 0
            while i < IterMax:
                StatError[i]   = 100*sqrt(StatError[i])/binContent[i]
                i=i+1

            x=np.linspace(1,IterMax,IterMax)
            pylab.plot(x, StatError,   marker='o', label= str(mTmin) + ' < ' + Var + ' < ' + str(mTmax) )
            pylab.xlabel('Iterations');
            pylab.ylabel('Bias [%]');
            pylab.legend(loc='center right')
            pylab.savefig("Output/"+channel+"/Bias_Optimisation.pdf")


        def EffSystematicStudy(self, inputFile, Syst_TrigSys, Syst_RecoSys, Syst_ElIsoSys, Syst_ElIDSys, IterMin, IterMax, mTmin, mTmax, channel, Var):

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

                jmin = data_Unfolded.FindBin(mTmin)
                jmax = data_Unfolded.FindBin(mTmax)
                kmin = data_Unfolded.FindBin(mTmin)
                kmax = data_Unfolded.FindBin(mTmax)

		#print("jmin: %d, jmax:%d, kmin:%d, kmax:%d"%(jmin, jmax, kmin, kmax))
                for j in range(jmin, jmax):
                    value   = value   + data_Unfolded.GetBinContent(j)
                    for k in range(kmin, kmax):
			print(" Iteration %d, bin j:%d, bin k:%d"%(compteur, j, k))
                        Error = Error + CovMatrix.GetBinContent(j,k)

                binContent.append(value)
                StatError.append(Error)
                compteur += 1

            i = 0
            while i < IterMax:
                StatError[i]   = 100*sqrt(StatError[i])/binContent[i]
                i=i+1

            x=np.linspace(1,IterMax,IterMax)
            pylab.plot(x, StatError,   marker='o', label= str(mTmin) + ' < ' + Var + ' < ' + str(mTmax) )
            pylab.xlabel('Iterations');
            pylab.ylabel('Eff Syst [%]');
            pylab.legend(loc='center right')
            pylab.savefig("Output/"+channel+"/EffSyst_Optimisation.pdf")


        def CalibSystematicStudy(self, inputFile, Calib, IterMin, IterMax, mTmin, mTmax, channel, Var):
            binContent  = []    # some bin content
            StatError   = []    # with correlation

            compteur  = IterMin
            while compteur <= IterMax:
                value   = 0
                Error   = 0

                CovMatrix        = Calib.Get("Calib_Covariance_Iter"+ str(compteur))
                data_Unfolded    = inputFile.Get("unfolded_data"+ str(compteur))

                jmin = data_Unfolded.FindBin(mTmin)
                jmax = data_Unfolded.FindBin(mTmax)
                kmin = data_Unfolded.FindBin(mTmin)
                kmax = data_Unfolded.FindBin(mTmax)


                #print("jmin: %d, jmax:%d, kmin:%d, kmax:%d"%(jmin, jmax, kmin, kmax))
                for j in range(jmin, jmax):
                    value   = value   + data_Unfolded.GetBinContent(j)
                    for k in range(kmin, kmax):
                        print(" Iteration %d, bin j:%d, bin k:%d"%(compteur, j, k))
                        Error = Error + CovMatrix.GetBinContent(j,k)

                binContent.append(value)
                StatError.append(Error)
                compteur += 1

            i = 0
            while i < IterMax:
                StatError[i]   = 100*sqrt(StatError[i])/binContent[i]
                i=i+1

            x=np.linspace(1,IterMax,IterMax)
            pylab.plot(x, StatError,   marker='o', label= str(mTmin) + ' < ' + Var + ' < ' + str(mTmax) )
            pylab.xlabel('Iterations');
            pylab.ylabel('Calib [%]');
            pylab.legend(loc='center right')
            pylab.savefig("Output/"+channel+"/Calib_Optimisation.pdf")


        def RecoilSystematicStudy(self, inputFile, Recoil, IterMin, IterMax, mTmin, mTmax, channel, Var):
            binContent  = []    # some bin content
            StatError   = []    # with correlation

            compteur  = IterMin
            while compteur <= IterMax:
                value   = 0
                Error   = 0

                CovMatrix        = Recoil.Get("Recoil_Covariance_Iter"+ str(compteur))
                data_Unfolded    = inputFile.Get("unfolded_data"+ str(compteur))

                jmin = data_Unfolded.FindBin(mTmin)
                jmax = data_Unfolded.FindBin(mTmax)
                kmin = data_Unfolded.FindBin(mTmin)
                kmax = data_Unfolded.FindBin(mTmax)


                #print("jmin: %d, jmax:%d, kmin:%d, kmax:%d"%(jmin, jmax, kmin, kmax))
                for j in range(jmin, jmax):
                    value   = value   + data_Unfolded.GetBinContent(j)
                    for k in range(kmin, kmax):
                        print(" Iteration %d, bin j:%d, bin k:%d"%(compteur, j, k))
                        Error = Error + CovMatrix.GetBinContent(j,k)

                binContent.append(value)
                StatError.append(Error)
                compteur += 1

            i = 0
            while i < IterMax:
                StatError[i]   = 100*sqrt(StatError[i])/binContent[i]
                i=i+1

            x=np.linspace(1,IterMax,IterMax)
            pylab.plot(x, StatError,   marker='o', label= str(mTmin) + ' < ' + Var + ' < ' + str(mTmax) )
            pylab.xlabel('Iterations');
            pylab.ylabel('Recoil [%]');
            pylab.legend(loc='center right')
            pylab.savefig("Output/"+channel+"/Recoil_Optimisation.pdf")


        def TotalSystematicStudy(self, inputFile, Bias, TrigSF, RecoSF, IsoSF, IdSF, Calib, Recoil, IterMin, IterMax, mTmin, mTmax, channel, Var):
            binContent  = []    # some bin content
            StatError   = []    # with correlation

            compteur  = IterMin
            while compteur <= IterMax:
                value   = 0
                Error   = 0

                if  (channel.find("enu")  != -1):
                    IDCovMatrix      = IdSF.Get(    "ElIDSys_Covariance_Iter"   + str(compteur))
                    TrigCovMatrix    = TrigSF.Get(  "ElTrigSys_Covariance_Iter"   + str(compteur))
                    RecoCovMatrix    = RecoSF.Get(  "ElRecoSys_Covariance_Iter"   + str(compteur))
                    IsoCovMatrix     = IsoSF.Get(   "ElIsoSys_Covariance_Iter"  + str(compteur))
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
                        	CovMatrix.SetBinContent(i, j, BiasCovMatrix.GetBinContent(i,j) + TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + IDCovMatrix.GetBinContent(i,j) + CalibCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) + StatCovMatrix.GetBinContent(i,j) )
                	if  (channel.find("munu")  != -1):
                                CovMatrix.SetBinContent(i, j, BiasCovMatrix.GetBinContent(i,j) + TrigCovMatrix.GetBinContent(i,j) + RecoCovMatrix.GetBinContent(i,j) + IsoCovMatrix.GetBinContent(i,j) + RecoilCovMatrix.GetBinContent(i,j) + StatCovMatrix.GetBinContent(i,j) )

                data_Unfolded    = inputFile.Get("unfolded_data"+ str(compteur))

                jmin = data_Unfolded.FindBin(mTmin)
                jmax = data_Unfolded.FindBin(mTmax)
                kmin = data_Unfolded.FindBin(mTmin)
                kmax = data_Unfolded.FindBin(mTmax)


                #print("jmin: %d, jmax:%d, kmin:%d, kmax:%d"%(jmin, jmax, kmin, kmax))
                for j in range(jmin, jmax):
                    value   = value   + data_Unfolded.GetBinContent(j)
                    for k in range(kmin, kmax):
                        print(" Iteration %d, bin j:%d, bin k:%d"%(compteur, j, k))
                        Error = Error + CovMatrix.GetBinContent(j,k)

                binContent.append(value)
                StatError.append(Error)
                compteur += 1

            i = 0
            while i < IterMax:
                StatError[i]   = 100*sqrt(StatError[i])/binContent[i]
                i=i+1

            x=np.linspace(1,IterMax,IterMax)
            pylab.plot(x, StatError,   marker='o', label= str(mTmin) + ' < ' + Var + ' < ' + str(mTmax) )
            pylab.xlabel('Iterations');
            pylab.ylabel('Bias + Stat + Systematics [%]');
            pylab.legend(loc='center right')
            pylab.savefig("Output/"+channel+"/Total_Optimisation.pdf")

        def OptimisationResum(self,inputFile, Bias, IterMin, IterMax, mTmin, mTmax, channel):

            binContent  = []    # some bin content
            StatErrorWC = []    # without correlation
            StatError   = []    # with correlation

            compteur  = IterMin
            while compteur <= IterMax:
                value   = 0
                Error   = 0
                ErrorWC = 0

                CovMatrix        = inputFile.Get("Covariance_Matrix"+ str(compteur))
                data_Unfolded    = inputFile.Get("Unfolded_data_iteration_"+ str(compteur))
                j = mTmin
                while j < mTmax:
                    value   = value   + data_Unfolded.GetBinContent(j)
                    ErrorWC = ErrorWC + CovMatrix.GetBinContent(j,j)
                    k = mTmin
                    while k < mTmax:
                        Error = Error + CovMatrix.GetBinContent(k,j)
                        k += 1
                    j += 1
                binContent.append(value)
                StatErrorWC.append(ErrorWC)
                StatError.append(Error)
                compteur += 1

            binContentB  = []    # some bin content
            StatErrorWCB = []    # without correlation
            StatErrorB   = []    # with correlation

            compteur  = IterMin
            while compteur <= IterMax:
                valueB   = 0
                ErrorB   = 0
                ErrorWCB = 0

                CovMatrixB        = Bias.Get("Covariance_Bias_Iteration"+ str(compteur))
                data_UnfoldedB    = inputFile.Get("Unfolded_data_iteration_"+ str(compteur))
                j = mTmin
                while j < mTmax:
                    valueB   = valueB   + data_UnfoldedB.GetBinContent(j)
                    ErrorWCB = ErrorWCB + CovMatrixB.GetBinContent(j,j)
                    k = mTmin
                    while k < mTmax:
                        ErrorB = ErrorB + CovMatrixB.GetBinContent(k,j)
                        k += 1
                    j += 1
                binContentB.append(valueB)
                StatErrorWCB.append(ErrorWCB)
                StatErrorB.append(ErrorB)
                compteur += 1

            i=0
            while i < 20:
                StatError[i]   = StatErrorB[i]  + StatError[i]
                StatErrorWC[i] = StatErrorWC[i] + StatErrorWCB[i]
                i=i+1

            i=0
            while i < 20:
                StatError[i]   = 100*sqrt(StatError[i])/binContent[i]
                StatErrorWC[i] = 100*sqrt(StatErrorWC[i])/binContent[i]
                i=i+1

            x=np.linspace(1,20,20)
            pylab.plot(x, StatErrorWC, marker='o', label=' Stat+Bias without correlation - '+channel)
            #pylab.plot(x, StatError,   marker='o', label=' Stat+Bias with correlation')
            pylab.xlabel('Iterations');
            pylab.ylabel('Uncertainties [%]');
            pylab.legend(loc='left left')
            pylab.savefig("Output/"+channel+"/Stat+Bias_Optimisation.pdf")
