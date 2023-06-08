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
    legend.SetTextSize(0.03)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetBorderSize(0)
    for hist in hists:
        legend.AddEntry(hist, hist.GetName())
    return legend


class SystVariations:
    """Classe représentant une personne"""

    def __init__(self):
        """Constructeur de notre classe"""


    def CompareSyst(self, Summarize, IdSF, IsoSF, RecoSF, TrigSF, Recoil, Calib, channel, Indice):

        hists = []
        Colorcompteur = 1

        if  (channel.find("enu") != -1):

            hist = IdSF.Get('ElIDSys_Systematics_Iter1')
            hist.SetName('Id Sf ')
            hists.append(hist)

            hist = IsoSF.Get('ElIsoSys_Systematics_Iter1')
            hist.SetName('Iso Sf ')
            hists.append(hist)

            hist = RecoSF.Get('ElRecoSys_Systematics_Iter1')
            hist.SetName('reco Sf ')
            hists.append(hist)

            hist = TrigSF.Get('ElTrigSys_Systematics_Iter1')
            hist.SetName('Trigger Sf ')
            hists.append(hist)

            Unfolded = Summarize.Get('unfolded_data1')
            CovCalib = Calib.Get('Calib_Covariance_Iter1')
            hist = Unfolded.Clone("hist")
            for k in range(1, 1 + hist.GetNbinsX()):
                if(Unfolded.GetBinContent(k) != 0):
                    hist.SetBinContent(k, 100*sqrt(CovCalib.GetBinContent(k,k))/Unfolded.GetBinContent(k))
                    hist.SetBinError(k, 0)
            hist.SetName('Calibration')
            hists.append(hist)

            RecoilCalib = Recoil.Get('Recoil_Covariance_Iter1')
            hist = Unfolded.Clone("hist")
            for k in range(1, 1 + hist.GetNbinsX()):
                if(Unfolded.GetBinContent(k) != 0):
                    hist.SetBinContent(k, 100*sqrt(RecoilCalib.GetBinContent(k,k))/Unfolded.GetBinContent(k))
                    hist.SetBinError(k, 0)
            hist.SetName('Recoil')
            hists.append(hist)

        if  (channel.find("munu") != -1):

            hist = IsoSF.Get('MuIsoSys_Systematics_Iter1')
            hist.SetName('Iso Sf ')
            hists.append(hist)

            hist = RecoSF.Get('MuRecoSys_Systematics_Iter1')
            hist.SetName('reco Sf ')
            hists.append(hist)

            hist = TrigSF.Get('MuTrigSys_Systematics_Iter1')
            hist.SetName('Trigger Sf ')
            hists.append(hist)

            Unfolded = Summarize.Get('unfolded_data1')
            RecoilCalib = Recoil.Get('Recoil_Covariance_Iter1')
            hist = Unfolded.Clone("hist")
            for k in range(1, 1 + hist.GetNbinsX()):
                if(Unfolded.GetBinContent(k) != 0):
                    hist.SetBinContent(k, 100*sqrt(RecoilCalib.GetBinContent(k,k))/Unfolded.GetBinContent(k))
                    hist.SetBinError(k, 0)
            hist.SetName('Recoil')
            hists.append(hist)

        colorline = 0
        for i in range(0, len(hists)):
            hists[i].GetXaxis().SetRangeUser(25,60)
            hists[i].GetYaxis().SetRangeUser(0.,2)
            hists[i].SetStats(0)
            hists[i].SetTitle("")
            if(Colorcompteur == 5):
                Colorcompteur = 28
            hists[i].SetMarkerColor(Colorcompteur)
            hists[i].GetXaxis().SetTitle("p^{l}_{T} [GeV]")
            hists[i].GetYaxis().SetTitle( "Systematics[%]")
            colorline = i+1
            if(colorline == 5): colorline = 28
            hists[i].SetLineColor(colorline)
            hists[i].SetLineWidth(3)
            #utils.SetHistogramLine(hist, Colorcompteur, 2, 2, 1) # def SetHistogramLine(hist, color=1, width=1, style=1, alpha=1):
            Colorcompteur +=1

        # Draw
        c1 = root.TCanvas("c1", "The FillRandom example", 0, 0, 800, 600)
        for hist in hists:
            hist.Draw("same")
        astyle.ATLASLabel(0.2, 0.82, "Internal")
        utils.DrawText(0.2, 0.76, Indice)
        legend = makeLegend(hists,0.2, 0.35,0.4,0.7)
        legend.Draw("same")
        c1.Update()
        c1.Print("Output/"+channel+"/"+channel+"_Systematic_Diff.pdf")


    def CompareSystId(self, inputFile, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Indice):
        hists = []
        compteur = NumberOfIterationMinimal
        Colorcompteur = 1
        while compteur < NumberOfIterationMaximal:
            global hist
            if  (channel.find("enu") != -1):
                hist = inputFile.Get('ElIDSys_Systematics_Iter' + str(compteur))
            else :
                hist = inputFile.Get('MuIDSys_Systematics_Iter' + str(compteur))
            hist.GetXaxis().SetRangeUser(25, 60)

            hist.GetYaxis().SetRangeUser(0.,1)
            hist.SetStats(0)
            hist.SetTitle("")
            hist.SetName('Id Sf Iteration ' +str(compteur))
            hist.SetMarkerColor(Colorcompteur)
            hist.GetXaxis().SetTitle("p^{l}_{T} [GeV]")
            hist.GetYaxis().SetTitle( "Id Sf[%]")
            utils.SetHistogramLine(hist, Colorcompteur, 2, compteur, 1) # def SetHistogramLine(hist, color=1, width=1, style=1, alpha=1):
            hists.append(hist)
            compteur +=1
            Colorcompteur +=1
        # Draw
        c1 = root.TCanvas("c1", "The FillRandom example", 0, 0, 800, 600)
        for hist in hists:
            hist.Draw("same")
        astyle.ATLASLabel(0.2, 0.82, "Internal")
        utils.DrawText(0.2, 0.76, Indice)
        legend = makeLegend(hists,0.2, 0.35,0.4,0.7)
        legend.Draw("same")
        c1.Update()
        c1.Print("Output/"+channel+"/"+channel+"_IdSF_from_Iter_"+str(NumberOfIterationMinimal)+"_to_"+str(NumberOfIterationMaximal)+".pdf")



    def CompareSystId(self, inputFile, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Indice):
        hists = []
        compteur = NumberOfIterationMinimal
        Colorcompteur = 1
        while compteur < NumberOfIterationMaximal:
	    global hist
	    if  (channel.find("enu") != -1):
            	hist = inputFile.Get('ElIDSys_Systematics_Iter' + str(compteur))
	    else :
                hist = inputFile.Get('MuIDSys_Systematics_Iter' + str(compteur))
            hist.GetXaxis().SetRangeUser(25, 60)

            hist.GetYaxis().SetRangeUser(0.,1)
            hist.SetStats(0)
            hist.SetTitle("")
            hist.SetName('Id Sf Iteration ' +str(compteur))
            hist.SetMarkerColor(Colorcompteur)
            hist.GetXaxis().SetTitle("p^{l}_{T} [GeV]")
            hist.GetYaxis().SetTitle( "Id Sf[%]")
            utils.SetHistogramLine(hist, Colorcompteur, 2, compteur, 1) # def SetHistogramLine(hist, color=1, width=1, style=1, alpha=1):
            hists.append(hist)
            compteur +=1
            Colorcompteur +=1
        # Draw
        c1 = root.TCanvas("c1", "The FillRandom example", 0, 0, 800, 600)
        for hist in hists:
            hist.Draw("same")
        astyle.ATLASLabel(0.2, 0.82, "Internal")
        utils.DrawText(0.2, 0.76, Indice)
        legend = makeLegend(hists,0.2, 0.35,0.4,0.7)
        legend.Draw("same")
        c1.Update()
        c1.Print("Output/"+channel+"/"+channel+"_IdSF_from_Iter_"+str(NumberOfIterationMinimal)+"_to_"+str(NumberOfIterationMaximal)+".pdf")


    def CompareSystIso(self, inputFile, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Indice):
        hists = []
        compteur = NumberOfIterationMinimal
        Colorcompteur = 1
        while compteur < NumberOfIterationMaximal:
            if  (channel.find("enu") != -1):
		print("Electron")
            	hist = inputFile.Get('ElIsoSys_Systematics_Iter' + str(compteur))
            if  (channel.find("munu") != -1):
		print("muon")
                hist = inputFile.Get('MuIsoSys_Systematics_Iter' + str(compteur))
            hist.GetXaxis().SetRangeUser(25, 60)
            hist.GetYaxis().SetRangeUser(0,0.2)
            hist.SetStats(0)
            hist.SetTitle("")
            hist.SetName('Iso Sf Iteration ' +str(compteur))
            hist.SetMarkerColor(Colorcompteur)
            hist.GetXaxis().SetTitle("p^{l}_{T} [GeV]")
            hist.GetYaxis().SetTitle( "Iso Sf[%]")
            utils.SetHistogramLine(hist, Colorcompteur, 2, compteur, 1) # def SetHistogramLine(hist, color=1, width=1, style=1, alpha=1):
            hists.append(hist)
            compteur +=1
            Colorcompteur +=1
        # Draw
        c1 = root.TCanvas("c1", "The FillRandom example", 0, 0, 800, 600)
        for hist in hists:
            hist.Draw("same")
        astyle.ATLASLabel(0.2, 0.82, "Internal")
        utils.DrawText(0.2, 0.76, Indice)
        legend = makeLegend(hists,0.2, 0.35,0.4,0.7)
        legend.Draw("same")
        c1.Update()
        c1.Print("Output/"+channel+"/"+channel+"_IsoSF_from_Iter_"+str(NumberOfIterationMinimal)+"_to_"+str(NumberOfIterationMaximal)+".pdf")



    def CompareSystReco(self, inputFile, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Indice):
        hists = []
        compteur = NumberOfIterationMinimal
        Colorcompteur = 1
        while compteur < NumberOfIterationMaximal:
            if  (channel.find("enu") != -1):
            	hist = inputFile.Get('ElRecoSys_Systematics_Iter' + str(compteur))
            if  (channel.find("munu") != -1):
                hist = inputFile.Get('MuRecoSys_Systematics_Iter' + str(compteur))
            hist.GetXaxis().SetRangeUser(25, 60)
            hist.GetYaxis().SetRangeUser(0.,1)
            hist.SetStats(0)
            hist.SetTitle("")
            hist.SetName('Reco Sf Iteration ' +str(compteur))
            hist.SetMarkerColor(Colorcompteur)
            hist.GetXaxis().SetTitle("p^{l}_{T} [GeV]")
            hist.GetYaxis().SetTitle( "Reco SF[%]")
            utils.SetHistogramLine(hist, Colorcompteur, 2, compteur, 1) # def SetHistogramLine(hist, color=1, width=1, style=1, alpha=1):
            hists.append(hist)
            compteur +=1
            Colorcompteur +=1
        # Draw
        c1 = root.TCanvas("c1", "The FillRandom example", 0, 0, 800, 600)
        for hist in hists:
            hist.Draw("same")
        astyle.ATLASLabel(0.2, 0.82, "Internal")
        utils.DrawText(0.2, 0.76, Indice)
        legend = makeLegend(hists,0.2, 0.35,0.4,0.7)
        legend.Draw("same")
        c1.Update()
        c1.Print("Output/"+channel+"/"+channel+"_RecoSF_from_Iter_"+str(NumberOfIterationMinimal)+"_to_"+str(NumberOfIterationMaximal)+".pdf")

    def CompareSystTrig(self, inputFile, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Indice):
        hists = []
        compteur = NumberOfIterationMinimal
        Colorcompteur = 1
        while compteur < NumberOfIterationMaximal:
            if  (channel.find("enu")  != -1):
            	hist = inputFile.Get('ElTrigSys_Systematics_Iter' + str(compteur))
            if  (channel.find("munu") != -1):
                hist = inputFile.Get('MuTrigSys_Systematics_Iter' + str(compteur))
            hist.GetXaxis().SetRangeUser(25, 60)
            hist.GetYaxis().SetRangeUser(0.,1)
            hist.SetStats(0)
            hist.SetTitle("")
            hist.SetName('Trig Sf Iteration ' +str(compteur))
            hist.SetMarkerColor(Colorcompteur)
            hist.GetXaxis().SetTitle("p^{l}_{T} [GeV]")
            hist.GetYaxis().SetTitle( "Trig Sf[%]")
            utils.SetHistogramLine(hist, Colorcompteur, 2, compteur, 1) # def SetHistogramLine(hist, color=1, width=1, style=1, alpha=1):
            hists.append(hist)
            compteur +=1
            Colorcompteur +=1
        # Draw
        c1 = root.TCanvas("c1", "The FillRandom example", 0, 0, 800, 600)
        for hist in hists:
            hist.Draw("same")
        astyle.ATLASLabel(0.2, 0.82, "Internal")
        utils.DrawText(0.2, 0.76, Indice)
        legend = makeLegend(hists,0.2, 0.35,0.4,0.7)
        legend.Draw("same")
        c1.Update()
        c1.Print("Output/"+channel+"/"+channel+"_TrigSF_from_Iter_"+str(NumberOfIterationMinimal)+"_to_"+str(NumberOfIterationMaximal)+".pdf")



    def CompareSystRecoil(self, inputFile, RecoilSyst, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Indice):
        hists = []
        compteur = NumberOfIterationMinimal
        Colorcompteur = 1
        while compteur < NumberOfIterationMaximal:
	    Unfolded   = inputFile.Get('unfolded_data'+str(compteur))
            Covariance = RecoilSyst.Get('Recoil_Covariance_Iter' + str(compteur))

 	    hist = Unfolded.Clone("hist")
	    for k in range(1, 1 + hist.GetNbinsX()):
		if (Unfolded.GetBinContent(k) != 0):
	           hist.SetBinContent(k, 100*sqrt(Covariance.GetBinContent(k,k))/Unfolded.GetBinContent(k))
	           hist.SetBinError(k, 0)

            hist.GetXaxis().SetRangeUser(25, 60)
            hist.GetYaxis().SetRangeUser(0,1)
            hist.SetStats(0)
            hist.SetTitle("")
            hist.SetName(' Recoil Iteration ' +str(compteur))
            hist.SetMarkerColor(Colorcompteur)
            hist.GetXaxis().SetTitle("p^{l}_{T} [GeV]")
            hist.GetYaxis().SetTitle( "Recoil[%]")
            utils.SetHistogramLine(hist, Colorcompteur, 2, compteur, 1) # def SetHistogramLine(hist, color=1, width=1, style=1, alpha=1):
            hists.append(hist)
            compteur +=1
            Colorcompteur +=1
        # Draw
        c1 = root.TCanvas("c1", "The FillRandom example", 0, 0, 800, 600)
        for hist in hists:
            hist.Draw("same")
        astyle.ATLASLabel(0.2, 0.82, "Internal")
        utils.DrawText(0.2, 0.76, Indice)
        legend = makeLegend(hists,0.2, 0.25,0.4,0.7)
        legend.Draw("same")
        c1.Update()
        c1.Print("Output/"+channel+"/"+channel+"_Recoil_from_Iter_"+str(NumberOfIterationMinimal)+"_to_"+str(NumberOfIterationMaximal)+".pdf")

    def CompareSystCalib(self, inputFile, CalibSyst, NumberOfIterationMinimal, NumberOfIterationMaximal, channel, Indice):
        hists = []
        compteur = NumberOfIterationMinimal
        Colorcompteur = 1
        while compteur < NumberOfIterationMaximal:

            Unfolded   = inputFile.Get('unfolded_data'+str(compteur))
            Covariance = CalibSyst.Get('Calib_Covariance_Iter' + str(compteur))

            hist = Unfolded.Clone("hist")
            for k in range(1, 1 + hist.GetNbinsX()):
		if( Unfolded.GetBinContent(k) != 0):
                	hist.SetBinContent(k, 100*sqrt(Covariance.GetBinContent(k,k))/Unfolded.GetBinContent(k))
                	hist.SetBinError(k, 0)

            hist.GetXaxis().SetRangeUser(25, 60)
            hist.GetYaxis().SetRangeUser(0,1)
            hist.SetStats(0)
            hist.SetTitle("")
            hist.SetName('Calib Iteration ' +str(compteur))
            hist.SetMarkerColor(Colorcompteur)
            hist.GetXaxis().SetTitle("p^{l}_{T} [GeV]")
            hist.GetYaxis().SetTitle( "Calib[%]")
            utils.SetHistogramLine(hist, Colorcompteur, 2, compteur, 1) # def SetHistogramLine(hist, color=1, width=1, style=1, alpha=1):
            hists.append(hist)
            compteur +=1
            Colorcompteur +=1
        # Draw
        c1 = root.TCanvas("c1", "The FillRandom example", 0, 0, 800, 600)
        for hist in hists:
            hist.Draw("same")
        astyle.ATLASLabel(0.2, 0.82, "Internal")
        utils.DrawText(0.2, 0.76, Indice)
        legend = makeLegend(hists,0.2, 0.25,0.4,0.7)
        legend.Draw("same")
        c1.Update()
        c1.Print("Output/"+channel+"/"+channel+"_Calib_from_Iter_"+str(NumberOfIterationMinimal)+"_to_"+str(NumberOfIterationMaximal)+".pdf")
