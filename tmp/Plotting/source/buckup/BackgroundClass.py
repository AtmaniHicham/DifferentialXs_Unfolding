#!/usr/bin/env python
# -*-coding:Latin-1 -*
import atlasplots
from   atlasplots import atlas_style as astyle
from   atlasplots import utils
from   atlasplots import config_reader as config

import ROOT
import ROOT as root
from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine,TChain

from array import array

def makeLegend(hists, xmin, ymin, xmax, ymax):
    legend = root.TLegend(xmin, ymin, xmax, ymax)
    legend.SetTextSize(0.03)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetBorderSize(0)
    for hist in hists:
        legend.AddEntry(hist, hist.GetName())
    return legend

def ColorParameter(Hist, FillStyle, FillColor, LineColor, LineWidth):
    Hist.SetFillStyle(FillStyle)
    Hist.SetFillColor(FillColor)
    Hist.SetLineColor(LineColor)
    Hist.SetLineWidth(LineWidth)

class BackgroundClass:
    """Classe représentant une personne"""

    def __init__(self):
        """Constructeur de notre classe"""


    def BackgroundPlotspTlepton(self, data, Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy):

        if(channel == "Wplusmunu" or channel == "Wminusmunu" ):
            lepton = "mu"
        else:
            lepton = "el"
        Nsignal                 = Signal.Get(channel+'Selection/'+lepton+'Pt_cut7')                # Signal
        Ndata                   = data.Get(channel+'Selection/'+lepton+'Pt_cut7')                  # Data

        print("Nombre de bins de MC : ",Nsignal.GetNbinsX())

        NBackgroundW            = Background_W.Get(channel+'Selection/'+lepton+'Pt_cut7')          # background "W"
        NBackgroundZ            = Background_Z.Get(channel+'Selection/'+lepton+'Pt_cut7')          # background "Z"
        NBackgroundDiboson      = Background_diboson.Get(channel+'Selection/'+lepton+'Pt_cut7')    # background "diboson"
        NBackgroundTop          = Background_Top.Get(channel+'Selection/'+lepton+'Pt_cut7')        # background "Top"
        NBackgroundMultijet     = Background_MiltiJet.Get('hist/'+lepton+'Pt')                     # background "Miltijet"

        # Make a Clone of hists
        Hdata               = TH1F("Hdata",                 "data",                 100, 0, 100)    # data
        Hsignal             = TH1F("Hsignal1",              "signal",               100, 0, 100)    # Signal
        HBackgroundW        = TH1F("NBackgroundW",          "Background_W",         100, 0, 100)    # W+-
        HBackgroundZ        = TH1F("HBackgroundZ",          "Background_Z",         100, 0, 100)    # Z+-
        HBackgroundDiboson  = TH1F("HBackgroundDiboson",    "Background_Diboson",   100, 0, 100)    # diboson
        HBackgroundMultijet = TH1F("HBackgroundMultijet",   "Background_Multijet",  100, 0, 100)    # MJ
        HBackgroundTop      = TH1F("HBackgroundTop",        "Background_Top",       100, 0, 100)    # Top
        Hratio1             = TH1F("Hratio1",               "ratio",                100, 0, 100)    # ratio

        i=1
        while i < 101:
            Hdata.SetBinContent(i,                  Ndata.GetBinContent(i))

            Hsignal.SetBinContent(i,                Nsignal.GetBinContent(i))
            Hsignal.SetBinError(i,                  Nsignal.GetBinError(i))

            HBackgroundW.SetBinContent(i,           NBackgroundW.GetBinContent(i))
            HBackgroundZ.SetBinContent(i,           NBackgroundZ.GetBinContent(i))
            HBackgroundDiboson.SetBinContent(i,     NBackgroundDiboson.GetBinContent(i))
            HBackgroundMultijet.SetBinContent(i,    NBackgroundMultijet.GetBinContent(i))
            HBackgroundTop.SetBinContent(i,         NBackgroundTop.GetBinContent(i))

            Hratio1.SetBinContent(i,                Ndata.GetBinContent(i))
            Hratio1.SetBinError(i,                  Ndata.GetBinError(i))
            i=i+1

        # MC Total
        Hsignal.SetLineColor(2)
        Hsignal.SetLineWidth(2)
        Hsignal.SetLineStyle(2)

        MCTotal = Hsignal.Clone("MCTotal")
        MCTotal.Add(HBackgroundDiboson)
        MCTotal.Add(HBackgroundTop)
        MCTotal.Add(HBackgroundMultijet)
        MCTotal.Add(HBackgroundW)
        MCTotal.Add(HBackgroundZ)
        MCTotal.SetMarkerSize(0)

        # Ratio Plot
        Nratio1 =  Hdata.Clone("Nratio1")

        print("Nombre de bins de data : ",Nratio1.GetNbinsX())
        print("Nombre de bins de MC : ",Hsignal.GetNbinsX())

        # Hist Syle 2, 4, 209, 93, 221
        ColorParameter(HBackgroundW, 1001, 2, 1, 1)
        ColorParameter(HBackgroundZ, 1001, 4, 1, 1)
        ColorParameter(HBackgroundDiboson, 1001, 209, 1, 1)
        ColorParameter(HBackgroundMultijet, 1001, 93, 1, 1)
        ColorParameter(HBackgroundTop, 1001, 53, 1, 1)

        Hdata.SetLineColor(1)
        Hdata.SetLineWidth(2)
        Hdata.SetMarkerStyle(1)
        Hdata.SetMarkerColor(1)
        Hdata.SetMarkerSize(1)

        Legend = ROOT.TLegend(0.6,0.45,0.88,0.85)
        Legend.AddEntry(Hdata,"data");
        Legend.AddEntry(MCTotal,"Signal+Bg");
        Legend.AddEntry(HBackgroundW,"W^{+-} #rightarrow l^{+-}v","f");
        Legend.AddEntry(HBackgroundZ,"Z #rightarrow ll","f");
        Legend.AddEntry(HBackgroundDiboson,"Diboson","f");
        Legend.AddEntry(HBackgroundMultijet,"Multijet","f");
        Legend.AddEntry(HBackgroundTop,"Top","f");
        Legend.SetBorderSize(0)

        # Tline
        line1 = ROOT.TLine(0,0.95,100,0.95)
        line2 = ROOT.TLine(0,1.05,100,1.05)
        line3 = ROOT.TLine(0,1.,100,1.)
        line1.SetLineStyle(2)
        line2.SetLineStyle(2)
        line3.SetLineStyle(2)

        # TStack
        BackgroundPlot = ROOT.THStack("ss","")
        BackgroundPlot.Add(HBackgroundDiboson)
        BackgroundPlot.Add(HBackgroundTop)
        BackgroundPlot.Add(HBackgroundMultijet)
        BackgroundPlot.Add(HBackgroundZ)
        BackgroundPlot.Add(HBackgroundW)
        BackgroundPlot.Add(Hsignal)

        # Draw
        astyle.SetAtlasStyle()
        c = TCanvas("c", "canvas", 800, 700)
        pad1 = TPad("pad1", "pad1", 0, 0.32, 1, 1.0)
        pad1.SetBottomMargin(0);
        pad1.Draw();
        pad1.SetLogy()
        pad1.cd();
        Hdata.SetStats(0)
        Hdata.SetName("")
        Hdata.SetTitle("")
        Hdata.GetYaxis().SetTitle("Event")
        Hdata.GetYaxis().SetTitleOffset(0.8)
        Hdata.GetYaxis().SetTitleSize(0.06)
        Hdata.GetXaxis().SetRangeUser(0,200)
        Hdata.SetMinimum(1.4)
        Hdata.SetMaximum(1000000)
        Hdata.Draw("same")
        BackgroundPlot.Draw("same")
        Legend.Draw("same")
        astyle.ATLASLabel(0.2, 0.82, "Internal")
        utils.DrawText(0.2, 0.76, Indice)

        c.Update()

        c.cd();
        pad2 = TPad("pad2", "pad2", 0, 0., 1, 0.3)
        pad2.SetTopMargin(0);
        pad2.SetBottomMargin(0.4);
        pad2.Draw();
        pad2.cd();
        Nratio1.GetXaxis().SetRangeUser(0,200)
        Nratio1.GetYaxis().SetRangeUser(0.899,1.101)
        Nratio1.GetXaxis().SetLabelSize(0.1)
        Nratio1.GetYaxis().SetLabelSize(0.1)
        Nratio1.GetXaxis().SetTitleSize(0.1)
        Nratio1.GetYaxis().SetTitleSize(0.1)
        Nratio1.GetYaxis().CenterTitle()
        Nratio1.GetYaxis().SetTitleOffset(0.5)
        Nratio1.GetYaxis().SetTitleSize(0.12)
        Nratio1.GetXaxis().SetTitleSize(0.12)
        Nratio1.GetXaxis().SetTitleOffset(1.2)
        Nratio1.SetTitle("")
        Nratio1.Divide(MCTotal)
        Nratio1.SetLineColor(1)
        Nratio1.SetMarkerStyle(10)
        Nratio1.SetMarkerSize(0.6)
        Nratio1.SetLineWidth(2)
        Nratio1.SetStats(0)
        Nratio1.GetXaxis().SetTitle("p_{T}^{lepton}[GeV]")
        Nratio1.GetYaxis().SetTitle("Data/MC")
        Nratio1.Draw("P")
        line1.Draw("same")
        line2.Draw("same")
        line3.Draw("same")
        c.Print("Output/Backgrounds/BackgroundPlot_ptLepton_"+channel+"_"+Energy+".pdf")


    def BackgroundPlotsetalepton(self, data, Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy):

        if(channel == "Wplusmunu" or channel == "Wminusmunu" ):
            lepton = "muEta"
        else:
            lepton = "elEta"

        Nsignal                 = Signal.Get(channel+'Selection/'+lepton+'_cut7')                # Signal
        Ndata                   = data.Get(channel+'Selection/'+lepton+'_cut7')                  # Data

        print("Nombre de bins de MC : ",Nsignal.GetNbinsX())

        NBackgroundW            = Background_W.Get(channel+'Selection/'+lepton+'_cut7')          # background "W"
        NBackgroundZ            = Background_Z.Get(channel+'Selection/'+lepton+'_cut7')          # background "Z"
        NBackgroundDiboson      = Background_diboson.Get(channel+'Selection/'+lepton+'_cut7')    # background "diboson"
        NBackgroundTop          = Background_Top.Get(channel+'Selection/'+lepton+'_cut7')        # background "Top"
        NBackgroundMultijet     = Background_MiltiJet.Get('hist/'+lepton+'Pt')                      # background "Miltijet"

        RecoBin = [-2.5, -2.18, -1.95, -1.74, -1.52, -1.37, -1.05, -0.84, -0.63, -0.42, -0.21, 0.0, 0.21, 0.42, 0.63, 0.84, 1.05, 1.37, 1.52, 1.74, 1.95, 2.18, 2.5]

        # Make a Clone of hists
        Hdata               = TH1F("Hdata",                 "data",                 len(RecoBin)-1, array('d',RecoBin))     # data
        Hsignal             = TH1F("Hsignal1",              "signal",               len(RecoBin)-1, array('d',RecoBin))   # Signal
        HBackgroundW        = TH1F("NBackgroundW",          "Background_W",         len(RecoBin)-1, array('d',RecoBin))    # W+-
        HBackgroundZ        = TH1F("HBackgroundZ",          "Background_Z",         len(RecoBin)-1, array('d',RecoBin))    # Z+-
        HBackgroundDiboson  = TH1F("HBackgroundDiboson",    "Background_Diboson",   len(RecoBin)-1, array('d',RecoBin))   # diboson
        HBackgroundMultijet = TH1F("HBackgroundMultijet",   "Background_Multijet",  len(RecoBin)-1, array('d',RecoBin))     # MJ
        HBackgroundTop      = TH1F("HBackgroundTop",        "Background_Top",       len(RecoBin)-1, array('d',RecoBin))     # Top
        Hratio1             = TH1F("Hratio1",               "ratio",                len(RecoBin)-1, array('d',RecoBin))     # ratio

        i=1
        while i < 101:
            Hdata.SetBinContent(i,                  Ndata.GetBinContent(i))

            Hsignal.SetBinContent(i,                Nsignal.GetBinContent(i))
            Hsignal.SetBinError(i,                  Nsignal.GetBinError(i))

            HBackgroundW.SetBinContent(i,           NBackgroundW.GetBinContent(i))
            HBackgroundZ.SetBinContent(i,           NBackgroundZ.GetBinContent(i))
            HBackgroundDiboson.SetBinContent(i,     NBackgroundDiboson.GetBinContent(i))
            HBackgroundMultijet.SetBinContent(i,    NBackgroundMultijet.GetBinContent(i))
            HBackgroundTop.SetBinContent(i,         NBackgroundTop.GetBinContent(i))

            Hratio1.SetBinContent(i,                Ndata.GetBinContent(i))
            Hratio1.SetBinError(i,                  Ndata.GetBinError(i))
            i=i+1

        # MC Total
        Hsignal.SetLineColor(2)
        Hsignal.SetLineWidth(2)
        Hsignal.SetLineStyle(2)

        MCTotal = Hsignal.Clone("MCTotal")
        MCTotal.Add(HBackgroundDiboson)
        MCTotal.Add(HBackgroundTop)
        MCTotal.Add(HBackgroundMultijet)
        MCTotal.Add(HBackgroundW)
        MCTotal.Add(HBackgroundZ)
        MCTotal.SetMarkerSize(0)

        # Ratio Plot
        Nratio1 =  Hdata.Clone("Nratio1")

        print("Nombre de bins de data : ",Nratio1.GetNbinsX())
        print("Nombre de bins de MC : ",Hsignal.GetNbinsX())

        # Hist Syle
        ColorParameter(HBackgroundW, 1001, 2, 1, 1)
        ColorParameter(HBackgroundZ, 1001, 4, 1, 1)
        ColorParameter(HBackgroundDiboson, 1001, 209, 1, 1)
        ColorParameter(HBackgroundMultijet, 1001, 93, 1, 1)
        ColorParameter(HBackgroundTop, 1001, 53, 1, 1)

        Hdata.SetLineColor(1)
        Hdata.SetLineWidth(2)
        Hdata.SetMarkerStyle(1)
        Hdata.SetMarkerColor(1)
        Hdata.SetMarkerSize(1)

        Legend = ROOT.TLegend(0.6,0.65,0.8,0.92)
        Legend.AddEntry(Hdata,"data");
        Legend.AddEntry(MCTotal,"Signal+Bg");
        Legend.AddEntry(HBackgroundW,"W^{+-} #rightarrow l^{+-}v","f");
        Legend.AddEntry(HBackgroundZ,"Z #rightarrow ll","f");
        Legend.AddEntry(HBackgroundDiboson,"Diboson","f");
        Legend.AddEntry(HBackgroundMultijet,"Multijet","f");
        Legend.AddEntry(HBackgroundTop,"Top","f");
        Legend.SetBorderSize(0)

        # Tline
        line1 = ROOT.TLine(-2.5,0.95,2.5,0.95)
        line2 = ROOT.TLine(-2.5,1.05,2.5,1.05)
        line3 = ROOT.TLine(-2.5,1.,2.5,1.)

        line1.SetLineStyle(2)
        line2.SetLineStyle(2)
        line3.SetLineStyle(2)

        # TStack
        BackgroundPlot = ROOT.THStack("ss","")
        BackgroundPlot.Add(HBackgroundDiboson)
        BackgroundPlot.Add(HBackgroundTop)
        BackgroundPlot.Add(HBackgroundMultijet)
        BackgroundPlot.Add(HBackgroundZ)
        BackgroundPlot.Add(HBackgroundW)
        BackgroundPlot.Add(Hsignal)

        # Draw
        astyle.SetAtlasStyle()
        c = TCanvas("c", "canvas", 800, 700)
        pad1 = TPad("pad1", "pad1", 0, 0.32, 1, 1.0)
        pad1.SetBottomMargin(0);
        pad1.Draw();
        pad1.SetLogy()
        pad1.cd();
        Hdata.SetStats(0)
        Hdata.SetName("")
        Hdata.SetTitle("")
        Hdata.GetYaxis().SetTitle("Event")
        Hdata.GetYaxis().SetTitleOffset(0.8)
        Hdata.GetYaxis().SetTitleSize(0.06)
        #Hdata.GetXaxis().SetRangeUser(0,200)
        Hdata.SetMinimum(1.4)
        Hdata.SetMaximum(1000000)
        Hdata.Draw("same")
        BackgroundPlot.Draw("same")
        Legend.Draw("same")
        astyle.ATLASLabel(0.2, 0.82, "Internal")
        utils.DrawText(0.2, 0.76, Indice)

        c.Update()

        c.cd();
        pad2 = TPad("pad2", "pad2", 0, 0., 1, 0.3)
        pad2.SetTopMargin(0);
        pad2.SetBottomMargin(0.4);
        pad2.Draw();
        pad2.cd();
        #Nratio1.GetXaxis().SetRangeUser(0,200)
        Nratio1.GetYaxis().SetRangeUser(0.899,1.101)
        Nratio1.GetXaxis().SetLabelSize(0.1)
        Nratio1.GetYaxis().SetLabelSize(0.1)
        Nratio1.GetXaxis().SetTitleSize(0.1)
        Nratio1.GetYaxis().SetTitleSize(0.1)
        Nratio1.GetYaxis().CenterTitle()
        Nratio1.GetYaxis().SetTitleOffset(0.5)
        Nratio1.GetYaxis().SetTitleSize(0.12)
        Nratio1.GetXaxis().SetTitleSize(0.12)
        Nratio1.GetXaxis().SetTitleOffset(1.2)
        Nratio1.SetTitle("")
        Nratio1.Divide(MCTotal)
        Nratio1.SetLineColor(1)
        Nratio1.SetMarkerStyle(10)
        Nratio1.SetMarkerSize(0.6)
        Nratio1.SetLineWidth(2)
        Nratio1.SetStats(0)
        Nratio1.GetXaxis().SetTitle("#eta^{lepton}[GeV]")
        Nratio1.GetYaxis().SetTitle("Data/MC")
        Nratio1.Draw("P")
        line1.Draw("same")
        line2.Draw("same")
        line3.Draw("same")
        c.Print("Output/Backgrounds/BackgroundPlot_eta_Lepton_"+channel+"_"+Energy+".pdf")


    def BackgroundPlotsmTW(self, data, Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy):

        Nsignal                 = Signal.Get(channel+'Selection/mT_cut7')                # Signal
        Ndata                   = data.Get(channel+'Selection/mT_cut7')                  # Data

        print("Nombre de bins de MC : ",Nsignal.GetNbinsX())

        NBackgroundW            = Background_W.Get(channel+'Selection/mT_cut7')          # background "W"
        NBackgroundZ            = Background_Z.Get(channel+'Selection/mT_cut7')          # background "Z"
        NBackgroundDiboson      = Background_diboson.Get(channel+'Selection/mT_cut7')    # background "diboson"
        NBackgroundTop          = Background_Top.Get(channel+'Selection/mT_cut7')        # background "Top"
        NBackgroundMultijet     = Background_MiltiJet.Get('hist/mT')                      # background "Miltijet"

        # Make a Clone of hists
        Hdata               = TH1F("Hdata",                 "data",                 100, 0, 200)    # data
        Hsignal             = TH1F("Hsignal1",              "signal",               100, 0, 200)   # Signal
        HBackgroundW        = TH1F("NBackgroundW",          "Background_W",         100, 0, 200)    # W+-
        HBackgroundZ        = TH1F("HBackgroundZ",          "Background_Z",         100, 0, 200)    # Z+-
        HBackgroundDiboson  = TH1F("HBackgroundDiboson",    "Background_Diboson",   100, 0, 200)  # diboson
        HBackgroundMultijet = TH1F("HBackgroundMultijet",   "Background_Multijet",  100, 0, 200)     # MJ
        HBackgroundTop      = TH1F("HBackgroundTop",        "Background_Top",       100, 0, 200)    # Top
        Hratio1             = TH1F("Hratio1",               "ratio",                100, 0, 200)     # ratio

        i=1
        while i < 101:
            Hdata.SetBinContent(i,                  Ndata.GetBinContent(i))

            Hsignal.SetBinContent(i,                Nsignal.GetBinContent(i))
            Hsignal.SetBinError(i,                  Nsignal.GetBinError(i))

            HBackgroundW.SetBinContent(i,           NBackgroundW.GetBinContent(i))
            HBackgroundZ.SetBinContent(i,           NBackgroundZ.GetBinContent(i))
            HBackgroundDiboson.SetBinContent(i,     NBackgroundDiboson.GetBinContent(i))
            HBackgroundMultijet.SetBinContent(i,    NBackgroundMultijet.GetBinContent(i))
            HBackgroundTop.SetBinContent(i,         NBackgroundTop.GetBinContent(i))

            Hratio1.SetBinContent(i,                Ndata.GetBinContent(i))
            Hratio1.SetBinError(i,                  Ndata.GetBinError(i))
            i=i+1

        # MC Total
        Hsignal.SetLineColor(2)
        Hsignal.SetLineWidth(2)
        Hsignal.SetLineStyle(2)

        MCTotal = Hsignal.Clone("MCTotal")
        MCTotal.Add(HBackgroundDiboson)
        MCTotal.Add(HBackgroundTop)
        MCTotal.Add(HBackgroundMultijet)
        MCTotal.Add(HBackgroundW)
        MCTotal.Add(HBackgroundZ)
        MCTotal.SetMarkerSize(0)

        # Ratio Plot
        Nratio1 =  Hdata.Clone("Nratio1")

        print("Nombre de bins de data : ",Nratio1.GetNbinsX())
        print("Nombre de bins de MC : ",Hsignal.GetNbinsX())

        # Hist Syle
        ColorParameter(HBackgroundW, 1001, 2, 1, 1)
        ColorParameter(HBackgroundZ, 1001, 4, 1, 1)
        ColorParameter(HBackgroundDiboson, 1001, 209, 1, 1)
        ColorParameter(HBackgroundMultijet, 1001, 93, 1, 1)
        ColorParameter(HBackgroundTop, 1001, 53, 1, 1)

        Hdata.SetLineColor(1)
        Hdata.SetLineWidth(2)
        Hdata.SetMarkerStyle(1)
        Hdata.SetMarkerColor(1)
        Hdata.SetMarkerSize(1)

        Legend = ROOT.TLegend(0.6,0.45,0.88,0.85)
        Legend.AddEntry(Hdata,"data");
        Legend.AddEntry(MCTotal,"Signal+Bg");
        Legend.AddEntry(HBackgroundW,"W^{+-} #rightarrow l^{+-}v","f");
        Legend.AddEntry(HBackgroundZ,"Z #rightarrow ll","f");
        Legend.AddEntry(HBackgroundDiboson,"Diboson","f");
        Legend.AddEntry(HBackgroundMultijet,"Multijet","f");
        Legend.AddEntry(HBackgroundTop,"Top","f");
        Legend.SetBorderSize(0)

        # Tline
        line1 = ROOT.TLine(0,0.95,200,0.95)
        line2 = ROOT.TLine(0,1.05,200,1.05)
        line3 = ROOT.TLine(0,1.,200,1.)
        line1.SetLineStyle(2)
        line2.SetLineStyle(2)
        line3.SetLineStyle(2)

        # TStack
        BackgroundPlot = ROOT.THStack("ss","")
        BackgroundPlot.Add(HBackgroundDiboson)
        BackgroundPlot.Add(HBackgroundTop)
        BackgroundPlot.Add(HBackgroundMultijet)
        BackgroundPlot.Add(HBackgroundZ)
        BackgroundPlot.Add(HBackgroundW)
        BackgroundPlot.Add(Hsignal)

        # Draw
        astyle.SetAtlasStyle()
        c = TCanvas("c", "canvas", 800, 700)
        pad1 = TPad("pad1", "pad1", 0, 0.32, 1, 1.0)
        pad1.SetBottomMargin(0);
        pad1.Draw();
        pad1.SetLogy()
        pad1.cd();
        Hdata.SetStats(0)
        Hdata.SetName("")
        Hdata.SetTitle("")
        Hdata.GetYaxis().SetTitle("Event")
        Hdata.GetYaxis().SetTitleOffset(0.8)
        Hdata.GetYaxis().SetTitleSize(0.06)
        Hdata.GetXaxis().SetRangeUser(0,200)
        Hdata.SetMinimum(1.4)
        Hdata.SetMaximum(1000000)
        Hdata.Draw("same")
        BackgroundPlot.Draw("same")
        Legend.Draw("same")
        astyle.ATLASLabel(0.2, 0.82, "Internal")
        utils.DrawText(0.2, 0.76, Indice)

        c.Update()

        c.cd();
        pad2 = TPad("pad2", "pad2", 0, 0., 1, 0.3)
        pad2.SetTopMargin(0);
        pad2.SetBottomMargin(0.4);
        pad2.Draw();
        pad2.cd();
        Nratio1.GetXaxis().SetRangeUser(0,200)
        Nratio1.GetYaxis().SetRangeUser(0.899,1.101)
        Nratio1.GetXaxis().SetLabelSize(0.1)
        Nratio1.GetYaxis().SetLabelSize(0.1)
        Nratio1.GetXaxis().SetTitleSize(0.1)
        Nratio1.GetYaxis().SetTitleSize(0.1)
        Nratio1.GetYaxis().CenterTitle()
        Nratio1.GetYaxis().SetTitleOffset(0.5)
        Nratio1.GetYaxis().SetTitleSize(0.12)
        Nratio1.GetXaxis().SetTitleSize(0.12)
        Nratio1.GetXaxis().SetTitleOffset(1.2)
        Nratio1.SetTitle("")
        Nratio1.Divide(MCTotal)
        Nratio1.SetLineColor(1)
        Nratio1.SetMarkerStyle(10)
        Nratio1.SetMarkerSize(0.6)
        Nratio1.SetLineWidth(2)
        Nratio1.SetStats(0)
        Nratio1.GetXaxis().SetTitle("m_{T}^{W}[GeV]")
        Nratio1.GetYaxis().SetTitle("Data/MC")
        Nratio1.Draw("P")
        line1.Draw("same")
        line2.Draw("same")
        line3.Draw("same")
        c.Print("Output/Backgrounds/BackgroundPlot_mTW_"+channel+"_"+Energy+".pdf")

        i=1
        while i < 101:
            print("ratio du bin %d  :  %f"%(i, Nratio1.GetBinContent(i)))
            i+=1


    def BackgroundPlotspTW(self, data, Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, Energy):

        Nsignal                 = Signal.Get(channel+'Selection/WpT_Reco_cut7')                # Signal
        Ndata                   = data.Get(channel+'Selection/WpT_Reco_cut7')                  # Data

        print("Nombre de bins de MC : ",Nsignal.GetNbinsX())

        NBackgroundW            = Background_W.Get(channel+'Selection/WpT_Reco_cut7')          # background "W"
        NBackgroundZ            = Background_Z.Get(channel+'Selection/WpT_Reco_cut7')          # background "Z"
        NBackgroundDiboson      = Background_diboson.Get(channel+'Selection/WpT_Reco_cut7')    # background "diboson"
        NBackgroundTop          = Background_Top.Get(channel+'Selection/WpT_Reco_cut7')        # background "Top"
        NBackgroundMultijet     = Background_MiltiJet.Get('hist/WpT_Reco_cut7')                      # background "Miltijet"

        RecoBin = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,105,110,115,120,125,130,135,140,145,150,155,160, 165,170,175,180,185,190,195,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550,560,570,580,590,600]

        # Make a Clone of hists
        Hdata               = TH1F("Hdata",                 "data",                 len(RecoBin)-1, array('d',RecoBin))     # data
        Hsignal             = TH1F("Hsignal1",              "signal",               len(RecoBin)-1, array('d',RecoBin))   # Signal
        HBackgroundW        = TH1F("NBackgroundW",          "Background_W",         len(RecoBin)-1, array('d',RecoBin))    # W+-
        HBackgroundZ        = TH1F("HBackgroundZ",          "Background_Z",         len(RecoBin)-1, array('d',RecoBin))    # Z+-
        HBackgroundDiboson  = TH1F("HBackgroundDiboson",    "Background_Diboson",   len(RecoBin)-1, array('d',RecoBin))   # diboson
        HBackgroundMultijet = TH1F("HBackgroundMultijet",   "Background_Multijet",  len(RecoBin)-1, array('d',RecoBin))     # MJ
        HBackgroundTop      = TH1F("HBackgroundTop",        "Background_Top",       len(RecoBin)-1, array('d',RecoBin))     # Top
        Hratio1             = TH1F("Hratio1",               "ratio",                len(RecoBin)-1, array('d',RecoBin))     # ratio

        i=1
        while i < 101:
            Hdata.SetBinContent(i,                  Ndata.GetBinContent(i))

            Hsignal.SetBinContent(i,                Nsignal.GetBinContent(i))
            Hsignal.SetBinError(i,                  Nsignal.GetBinError(i))

            HBackgroundW.SetBinContent(i,           NBackgroundW.GetBinContent(i))
            HBackgroundZ.SetBinContent(i,           NBackgroundZ.GetBinContent(i))
            HBackgroundDiboson.SetBinContent(i,     NBackgroundDiboson.GetBinContent(i))
            HBackgroundMultijet.SetBinContent(i,    NBackgroundMultijet.GetBinContent(i))
            HBackgroundTop.SetBinContent(i,         NBackgroundTop.GetBinContent(i))

            Hratio1.SetBinContent(i,                Ndata.GetBinContent(i))
            Hratio1.SetBinError(i,                  Ndata.GetBinError(i))
            i=i+1

        # MC Total
        Hsignal.SetLineColor(2)
        Hsignal.SetLineWidth(2)
        Hsignal.SetLineStyle(2)

        MCTotal = Hsignal.Clone("MCTotal")
        MCTotal.Add(HBackgroundDiboson)
        MCTotal.Add(HBackgroundTop)
        MCTotal.Add(HBackgroundMultijet)
        MCTotal.Add(HBackgroundW)
        MCTotal.Add(HBackgroundZ)
        MCTotal.SetMarkerSize(0)


        # Ratio Plot
        Nratio1 =  Hdata.Clone("Nratio1")

        print("Nombre de bins de data : ",Nratio1.GetNbinsX())
        print("Nombre de bins de MC : ",Hsignal.GetNbinsX())

        # Hist Syle
        ColorParameter(HBackgroundW, 1001, 2, 1, 1)
        ColorParameter(HBackgroundZ, 1001, 4, 1, 1)
        ColorParameter(HBackgroundDiboson, 1001, 209, 1, 1)
        ColorParameter(HBackgroundMultijet, 1001, 93, 1, 1)
        ColorParameter(HBackgroundTop, 1001, 53, 1, 1)

        Hdata.SetLineColor(1)
        Hdata.SetLineWidth(2)
        Hdata.SetMarkerStyle(1)
        Hdata.SetMarkerColor(1)
        Hdata.SetMarkerSize(1)

        Legend = ROOT.TLegend(0.6,0.45,0.88,0.85)
        Legend.AddEntry(Hdata,"data");
        Legend.AddEntry(MCTotal,"Signal+Bg");
        Legend.AddEntry(HBackgroundW,"W^{+-} #rightarrow l^{+-}v","f");
        Legend.AddEntry(HBackgroundZ,"Z #rightarrow ll","f");
        Legend.AddEntry(HBackgroundDiboson,"Diboson","f");
        Legend.AddEntry(HBackgroundMultijet,"Multijet","f");
        Legend.AddEntry(HBackgroundTop,"Top","f");
        Legend.SetBorderSize(0)

        # Tline
        line1 = ROOT.TLine(0,0.95,100,0.95)
        line2 = ROOT.TLine(0,1.05,100,1.05)
        line3 = ROOT.TLine(0,1.,100,1.)
        line1.SetLineStyle(2)
        line2.SetLineStyle(2)
        line3.SetLineStyle(2)

        # TStack
        BackgroundPlot = ROOT.THStack("ss","")
        BackgroundPlot.Add(HBackgroundDiboson)
        BackgroundPlot.Add(HBackgroundTop)
        BackgroundPlot.Add(HBackgroundMultijet)
        BackgroundPlot.Add(HBackgroundZ)
        BackgroundPlot.Add(HBackgroundW)
        BackgroundPlot.Add(Hsignal)

        # Draw
        astyle.SetAtlasStyle()
        c = TCanvas("c", "canvas", 800, 700)
        pad1 = TPad("pad1", "pad1", 0, 0.32, 1, 1.0)
        pad1.SetBottomMargin(0);
        pad1.Draw();
        pad1.SetLogy()
        pad1.cd();
        Hdata.SetStats(0)
        Hdata.SetName("")
        Hdata.SetTitle("")
        Hdata.GetYaxis().SetTitle("Event")
        Hdata.GetYaxis().SetTitleOffset(0.8)
        Hdata.GetYaxis().SetTitleSize(0.06)

        Hdata.GetXaxis().SetRangeUser(0,100)
        Hdata.SetMinimum(1.4)
        Hdata.SetMaximum(1000000)
        Hdata.Draw("same")
        BackgroundPlot.Draw("same")
        Legend.Draw("same")
        astyle.ATLASLabel(0.2, 0.82, "Internal")
        utils.DrawText(0.2, 0.76, Indice)

        c.Update()

        c.cd();
        pad2 = TPad("pad2", "pad2", 0, 0., 1, 0.3)
        pad2.SetTopMargin(0);
        pad2.SetBottomMargin(0.4);
        pad2.Draw();
        pad2.cd();
        Nratio1.GetXaxis().SetRangeUser(0,100)
        Nratio1.GetYaxis().SetRangeUser(0.899,1.101)
        Nratio1.GetXaxis().SetLabelSize(0.1)
        Nratio1.GetYaxis().SetLabelSize(0.1)
        Nratio1.GetXaxis().SetTitleSize(0.1)
        Nratio1.GetYaxis().SetTitleSize(0.1)
        Nratio1.GetYaxis().CenterTitle()
        Nratio1.GetYaxis().SetTitleOffset(0.5)
        Nratio1.GetYaxis().SetTitleSize(0.12)
        Nratio1.GetXaxis().SetTitleSize(0.12)
        Nratio1.GetXaxis().SetTitleOffset(1.2)
        Nratio1.SetTitle("")
        Nratio1.Divide(MCTotal)
        Nratio1.SetLineColor(1)
        Nratio1.SetLineWidth(2)
        Nratio1.SetStats(0)
        Nratio1.SetMarkerStyle(10)
        Nratio1.SetMarkerSize(0.6)
        Nratio1.GetXaxis().SetTitle("p_{T}^{W}[GeV]")
        Nratio1.GetYaxis().SetTitle("Data/MC")
        Nratio1.Draw()
        line1.Draw("same")
        line2.Draw("same")
        line3.Draw("same")
        c.Print("Output/Backgrounds/BackgroundPlot_pTW_"+channel+"_"+Energy+".pdf")
