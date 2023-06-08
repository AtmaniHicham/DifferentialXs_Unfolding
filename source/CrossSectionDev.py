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
import numpy as np
from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TAttFill, TMatrixD, THilbertMatrixD, TDecompSVD

def makeLegend(hists, xmin, ymin, xmax, ymax):
    legend = root.TLegend(xmin, ymin, xmax, ymax)
    legend.SetTextSize(0.03)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetBorderSize(0)
    for hist in hists:
        legend.AddEntry(hist, hist.GetName())
    return legend


class CrossSectionDev:
    """Classe représentant une personne"""

    def __init__(self):
        """Constructeur de notre classe"""

    def CalculateFidXsElectrons(self, Summarize_minusenu5, Summarize_plusenu5, Niter, Lum, Energy, SystLum, MCsamples_minusenu5, MCsamples_plusenu5, SystTocal_minusenu5, SystTocal_plusenu5, WminusenuByB, WplusenuByB):

        # **********************************************************************************************************************************************
        # **************************************************** calculate the Fiducial cross-Section ****************************************************
        # **********************************************************************************************************************************************

        Reco_MC_minusenu        = Summarize_minusenu5.Get("Truth_MC_Raw")
        Reco_MC_plusenu         = Summarize_plusenu5.Get("Truth_MC_Raw")

        hreco_noFakes_minusenu  = Summarize_minusenu5.Get("htruth_noMiss")
        hreco_noFakes_plusenu   = Summarize_plusenu5.Get("htruth_noMiss")

        Unfolded1               = Summarize_minusenu5.Get("Unfolded_data_iteration_"  + str(Niter))
        Unfolded3               = Summarize_plusenu5.Get("Unfolded_data_iteration_"   + str(Niter))

        sum1f = Unfolded1.Integral() * ( Reco_MC_minusenu.Integral() / hreco_noFakes_minusenu.Integral() )
        sum3f = Unfolded3.Integral() * ( Reco_MC_plusenu.Integral()  / hreco_noFakes_plusenu.Integral() )

        # bin by bin pT_Unfolding
        Recominusenu               =  Summarize_minusenu5.Get("Reco_MC")
        Truthminusenu              =  Summarize_minusenu5.Get("Truth_MC_Raw")
        data_BgSubtminusenu        =  Summarize_minusenu5.Get("data_BgSubtract")
        TruthSelectionHminusenu    =  MCsamples_minusenu5.Get("WminusenuSelectionCutFlow")
        AcceptanceFactorminusenu   =  TruthSelectionHminusenu.GetBinContent(5) / TruthSelectionHminusenu.GetBinContent(1)
        secminusenu                = (data_BgSubtminusenu.Integral() / Lum)*(Truthminusenu.Integral() / Recominusenu.Integral())

        Recoplusenu               =  Summarize_plusenu5.Get("Reco_MC")
        Truthplusenu              =  Summarize_plusenu5.Get("Truth_MC_Raw")
        data_BgSubtplusenu        =  Summarize_plusenu5.Get("data_BgSubtract")
        TruthSelectionHplusenu    =  MCsamples_plusenu5.Get("WplusenuSelectionCutFlow")
        AcceptanceFactorplusenu   =  TruthSelectionHplusenu.GetBinContent(5) / TruthSelectionHplusenu.GetBinContent(1)
        secplusenu                =  (data_BgSubtplusenu.Integral() / Lum)*(Truthplusenu.Integral() / Recoplusenu.Integral())

        print("Nombre of data     +:  ", data_BgSubtplusenu.Integral())
        print("Rapport Truth/Reco +:  ", (Truthplusenu.Integral() / Recoplusenu.Integral()))
        print("Xs Byb             +:  ", (data_BgSubtplusenu.Integral()/Lum)*(Truthplusenu.Integral() / Recoplusenu.Integral()))
        print("Xs Unfolding       +:  ", sum3f/Lum)


        """
        # **********************************************************************************************************************************************
        # ************************************************ calculate the cross-Section Uncertainties ***************************************************
        # **********************************************************************************************************************************************


        # Define the number of Events
        Nevents_Wminusenu  = 0
        Nevents_Wplusenu   = 0

        i=0
        while i < Summarize_minusenu5.Get("Unfolded_data_iteration_"+str(Niter)).GetNbinsX():
            Nevents_Wminusenu  =  Nevents_Wminusenu  +  Summarize_minusenu5.Get("Unfolded_data_iteration_"+str(Niter)).GetBinContent(i+1)
            Nevents_Wplusenu   =  Nevents_Wplusenu   +  Summarize_plusenu5.Get("Unfolded_data_iteration_"+str(Niter)).GetBinContent(i+1)
            i=i+1

        # Define the Systematics
        Bias_minusenu             = 0
        Bias_plusenu              = 0
        StatError_minusenu        = 0
        StatError_plusenu         = 0
        SystTotal_minusenu        = 0
        SystTotal_plusenu         = 0

        for i in range(0, SystTocal_minusenu5[1].GetNbinsX() ):
            for j in range(0, SystTocal_minusenu5[1].GetNbinsX()):

                StatError_minusenu  =  StatError_minusenu  +  SystTocal_minusenu5[1].GetBinContent(i+1,j+1)
                Bias_minusenu       =  Bias_minusenu       +  SystTocal_minusenu5[2].GetBinContent(i+1,j+1)

                StatError_plusenu   =  StatError_plusenu   +  SystTocal_plusenu5[1].GetBinContent(i+1,j+1)
                Bias_plusenu        =  Bias_plusenu        +  SystTocal_plusenu5[2].GetBinContent(i+1,j+1)


        SystTotal_minusenu = sqrt( SystTocal_minusenu5[0].GetBinContent(1)*SystTocal_minusenu5[0].GetBinContent(1) + SystTocal_minusenu5[0].GetBinContent(2)*SystTocal_minusenu5[0].GetBinContent(2) + SystTocal_minusenu5[0].GetBinContent(3)*SystTocal_minusenu5[0].GetBinContent(3) + SystTocal_minusenu5[0].GetBinContent(4)*SystTocal_minusenu5[0].GetBinContent(4) + SystTocal_minusenu5[0].GetBinContent(5)*SystTocal_minusenu5[0].GetBinContent(5) + SystTocal_minusenu5[0].GetBinContent(6)*SystTocal_minusenu5[0].GetBinContent(6)  )

        SystTotal_plusenu  = sqrt( SystTocal_plusenu5[0].GetBinContent(1)*SystTocal_plusenu5[0].GetBinContent(1) + SystTocal_plusenu5[0].GetBinContent(2)*SystTocal_plusenu5[0].GetBinContent(2) + SystTocal_plusenu5[0].GetBinContent(3)*SystTocal_plusenu5[0].GetBinContent(3) + SystTocal_plusenu5[0].GetBinContent(4)*SystTocal_plusenu5[0].GetBinContent(4) + SystTocal_plusenu5[0].GetBinContent(5)*SystTocal_plusenu5[0].GetBinContent(5) + SystTocal_plusenu5[0].GetBinContent(6)*SystTocal_plusenu5[0].GetBinContent(6)   )

        StatError_minusenu = 100*sqrt(StatError_minusenu)/Nevents_Wminusenu
        StatError_plusenu  = 100*sqrt(StatError_plusenu)/Nevents_Wplusenu

        Bias_minusenu = (100*Bias_minusenu/Nevents_Wminusenu)
        Bias_plusenu = (100*Bias_plusenu/Nevents_Wplusenu)

        print("Id SF    " , SystTocal_plusenu5[0].GetBinContent(1)  )
        print("Reco SF  " , SystTocal_plusenu5[0].GetBinContent(2)  )
        print("Iso SF   " , SystTocal_plusenu5[0].GetBinContent(3)  )
        print("Trig SF  " , SystTocal_plusenu5[0].GetBinContent(4)  )
        print("Recoil   " , SystTocal_plusenu5[0].GetBinContent(5)  )
        print("Calib    " , SystTocal_plusenu5[0].GetBinContent(6)  )
        print("Bias     " , Bias_plusenu      )
        print("Stat     " , StatError_plusenu )
        print("Syst     " , SystTotal_plusenu )


        #latexFile = open("Output/LatexTableau/latex5TeV.tex","w+")
        latexFile = open("Output/LatexTableau/FiducialCross_Section_Elec_"+Energy+".tex","w+")
        latexFile.write("\\documentclass[12pt]{article} \n")
        latexFile.write("\\usepackage{amsmath}\n")
        latexFile.write("\\usepackage{graphicx}\n")
        latexFile.write("\\usepackage{hyperref}\n")
        latexFile.write("\\usepackage{hyperref}\n")
        latexFile.write("\\usepackage[latin1]{inputenc}\n")
        latexFile.write("\\begin{document}\n")

        latexFile.write("\\begin{table}[ht]\n")
        latexFile.write("\\begin{tabular}{l|l|}\n")
        latexFile.write("\\cline{2-2}\n")
        latexFile.write("                                                           &    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ W$^{-}$ $\\rightarrow$ e$^{-} \\nu $  \\\ \\hline \n")
        latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ (Unfolding)}         &    %5.1f   $\\pm$ %5.3f(Stat) $\\pm$ %5.3f(Syst) $\\pm$ %5.3f(Unf) $\\pm$ %5.f(Lum)     \\\ \\hline \n" %(sum1f/Lum,    StatError_minusenu,  SystTotal_minusenu, Bias_minusenu, 5,))
        latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ $(bin\\_by\\_bin)$}  &    %5.1f   $\\pm$ %5.3f(Stat) $\\pm$ %5.3f(Syst) $\\pm$ %5.3f(Unf) $\\pm$ %5.f(Lum)     \\\ \\hline \n" %(secminusenu,  WminusenuByB[0],     WminusenuByB[1],     00,        5,))
        latexFile.write("\\end{tabular}\n")
        latexFile.write("\\end{table}\n")


        latexFile.write("\\begin{table}[ht]\n")
        latexFile.write("\\begin{tabular}{l|l|}\n")
        latexFile.write("\\cline{2-2}\n")
        latexFile.write("                                                           &    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ W$^{+}$ $\\rightarrow$ e$^{+} \\nu $  \\\ \\hline \n")
        latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ (Unfolding)}         &    %5.1f   $\\pm$ %5.3f(Stat) $\\pm$ %5.3f(Syst) $\\pm$ %5.3f(Unf) $\\pm$ %5.f(Lum)     \\\ \\hline \n" %(sum3f/Lum,   StatError_plusenu, SystTotal_plusenu, Bias_plusenu, 5,))
        latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ $(bin\\_by\\_bin)$}  &    %5.1f   $\\pm$ %5.3f(Stat) $\\pm$ %5.3f(Syst) $\\pm$ %5.3f(Unf) $\\pm$ %5.f(Lum)     \\\ \\hline \n" %(secplusenu,  WplusenuByB[0],    WplusenuByB[1],    00,        5,))
        latexFile.write("\\end{tabular}\n")
        latexFile.write("\\end{table}\n")


        latexFile.write("\\end{document}\n")
        latexFile.close()
        """

    def CalculateFidXsMuons(self, Summarize_minusenu5, Summarize_plusenu5, Niter, Lum, Energy, SystLum, MCsamples_minusenu5, MCsamples_plusenu5, SystTocal_minusenu5, SystTocal_plusenu5, WminusenuByB, WplusenuByB):

        # **********************************************************************************************************************************************
        # **************************************************** calculate the Fiducial cross-Section ****************************************************
        # **********************************************************************************************************************************************

        Reco_MC_minusenu       = Summarize_minusenu5.Get("Truth_MC_Raw")
        hreco_noFakes_minusenu = Summarize_minusenu5.Get("htruth_noMiss")
        hreco_noFakes_minusenu.Divide(Reco_MC_minusenu)

        Reco_MC_plusenu         = Summarize_plusenu5.Get("Truth_MC_Raw")
        hreco_noFakes_plusenu   = Summarize_plusenu5.Get("htruth_noMiss")
        hreco_noFakes_plusenu.Divide(Reco_MC_plusenu)

        Unfolded1    = Summarize_minusenu5.Get("Unfolded_data_iteration_"  + str(Niter))
        Unfolded3    = Summarize_plusenu5.Get("Unfolded_data_iteration_"   + str(Niter))

        sum1 = 0
        sum3 = 0
        sum1f = 0
        sum3f = 0
        sum1Error = 0
        sum3Error = 0

        i=0
        while i < Unfolded1.GetNbinsX():
            sum1      = sum1 + Unfolded1.GetBinContent(i+1)
            sum1Error = sum1Error +Unfolded1.GetBinError(i+1)
            if ( hreco_noFakes_minusenu.GetBinContent(i+1) != 0 ):
                sum1f      = sum1f + Unfolded1.GetBinContent(i+1)*(1/hreco_noFakes_minusenu.GetBinContent(i+1))
            sum3      = sum3 + Unfolded3.GetBinContent(i+1)
            sum3Error = sum3Error +Unfolded3.GetBinError(i+1)
            if ( hreco_noFakes_plusenu.GetBinContent(i+1) != 0 ):
                    sum3f      = sum3f + Unfolded3.GetBinContent(i+1)*(1/hreco_noFakes_plusenu.GetBinContent(i+1))
            i=i+1

        # bin by bin pT_Unfolding
        Recominusenu               =  Summarize_minusenu5.Get("Reco_MC")
        Truthminusenu              =  Summarize_minusenu5.Get("Truth_MC_Raw")
        data_BgSubtminusenu        =  Summarize_minusenu5.Get("data_BgSubtract")
        TruthSelectionHminusenu    =  MCsamples_minusenu5.Get("WminusmunuSelectionCutFlow")
        AcceptanceFactorminusenu   =  TruthSelectionHminusenu.GetBinContent(5) / TruthSelectionHminusenu.GetBinContent(1)
        i=0
        secminusenu    = 0
        secminusenu    = (data_BgSubtminusenu.Integral() / Lum)*(Truthminusenu.Integral() / Recominusenu.Integral())
        print("Section efficace fiducial:", secminusenu)
        while i < Unfolded1.GetNbinsX():
            if( Recominusenu.GetBinContent(i+1) != 0):
                secminusenu         = secminusenu #+ (data_BgSubtminusenu.GetBinContent(i+1) / Lum)*( Truthminusenu.GetBinContent(i+1) / Recominusenu.GetBinContent(i+1))
            i=i+1

        Recoplusenu               =  Summarize_plusenu5.Get("Reco_MC")
        Truthplusenu              =  Summarize_plusenu5.Get("Truth_MC_Raw")
        data_BgSubtplusenu        =  Summarize_plusenu5.Get("data_BgSubtract")
        TruthSelectionHplusenu    =  MCsamples_plusenu5.Get("WplusmunuSelectionCutFlow")
        AcceptanceFactorplusenu   =  TruthSelectionHplusenu.GetBinContent(5) / TruthSelectionHplusenu.GetBinContent(1)
        i=0
        secplusenu    = 0
        secplusenu    = (data_BgSubtplusenu.Integral() / Lum)*(Truthplusenu.Integral() / Recoplusenu.Integral())
        print("Section efficace fiducial:", secplusenu)
        while i < 200:
            if( Recoplusenu.GetBinContent(i+1) != 0):
                secplusenu         = secplusenu #+ (data_BgSubtplusenu.GetBinContent(i+1) / Lum)*( Truthplusenu.GetBinContent(i+1) / Recoplusenu.GetBinContent(i+1))
            i=i+1

        # **********************************************************************************************************************************************
        # ************************************************ calculate the cross-Section Uncertainties ***************************************************
        # **********************************************************************************************************************************************


        # Define the number of Events
        Nevents_Wminusenu  = 0
        Nevents_Wplusenu   = 0

        i=0
        while i < 15:
            Nevents_Wminusenu  =  Nevents_Wminusenu  +  Summarize_minusenu5.Get("Unfolded_data_iteration_"+str(Niter)).GetBinContent(i+1)
            Nevents_Wplusenu   =  Nevents_Wplusenu   +  Summarize_plusenu5.Get("Unfolded_data_iteration_"+str(Niter)).GetBinContent(i+1)
            i=i+1

        # Define the Systematics
        Bias_minusenu             = 0
        Bias_plusenu              = 0
        StatError_minusenu        = 0
        StatError_plusenu         = 0
        SystTotal_minusenu        = 0
        SystTotal_plusenu         = 0

        for i in range(0, 15):
            for j in range(0, 15):

                StatError_minusenu  =  StatError_minusenu  +  SystTocal_minusenu5[1].GetBinContent(i+1,j+1)
                Bias_minusenu       =  Bias_minusenu       +  SystTocal_minusenu5[2].GetBinContent(i+1,j+1)

                StatError_plusenu   =  StatError_plusenu   +  SystTocal_plusenu5[1].GetBinContent(i+1,j+1)
                Bias_plusenu        =  Bias_plusenu        +  SystTocal_plusenu5[2].GetBinContent(i+1,j+1)

        SystTotal_minusenu = sqrt( SystTocal_minusenu5[0].GetBinContent(1)*SystTocal_minusenu5[0].GetBinContent(1) + SystTocal_minusenu5[0].GetBinContent(2)*SystTocal_minusenu5[0].GetBinContent(2) + SystTocal_minusenu5[0].GetBinContent(3)*SystTocal_minusenu5[0].GetBinContent(3) + SystTocal_minusenu5[0].GetBinContent(4)*SystTocal_minusenu5[0].GetBinContent(4) + SystTocal_minusenu5[0].GetBinContent(5)*SystTocal_minusenu5[0].GetBinContent(5) + SystTocal_minusenu5[0].GetBinContent(6)*SystTocal_minusenu5[0].GetBinContent(6) + (100*Bias_minusenu/Nevents_Wminusenu)*(100*Bias_minusenu/Nevents_Wminusenu)  )

        SystTotal_plusenu  = sqrt( SystTocal_plusenu5[0].GetBinContent(1)*SystTocal_plusenu5[0].GetBinContent(1) + SystTocal_plusenu5[0].GetBinContent(2)*SystTocal_plusenu5[0].GetBinContent(2) + SystTocal_plusenu5[0].GetBinContent(3)*SystTocal_plusenu5[0].GetBinContent(3) + SystTocal_plusenu5[0].GetBinContent(4)*SystTocal_plusenu5[0].GetBinContent(4) + SystTocal_plusenu5[0].GetBinContent(5)*SystTocal_plusenu5[0].GetBinContent(5) + SystTocal_plusenu5[0].GetBinContent(6)*SystTocal_plusenu5[0].GetBinContent(6) + (100*Bias_plusenu/Nevents_Wplusenu)*(100*Bias_plusenu/Nevents_Wplusenu)  )

        StatError_minusenu = 100*sqrt(StatError_minusenu)/Nevents_Wminusenu
        StatError_plusenu  = 100*sqrt(StatError_plusenu)/Nevents_Wplusenu

        print("Id SF    " , SystTocal_plusenu5[0].GetBinContent(1)  )
        print("Reco SF  " , SystTocal_plusenu5[0].GetBinContent(2)  )
        print("Iso SF   " , SystTocal_plusenu5[0].GetBinContent(3)  )
        print("Trig SF  " , SystTocal_plusenu5[0].GetBinContent(4)  )
        print("Recoil   " , SystTocal_plusenu5[0].GetBinContent(5)  )
        print("Calib    " , SystTocal_plusenu5[0].GetBinContent(6)  )
        print("Bias     " , 100*Bias_plusenu/Nevents_Wplusenu      )
        print("Stat     " , StatError_plusenu )
        print("Syst     " , SystTotal_plusenu )


        #latexFile = open("Output/LatexTableau/latex5TeV.tex","w+")
        latexFile = open("Output/LatexTableau/FiducialCross_Section_Muon_"+Energy+".tex","w+")
        latexFile.write("\\documentclass[12pt]{article} \n")
        latexFile.write("\\usepackage{amsmath}\n")
        latexFile.write("\\usepackage{graphicx}\n")
        latexFile.write("\\usepackage{hyperref}\n")
        latexFile.write("\\usepackage{hyperref}\n")
        latexFile.write("\\usepackage[latin1]{inputenc}\n")
        latexFile.write("\\begin{document}\n")

        latexFile.write("\\begin{table}[ht]\n")
        latexFile.write("\\begin{tabular}{l|l|}\n")
        latexFile.write("\\cline{2-2}\n")
        latexFile.write("                                                           &    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ W$^{-}$ $\\rightarrow$ $\\mu$$^{-} \\nu $  \\\ \\hline \n")
        latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ (Unfolding)}         &    %5.1f   $\\pm$ %5.3f(Stat) $\\pm$ %5.3f(Syst) $\\pm$ %5.3f(Unf) $\\pm$ %5.f(Lum)     \\\ \\hline \n" %(sum1f/Lum,    StatError_minusenu,  SystTotal_minusenu, Bias_minusenu, 5,))
        latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ $(bin\\_by\\_bin)$}  &    %5.1f   $\\pm$ %5.3f(Stat) $\\pm$ %5.3f(Syst) $\\pm$ %5.3f(Unf) $\\pm$ %5.f(Lum)     \\\ \\hline \n" %(secminusenu,  WminusenuByB[0],     WminusenuByB[1],     00,        5,))
        latexFile.write("\\end{tabular}\n")
        latexFile.write("\\end{table}\n")


        latexFile.write("\\begin{table}[ht]\n")
        latexFile.write("\\begin{tabular}{l|l|}\n")
        latexFile.write("\\cline{2-2}\n")
        latexFile.write("                                                           &    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ W$^{+}$ $\\rightarrow$ $\\mu$$^{+} \\nu $  \\\ \\hline \n")
        latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ (Unfolding)}         &    %5.1f   $\\pm$ %5.3f(Stat) $\\pm$ %5.3f(Syst) $\\pm$ %5.3f(Unf) $\\pm$ %5.f(Lum)     \\\ \\hline \n" %(sum3f/Lum,   StatError_plusenu, SystTotal_plusenu, Bias_plusenu, 5,))
        latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ $(bin\\_by\\_bin)$}  &    %5.1f   $\\pm$ %5.3f(Stat) $\\pm$ %5.3f(Syst) $\\pm$ %5.3f(Unf) $\\pm$ %5.f(Lum)     \\\ \\hline \n" %(secplusenu,  WplusenuByB[0],    WplusenuByB[1],    00,        5,))
        latexFile.write("\\end{tabular}\n")
        latexFile.write("\\end{table}\n")

        latexFile.write("\\end{document}\n")
        latexFile.close()
