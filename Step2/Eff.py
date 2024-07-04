import ROOT,math,os,sys
ROOT.gROOT.LoadMacro('./RooCMSShape.cc+')
ROOT.gROOT.LoadMacro('./RooCBExGaussShape.cc+')

from math import sqrt
from ROOT import RooCMSShape, RooCBExGaussShape, TCanvas, TPad
import CMSTDRStyle
CMSTDRStyle.setTDRStyle().cd()
import CMSstyle
from array import array

# shape chose could be:
# for MC (DY): nominal--use gauss to smear the MC HIST
#         altsig--use CB+Gauss to smear the MC HIST

# for DATA fit, 
# the signal shape: nominal--use gauss to smear the MC HIST
#                   altsig--use CB+Gauss to smear the MC HIST
#for bkg shape: nominal--use RooCMSShape
#               altbkg--use exponential function
# and for DATA fit, the signal shape should be constrained from the MC DY fit

ismc=sys.argv[1]
sig_shape=sys.argv[2] # gauss or cbexgauss 
bkg_shape=sys.argv[3] # cmsshape or expo

# unc of MC: NLO vs LO, pileup Up vs Down, so the choices are: 1. nominal, 2. LO, 3. pileup
unc=sys.argv[4] # nominal, LO, puUp, puDo


def Eff(filename,ismc,sig_shape,bkg_shape,unc):

  path='./'
  if unc=='LO':
    path='../Inputs/LO/'
    path_data='../Inputs/NLO/'
  else:
    path='../Inputs/NLO/'

  filein=ROOT.TFile(path+filename,"READ")
  
 
  # define valuable to be fitted
  x = ROOT.RooRealVar("x", "x", 60, 120)
  x.setRange("FitRange", 61, 119)

  # for MC fit, no bkg shape needed. The signal shape is obtained using DY hist smeared with another distribution, e.g., a gaussian or Crystal ball extended gaussian (CBExGauss)

  if int(ismc):
    dy_pass=ROOT.TH1D()
    dy_fail=ROOT.TH1D()
    if unc=='puUp':
      filein.GetObject("TnP_mass_DYpass_pileUp",dy_pass)
      filein.GetObject("TnP_mass_DYfail_pileUp",dy_fail)
    elif unc=='puDo':
      filein.GetObject("TnP_mass_DYpass_pileDo",dy_pass)
      filein.GetObject("TnP_mass_DYfail_pileDo",dy_fail)
    else:
      filein.GetObject("TnP_mass_DYpass",dy_pass)
      filein.GetObject("TnP_mass_DYfail",dy_fail)

    DY_pass_error=ROOT.Double(0.)
    DY_fail_error=ROOT.Double(0.)
    DY_pass_total=dy_pass.IntegralAndError(1,60,DY_pass_error)
    DY_fail_total=dy_fail.IntegralAndError(1,60,DY_fail_error)
  
    # make RootDataHist for DY
    GenPass = ROOT.RooDataHist("GenPass","GenPass",ROOT.RooArgList(x),dy_pass)
    GenFail = ROOT.RooDataHist("GenFail","GenFail",ROOT.RooArgList(x),dy_fail)
  
    # make pdf for signal shape costruction 
    ZPassShape = ROOT.RooHistPdf("ZPassShape","ZPassShape",ROOT.RooArgSet(x), GenPass)
    ZFailShape = ROOT.RooHistPdf("ZFailShape","ZFailShape",ROOT.RooArgSet(x), GenFail)
  
    sig_pass = ROOT.RooFFTConvPdf()
    sig_fail = ROOT.RooFFTConvPdf()

    if sig_shape=='gauss':
      # gauss smearing of DY shape
      meanPass = ROOT.RooRealVar("meanPass", "meanPass", -0.0, -5.0, 5.0)
      sigmaPass = ROOT.RooRealVar("sigmaPass", "sigmaPass", 0.9,0.5,5.0)
      gaussPass = ROOT.RooGaussian("gaussPass", "gaussPass", x, meanPass, sigmaPass)

      meanFail = ROOT.RooRealVar("meanFail", "meanFail", -0.0, -5.0, 5.0)
      sigmaFail = ROOT.RooRealVar("sigmaFail", "sigmaFail", 0.9,0.5,5.0)
      gaussFail = ROOT.RooGaussian("gaussFail", "gaussFail", x, meanFail, sigmaFail)
      
      sig_pass = ROOT.RooFFTConvPdf("sigP", "signal shape", x, ZPassShape, gaussPass)
      sig_fail = ROOT.RooFFTConvPdf("sigF", "signal shape", x, ZFailShape, gaussFail)
  
    if sig_shape=='cbexgauss':
      # parameter of RooCBEXGaussian Shape, alternative pdf of signal
      meanP_altsig = ROOT.RooRealVar("meanP_altsig", "meanP_altsig", -0.0, -5.0, 5.0)
      sigmaP_altsig = ROOT.RooRealVar("sigmaP_altsig", "sigmaP_altsig", 1.0,0.7,6.0)
      alphaP_altsig = ROOT.RooRealVar("alphaP_altsig", "alphaP_altsig", 2.0,1.2,3.5)
      nP_altsig = ROOT.RooRealVar("nP_altsig", "nP_altsig", 3,-5,5)
      sigmaP2_altsig = ROOT.RooRealVar("sigmaP2_altsig", "sigmaP2_altsig", 1.5,0.5,6.0)
      sosP_altsig = ROOT.RooRealVar("sosP_altsig", "sosP_altsig", 1,0.5,5.0)
      sigmaExprP = ROOT.RooFormulaVar("sigmaExprP", "sqrt(sigmaP_altsig*sigmaP_altsig+sosP_altsig*sosP_altsig)", ROOT.RooArgList(sigmaP_altsig,sosP_altsig))
      sigmaExprP2 = ROOT.RooFormulaVar("sigmaExprP2","sqrt(sigmaP2_altsig*sigmaP2_altsig+sosP_altsig*sosP_altsig)", ROOT.RooArgList(sigmaP2_altsig,sosP_altsig))

      meanF_altsig = ROOT.RooRealVar("meanF_altsig", "meanF_altsig", -0.0, -5.0, 5.0)
      sigmaF_altsig = ROOT.RooRealVar("sigmaF_altsig", "sigmaF_altsig", 2.0,0.7,15.0)
      alphaF_altsig = ROOT.RooRealVar("alphaF_altsig", "alphaF_altsig", 2.0,1.2,3.5)
      nF_altsig = ROOT.RooRealVar("nF_altsig", "nF_altsig", 3,-5,5)
      sigmaF2_altsig = ROOT.RooRealVar("sigmaF2_altsig", "sigmaF2_altsig", 2.0,0.5,6.0)
      sosF_altsig = ROOT.RooRealVar("sosF_altsig", "sosF_altsig", 1,0.5,5.0)
      sigmaExprF = ROOT.RooFormulaVar("sigmaExprF","sqrt(sigmaF_altsig*sigmaF_altsig+sosF_altsig*sosF_altsig)", ROOT.RooArgList(sigmaF_altsig,sosF_altsig))
      sigmaExprF2 = ROOT.RooFormulaVar("sigmaExprF2","sqrt(sigmaF2_altsig*sigmaF2_altsig+sosF_altsig*sosF_altsig)", ROOT.RooArgList(sigmaF2_altsig,sosF_altsig))
  
      tailLeft = ROOT.RooRealVar("tailLeft", "tailLeft", 0.5,0.4,0.6)

      cbex_pass_Smear = RooCBExGaussShape("cbex_pass_Smear", "cbex_pass_Smear", x, meanP_altsig, sigmaExprP, alphaP_altsig, nP_altsig, sigmaExprP2, tailLeft)
      sig_pass = ROOT.RooFFTConvPdf("sigP", "signal shape", x, ZPassShape,cbex_pass_Smear)
      cbex_fail_Smear = RooCBExGaussShape("cbex_fail_Smear", "cbex_fail_Smear", x, meanF_altsig, sigmaExprF, alphaF_altsig, nF_altsig, sigmaExprF2, tailLeft)
      sig_fail = ROOT.RooFFTConvPdf("sigF", "signal shape", x, ZFailShape,cbex_fail_Smear)

    Pframe = x.frame(ROOT.RooFit.Title("passing probe"))
    Fframe = x.frame(ROOT.RooFit.Title("failing probe"))
  
    GenPass.plotOn(Pframe)
    rPass = sig_pass.fitTo(GenPass, ROOT.RooFit.Range("FitRange"), ROOT.RooFit.Save())
    sig_pass.plotOn(Pframe, ROOT.RooFit.LineColor(ROOT.kRed))

    GenFail.plotOn(Fframe)
    rFail = sig_fail.fitTo(GenFail, ROOT.RooFit.Range("FitRange"), ROOT.RooFit.Save())
    sig_fail.plotOn(Fframe, ROOT.RooFit.LineColor(ROOT.kRed))

    nTot=DY_pass_total+DY_fail_total
    eff=DY_pass_total/nTot
    e_eff = 1./(nTot*nTot) * sqrt( DY_pass_total*DY_pass_total* DY_fail_error*DY_fail_error + DY_fail_total*DY_fail_total * DY_pass_error*DY_pass_error )

#    c1=TCanvas("TnP","TnP",1200,600)
#    c1.Divide(3,1)
#    c1.cd(2)
#    Pframe.Draw()
#    c1.cd(3)
#    Fframe.Draw()
#
#    # Add text1, results, Eff and err
#    text1 = ROOT.TPaveText(0.05,0.75,0.3,0.95)
#    text1.SetFillColor(0)
#    text1.SetBorderSize(0)
#    text1.SetTextAlign(12)
#    text1.AddText('* MC Fit status:')
#    text1.AddText('passing: '+str(rPass.status())+', '+'failing: '+str(rFail.status()))
#    text1.AddText('* Eff = '+str('%1.4f'%eff)+' #pm '+str('%1.4f'%e_eff))
#    text1.SetTextSize(0.08)

#    # Add text2, results of fit parameters
#    text2 = ROOT.TPaveText(0.05,0.05,0.3,0.72)
#    text2.SetFillColor(0)
#    text2.SetBorderSize(0)
#    text2.SetTextAlign(12)

#    if sig_shape=='gauss':
#      text2.AddText('  --- parameters ')
#      text2.AddText('-meanP = '+str('%1.3f'%meanPass.getVal())+' #pm '+str('%1.3f'%meanPass.getError()))
#      text2.AddText('-sigmaP = '+str('%1.3f'%sigmaPass.getVal())+' #pm '+str('%1.3f'%sigmaPass.getError()))
#      text2.AddText('-nSigP = '+str('%1.3f'%DY_pass_total)+' #pm '+str('%1.3f'%DY_pass_error))
#      text2.AddText('-meanF = '+str('%1.3f'%meanFail.getVal())+' #pm '+str('%1.3f'%meanFail.getError()))
#      text2.AddText('-sigmaF = '+str('%1.3f'%sigmaFail.getVal())+' #pm '+str('%1.3f'%sigmaFail.getError()))
#      text2.AddText('-nSigF = '+str('%1.3f'%DY_fail_total)+' #pm '+str('%1.3f'%DY_fail_error))
#      text2.SetTextSize(0.06)
#
#    if sig_shape=='cbexgauss':
#      text2.AddText('  --- parameters ')
#      text2.AddText('-meanP = '+str('%1.3f'%meanP_altsig.getVal())+' #pm '+str('%1.3f'%meanP_altsig.getError()))
#      text2.AddText('-sigmaP = '+str('%1.3f'%sigmaP_altsig.getVal())+' #pm '+str('%1.3f'%sigmaP_altsig.getError()))
#      text2.AddText('-alphaP = '+str('%1.3f'%alphaP_altsig.getVal())+' #pm '+str('%1.3f'%alphaP_altsig.getError()))
#      text2.AddText('-nP = '+str('%1.3f'%nP_altsig.getVal())+' #pm '+str('%1.3f'%nP_altsig.getError()))
#      text2.AddText('-sigmaP2 = '+str('%1.3f'%sigmaP2_altsig.getVal())+' #pm '+str('%1.3f'%sigmaP2_altsig.getError()))
#      text2.AddText('-sosP = '+str('%1.3f'%sosP_altsig.getVal())+' #pm '+str('%1.3f'%sosP_altsig.getError()))
#      text2.AddText('-nSigP = '+str('%1.3f'%DY_pass_total)+' #pm '+str('%1.3f'%DY_pass_error))
#      text2.AddText('-meanF = '+str('%1.3f'%meanF_altsig.getVal())+' #pm '+str('%1.3f'%meanF_altsig.getError()))
#      text2.AddText('-sigmaF = '+str('%1.3f'%sigmaF_altsig.getVal())+' #pm '+str('%1.3f'%sigmaF_altsig.getError()))
#      text2.AddText('-alphaF = '+str('%1.3f'%alphaF_altsig.getVal())+' #pm '+str('%1.3f'%alphaF_altsig.getError()))
#      text2.AddText('-nF = '+str('%1.3f'%nF_altsig.getVal())+' #pm '+str('%1.3f'%nF_altsig.getError()))
#      text2.AddText('-sigmaF2 = '+str('%1.3f'%sigmaF2_altsig.getVal())+' #pm '+str('%1.3f'%sigmaF2_altsig.getError()))
#      text2.AddText('-sosF = '+str('%1.3f'%sosF_altsig.getVal())+' #pm '+str('%1.3f'%sosF_altsig.getError()))
#      text2.AddText('-nSigF = '+str('%1.3f'%DY_fail_total)+' #pm '+str('%1.3f'%DY_fail_error))
#      text2.SetTextSize(0.06)

#    c1.cd(1)
#    text1.Draw()
#    text2.Draw()
#    c1.SaveAs("mc_"+filename.strip('.root')+".png")

  else:
    # get data hist
    Mu_pass=ROOT.TH1D()
    Mu_fail=ROOT.TH1D()
    filein_data=ROOT.TFile()
    if unc=='LO':
      filein_data=ROOT.TFile(path_data+filename,"READ")
      filein_data.GetObject("TnP_mass_Mupass",Mu_pass)
      filein_data.GetObject("TnP_mass_Mufail",Mu_fail)
    else: 
      filein.GetObject("TnP_mass_Mupass",Mu_pass)
      filein.GetObject("TnP_mass_Mufail",Mu_fail)

    Mu_pass_error=ROOT.Double(0.)
    Mu_fail_error=ROOT.Double(0.)
    Mu_pass_total=Mu_pass.IntegralAndError(1,60,Mu_pass_error)
    Mu_fail_total=Mu_fail.IntegralAndError(1,60,Mu_fail_error)

    nSigP = ROOT.RooRealVar("nSigP","nSigP",0.9*Mu_pass_total,0.5*Mu_pass_total,1.5*Mu_pass_total)
    nBkgP = ROOT.RooRealVar("nBkgP","nBkgP",0.1*Mu_pass_total,0.,1.5*Mu_pass_total)
  
    nSigF = ROOT.RooRealVar("nSigF","nSigF",0.9*Mu_fail_total,0.5*Mu_fail_total,1.5*Mu_fail_total)
    nBkgF = ROOT.RooRealVar("nBkgF","nBkgF",0.1*Mu_fail_total,0.,1.5*Mu_fail_total)

    dy_pass=ROOT.TH1D()
    dy_fail=ROOT.TH1D()
    if unc=='puUp':
      filein.GetObject("TnP_mass_DYpass_pileUp",dy_pass)
      filein.GetObject("TnP_mass_DYfail_pileUp",dy_fail)
    elif unc=='puDo':
      filein.GetObject("TnP_mass_DYpass_pileDo",dy_pass)
      filein.GetObject("TnP_mass_DYfail_pileDo",dy_fail)
    else:
      filein.GetObject("TnP_mass_DYpass",dy_pass)
      filein.GetObject("TnP_mass_DYfail",dy_fail)

    DY_pass_error=ROOT.Double(0.)
    DY_fail_error=ROOT.Double(0.)
    DY_pass_total=dy_pass.IntegralAndError(1,60,DY_pass_error)
    DY_fail_total=dy_fail.IntegralAndError(1,60,DY_fail_error)
    sig_pass = ROOT.RooFFTConvPdf()
    sig_fail = ROOT.RooFFTConvPdf()
    bkgP = RooCMSShape()
    bkgF = RooCMSShape()

    GenPass = ROOT.RooDataHist("GenPass","GenPass",ROOT.RooArgList(x),dy_pass)
    GenFail = ROOT.RooDataHist("GenFail","GenFail",ROOT.RooArgList(x),dy_fail)

    DataPass = ROOT.RooDataHist("DataPass","DataPass",ROOT.RooArgList(x),Mu_pass)
    DataFail = ROOT.RooDataHist("DataFail","DataFail",ROOT.RooArgList(x),Mu_fail)

    ZPassShape = ROOT.RooHistPdf("ZPassShape","ZPassShape",ROOT.RooArgSet(x), GenPass)
    ZFailShape = ROOT.RooHistPdf("ZFailShape","ZFailShape",ROOT.RooArgSet(x), GenFail)

    if sig_shape=='gauss':
      # gauss smearing of DY shape
      meanPass = ROOT.RooRealVar("meanPass", "meanPass", -0.0, -5.0, 5.0)
      sigmaPass = ROOT.RooRealVar("sigmaPass", "sigmaPass", 0.9,0.5,5.0)
      gaussPass = ROOT.RooGaussian("gaussPass", "gaussPass", x, meanPass, sigmaPass)

      meanFail = ROOT.RooRealVar("meanFail", "meanFail", -0.0, -5.0, 5.0)
      sigmaFail = ROOT.RooRealVar("sigmaFail", "sigmaFail", 0.9,0.5,5.0)
      gaussFail = ROOT.RooGaussian("gaussFail", "gaussFail", x, meanFail, sigmaFail)
      
      sig_pass = ROOT.RooFFTConvPdf("sigP", "signal shape", x, ZPassShape, gaussPass)
      sig_fail = ROOT.RooFFTConvPdf("sigF", "signal shape", x, ZFailShape, gaussFail)

      if bkg_shape=='cmsshape':
        # parameter of RooCMSShape, pdf of background
        acmsP = ROOT.RooRealVar("acmsP", "acms", 60.,50.,80.)
        betaP = ROOT.RooRealVar("betaP", "beta", 0.05,0.01,0.08)
        gammaP = ROOT.RooRealVar("gammaP", "gamma", 0.1, -2, 2)
        peakP = ROOT.RooRealVar("peakP", "peak", 90.0)
      
        bkgP = RooCMSShape("bkgP", "bkg shape", x, acmsP, betaP, gammaP, peakP)
      
        acmsF = ROOT.RooRealVar("acmsF", "acms", 60.,50.,80.)
        betaF = ROOT.RooRealVar("betaF", "beta", 0.05,0.01,0.08)
        gammaF = ROOT.RooRealVar("gammaF", "gamma", 0.1, -2, 2)
        peakF = ROOT.RooRealVar("peakF", "peak", 90.0)
      
        bkgF = RooCMSShape("bkgF", "bkg shape", x, acmsF, betaF, gammaF, peakF)

      if bkg_shape=='expo':
        # parameter of Exponential, alternative pdf of background
        alphaP = ROOT.RooRealVar("alphaP", "alphaP", 0.,-5.,5.)
        alphaF = ROOT.RooRealVar("alphaF", "alphaF", 0.,-5.,5.)
        bkgP = ROOT.RooExponential("bkgP", "bkg shape", x, alphaP)
        bkgF = ROOT.RooExponential("bkgF", "bkg shape", x, alphaF)

      modelP=ROOT.RooAddPdf("modelP","modelP", ROOT.RooArgList(sig_pass,bkgP), ROOT.RooArgList(nSigP,nBkgP))
      modelF=ROOT.RooAddPdf("modelF","modelF", ROOT.RooArgList(sig_fail,bkgF), ROOT.RooArgList(nSigF,nBkgF))

    if sig_shape=='cbexgauss':
      # parameter of RooCBEXGaussian Shape, alternative pdf of signal
      meanP_altsig = ROOT.RooRealVar("meanP_altsig", "meanP_altsig", -0.0, -5.0, 5.0)
      sigmaP_altsig = ROOT.RooRealVar("sigmaP_altsig", "sigmaP_altsig", 1.0,0.7,6.0)
      alphaP_altsig = ROOT.RooRealVar("alphaP_altsig", "alphaP_altsig", 2.0,1.2,3.5)
      nP_altsig = ROOT.RooRealVar("nP_altsig", "nP_altsig", 3,-5,5)
      sigmaP2_altsig = ROOT.RooRealVar("sigmaP2_altsig", "sigmaP2_altsig", 1.5,0.5,6.0)
      sosP_altsig = ROOT.RooRealVar("sosP_altsig", "sosP_altsig", 1,0.5,5.0)
      sigmaExprP = ROOT.RooFormulaVar("sigmaExprP", "sqrt(sigmaP_altsig*sigmaP_altsig+sosP_altsig*sosP_altsig)", ROOT.RooArgList(sigmaP_altsig,sosP_altsig))
      sigmaExprP2 = ROOT.RooFormulaVar("sigmaExprP2","sqrt(sigmaP2_altsig*sigmaP2_altsig+sosP_altsig*sosP_altsig)", ROOT.RooArgList(sigmaP2_altsig,sosP_altsig))

      meanF_altsig = ROOT.RooRealVar("meanF_altsig", "meanF_altsig", -0.0, -5.0, 5.0)
      sigmaF_altsig = ROOT.RooRealVar("sigmaF_altsig", "sigmaF_altsig", 2.0,0.7,15.0)
      alphaF_altsig = ROOT.RooRealVar("alphaF_altsig", "alphaF_altsig", 2.0,1.2,3.5)
      nF_altsig = ROOT.RooRealVar("nF_altsig", "nF_altsig", 3,-5,5)
      sigmaF2_altsig = ROOT.RooRealVar("sigmaF2_altsig", "sigmaF2_altsig", 2.0,0.5,6.0)
      sosF_altsig = ROOT.RooRealVar("sosF_altsig", "sosF_altsig", 1,0.5,5.0)
      sigmaExprF = ROOT.RooFormulaVar("sigmaExprF","sqrt(sigmaF_altsig*sigmaF_altsig+sosF_altsig*sosF_altsig)", ROOT.RooArgList(sigmaF_altsig,sosF_altsig))
      sigmaExprF2 = ROOT.RooFormulaVar("sigmaExprF2","sqrt(sigmaF2_altsig*sigmaF2_altsig+sosF_altsig*sosF_altsig)", ROOT.RooArgList(sigmaF2_altsig,sosF_altsig))
  
      tailLeft = ROOT.RooRealVar("tailLeft", "tailLeft", 0.5,0.4,0.6)

      cbex_pass_Smear = RooCBExGaussShape("cbex_pass_Smear", "cbex_pass_Smear", x, meanP_altsig, sigmaExprP, alphaP_altsig, nP_altsig, sigmaExprP2, tailLeft)
      sig_pass = ROOT.RooFFTConvPdf("sigP", "signal shape", x, ZPassShape,cbex_pass_Smear)
      cbex_fail_Smear = RooCBExGaussShape("cbex_fail_Smear", "cbex_fail_Smear", x, meanF_altsig, sigmaExprF, alphaF_altsig, nF_altsig, sigmaExprF2, tailLeft)
      sig_fail = ROOT.RooFFTConvPdf("sigF", "signal shape", x, ZFailShape,cbex_fail_Smear)

      acmsP = ROOT.RooRealVar("acmsP", "acms", 60.,50.,80.)
      betaP = ROOT.RooRealVar("betaP", "beta", 0.05,0.01,0.08)
      gammaP = ROOT.RooRealVar("gammaP", "gamma", 0.1, -2, 2)
      peakP = ROOT.RooRealVar("peakP", "peak", 90.0)
      
      bkgP = RooCMSShape("bkgP", "bkg shape", x, acmsP, betaP, gammaP, peakP)
      
      acmsF = ROOT.RooRealVar("acmsF", "acms", 60.,50.,80.)
      betaF = ROOT.RooRealVar("betaF", "beta", 0.05,0.01,0.08)
      gammaF = ROOT.RooRealVar("gammaF", "gamma", 0.1, -2, 2)
      peakF = ROOT.RooRealVar("peakF", "peak", 90.0)
      
      bkgF = RooCMSShape("bkgF", "bkg shape", x, acmsF, betaF, gammaF, peakF)

      modelP=ROOT.RooAddPdf("modelP","modelP", ROOT.RooArgList(sig_pass,bkgP), ROOT.RooArgList(nSigP,nBkgP))
      modelF=ROOT.RooAddPdf("modelF","modelF", ROOT.RooArgList(sig_fail,bkgF), ROOT.RooArgList(nSigF,nBkgF))

    Pframe = x.frame(ROOT.RooFit.Title("passing probe"))
    Fframe = x.frame(ROOT.RooFit.Title("failing probe"))

    DataPass.plotOn(Pframe)
    rPass = modelP.fitTo(DataPass, ROOT.RooFit.Range("FitRange"), ROOT.RooFit.Save())
    modelP.plotOn(Pframe, ROOT.RooFit.LineColor(ROOT.kRed))
    modelP.plotOn(Pframe, ROOT.RooFit.Components("bkgP"),ROOT.RooFit.LineColor(ROOT.kBlue),ROOT.RooFit.LineStyle(ROOT.kDashed))
    DataFail.plotOn(Fframe)
    rFail = modelF.fitTo(DataFail, ROOT.RooFit.Range("FitRange"), ROOT.RooFit.Save())
    modelF.plotOn(Fframe, ROOT.RooFit.LineColor(ROOT.kRed))
    modelF.plotOn(Fframe, ROOT.RooFit.Components("bkgF"),ROOT.RooFit.LineColor(ROOT.kBlue),ROOT.RooFit.LineStyle(ROOT.kDashed))
    nTot=nSigP.getVal()+nSigF.getVal()
    eff=nSigP.getVal()/nTot
    e_eff = 1./(nTot*nTot)*sqrt(nSigP.getVal()*nSigP.getVal()*nSigF.getError()*nSigF.getError() + nSigF.getVal()*nSigF.getVal() * nSigP.getError()*nSigP.getError() )

  c1=TCanvas("TnP","TnP",1200,600)
  c1.Divide(3,1)
  c1.cd(2)
  Pframe.Draw()
  c1.cd(3)
  Fframe.Draw()

  # Add text1, results, Eff and err
  text1 = ROOT.TPaveText(0.05,0.75,0.3,0.95)
  text1.SetFillColor(0)
  text1.SetBorderSize(0)
  text1.SetTextAlign(12)
  if int(ismc):
    text1.AddText('* MC Fit status:')
    text1.AddText('passing: '+str(rPass.status())+', '+'failing: '+str(rFail.status()))
    text1.AddText('* Eff = '+str('%1.4f'%eff)+' #pm '+str('%1.4f'%e_eff))
    text1.SetTextSize(0.08)
  else:
    text1.AddText('* Data Fit status:')
    text1.AddText('passing: '+str(rPass.status())+', '+'failing: '+str(rFail.status()))
    text1.AddText('* Eff = '+str('%1.4f'%eff)+' #pm '+str('%1.4f'%e_eff))
    text1.SetTextSize(0.08)
   
  # Add text2, results of fit parameters
  text2 = ROOT.TPaveText(0.05,0.05,0.3,0.72)
  text2.SetFillColor(0)
  text2.SetBorderSize(0)
  text2.SetTextAlign(12)

  if int(ismc):
    if sig_shape=='gauss':
      text2.AddText('  --- parameters ')
      text2.AddText('-meanP = '+str('%1.3f'%meanPass.getVal())+' #pm '+str('%1.3f'%meanPass.getError()))
      text2.AddText('-sigmaP = '+str('%1.3f'%sigmaPass.getVal())+' #pm '+str('%1.3f'%sigmaPass.getError()))
      text2.AddText('-nSigP = '+str('%1.3f'%DY_pass_total)+' #pm '+str('%1.3f'%DY_pass_error))
      text2.AddText('-meanF = '+str('%1.3f'%meanFail.getVal())+' #pm '+str('%1.3f'%meanFail.getError()))
      text2.AddText('-sigmaF = '+str('%1.3f'%sigmaFail.getVal())+' #pm '+str('%1.3f'%sigmaFail.getError()))
      text2.AddText('-nSigF = '+str('%1.3f'%DY_fail_total)+' #pm '+str('%1.3f'%DY_fail_error))
      text2.SetTextSize(0.06)

    if sig_shape=='cbexgauss':
      text2.AddText('  --- parameters ')
      text2.AddText('-meanP = '+str('%1.3f'%meanP_altsig.getVal())+' #pm '+str('%1.3f'%meanP_altsig.getError()))
      text2.AddText('-sigmaP = '+str('%1.3f'%sigmaP_altsig.getVal())+' #pm '+str('%1.3f'%sigmaP_altsig.getError()))
      text2.AddText('-alphaP = '+str('%1.3f'%alphaP_altsig.getVal())+' #pm '+str('%1.3f'%alphaP_altsig.getError()))
      text2.AddText('-nP = '+str('%1.3f'%nP_altsig.getVal())+' #pm '+str('%1.3f'%nP_altsig.getError()))
      text2.AddText('-sigmaP2 = '+str('%1.3f'%sigmaP2_altsig.getVal())+' #pm '+str('%1.3f'%sigmaP2_altsig.getError()))
      text2.AddText('-sosP = '+str('%1.3f'%sosP_altsig.getVal())+' #pm '+str('%1.3f'%sosP_altsig.getError()))
      text2.AddText('-nSigP = '+str('%1.3f'%DY_pass_total)+' #pm '+str('%1.3f'%DY_pass_error))
      text2.AddText('-meanF = '+str('%1.3f'%meanF_altsig.getVal())+' #pm '+str('%1.3f'%meanF_altsig.getError()))
      text2.AddText('-sigmaF = '+str('%1.3f'%sigmaF_altsig.getVal())+' #pm '+str('%1.3f'%sigmaF_altsig.getError()))
      text2.AddText('-alphaF = '+str('%1.3f'%alphaF_altsig.getVal())+' #pm '+str('%1.3f'%alphaF_altsig.getError()))
      text2.AddText('-nF = '+str('%1.3f'%nF_altsig.getVal())+' #pm '+str('%1.3f'%nF_altsig.getError()))
      text2.AddText('-sigmaF2 = '+str('%1.3f'%sigmaF2_altsig.getVal())+' #pm '+str('%1.3f'%sigmaF2_altsig.getError()))
      text2.AddText('-sosF = '+str('%1.3f'%sosF_altsig.getVal())+' #pm '+str('%1.3f'%sosF_altsig.getError()))
      text2.AddText('-nSigF = '+str('%1.3f'%DY_fail_total)+' #pm '+str('%1.3f'%DY_fail_error))
      text2.SetTextSize(0.06)

  else:
    if sig_shape=='gauss':
      if bkg_shape=='cmsshape':
        text2.AddText('  --- parameters ')
        text2.AddText('-nBkgP = '+str('%1.3f'%nBkgP.getVal())+' #pm '+str('%1.3f'%nBkgP.getError()))
        text2.AddText('-nSigP = '+str('%1.3f'%nSigP.getVal())+' #pm '+str('%1.3f'%nSigP.getError()))
        text2.AddText('-meanP = '+str('%1.3f'%meanPass.getVal())+' #pm '+str('%1.3f'%meanPass.getError()))
        text2.AddText('-sigmaP = '+str('%1.3f'%sigmaPass.getVal())+' #pm '+str('%1.3f'%sigmaPass.getError()))
        text2.AddText('-acmsP = '+str('%1.3f'%acmsP.getVal())+' #pm '+str('%1.3f'%acmsP.getError()))
        text2.AddText('-betaP = '+str('%1.3f'%betaP.getVal())+' #pm '+str('%1.3f'%betaP.getError()))
        text2.AddText('-gammaP = '+str('%1.3f'%gammaP.getVal())+' #pm '+str('%1.3f'%gammaP.getError()))
        text2.AddText('-nBkgF = '+str('%1.3f'%nBkgF.getVal())+' #pm '+str('%1.3f'%nBkgF.getError()))
        text2.AddText('-nSigF = '+str('%1.3f'%nSigF.getVal())+' #pm '+str('%1.3f'%nSigF.getError()))
        text2.AddText('-meanF = '+str('%1.3f'%meanFail.getVal())+' #pm '+str('%1.3f'%meanFail.getError()))
        text2.AddText('-sigmaF = '+str('%1.3f'%sigmaFail.getVal())+' #pm '+str('%1.3f'%sigmaFail.getError()))
        text2.AddText('-acmsF = '+str('%1.3f'%acmsF.getVal())+' #pm '+str('%1.3f'%acmsF.getError()))
        text2.AddText('-betaF = '+str('%1.3f'%betaF.getVal())+' #pm '+str('%1.3f'%betaF.getError()))
        text2.AddText('-gammaF = '+str('%1.3f'%gammaF.getVal())+' #pm '+str('%1.3f'%gammaF.getError()))
        text2.SetTextSize(0.05)
      if bkg_shape=='expo':
        text2.AddText('  --- parameters ')
        text2.AddText('-nBkgP = '+str('%1.3f'%nBkgP.getVal())+' #pm '+str('%1.3f'%nBkgP.getError()))
        text2.AddText('-nSigP = '+str('%1.3f'%nSigP.getVal())+' #pm '+str('%1.3f'%nSigP.getError()))
        text2.AddText('-meanP = '+str('%1.3f'%meanPass.getVal())+' #pm '+str('%1.3f'%meanPass.getError()))
        text2.AddText('-sigmaP = '+str('%1.3f'%sigmaPass.getVal())+' #pm '+str('%1.3f'%sigmaPass.getError()))
        text2.AddText('-alphaP = '+str('%1.3f'%alphaP.getVal())+' #pm '+str('%1.3f'%alphaP.getError()))
        text2.AddText('-nBkgF = '+str('%1.3f'%nBkgF.getVal())+' #pm '+str('%1.3f'%nBkgF.getError()))
        text2.AddText('-nSigF = '+str('%1.3f'%nSigF.getVal())+' #pm '+str('%1.3f'%nSigF.getError()))
        text2.AddText('-meanF = '+str('%1.3f'%meanFail.getVal())+' #pm '+str('%1.3f'%meanFail.getError()))
        text2.AddText('-sigmaF = '+str('%1.3f'%sigmaFail.getVal())+' #pm '+str('%1.3f'%sigmaFail.getError()))
        text2.AddText('-alphaF = '+str('%1.3f'%alphaF.getVal())+' #pm '+str('%1.3f'%alphaF.getError()))
        text2.SetTextSize(0.05)
    if sig_shape=='cbexgauss':
      text2.AddText('  --- parameters ')
      text2.AddText('-nBkgP = '+str('%1.3f'%nBkgP.getVal())+' #pm '+str('%1.3f'%nBkgP.getError()))
      text2.AddText('-nSigP = '+str('%1.3f'%nSigP.getVal())+' #pm '+str('%1.3f'%nSigP.getError()))
      text2.AddText('-meanP = '+str('%1.3f'%meanP_altsig.getVal())+' #pm '+str('%1.3f'%meanP_altsig.getError()))
      text2.AddText('-sigmaP = '+str('%1.3f'%sigmaP_altsig.getVal())+' #pm '+str('%1.3f'%sigmaP_altsig.getError()))
      text2.AddText('-alphaP = '+str('%1.3f'%alphaP_altsig.getVal())+' #pm '+str('%1.3f'%alphaP_altsig.getError()))
      text2.AddText('-nP = '+str('%1.3f'%nP_altsig.getVal())+' #pm '+str('%1.3f'%nP_altsig.getError()))
      text2.AddText('-sigmaP2 = '+str('%1.3f'%sigmaP2_altsig.getVal())+' #pm '+str('%1.3f'%sigmaP2_altsig.getError()))
      text2.AddText('-sosP = '+str('%1.3f'%sosP_altsig.getVal())+' #pm '+str('%1.3f'%sosP_altsig.getError()))
      text2.AddText('-nSigP = '+str('%1.3f'%DY_pass_total)+' #pm '+str('%1.3f'%DY_pass_error))
      text2.AddText('-acmsP = '+str('%1.3f'%acmsP.getVal())+' #pm '+str('%1.3f'%acmsP.getError()))
      text2.AddText('-betaP = '+str('%1.3f'%betaP.getVal())+' #pm '+str('%1.3f'%betaP.getError()))
      text2.AddText('-gammaP = '+str('%1.3f'%gammaP.getVal())+' #pm '+str('%1.3f'%gammaP.getError()))
      text2.AddText('-nBkgF = '+str('%1.3f'%nBkgF.getVal())+' #pm '+str('%1.3f'%nBkgF.getError()))
      text2.AddText('-nSigF = '+str('%1.3f'%nSigF.getVal())+' #pm '+str('%1.3f'%nSigF.getError()))
      text2.AddText('-meanF = '+str('%1.3f'%meanF_altsig.getVal())+' #pm '+str('%1.3f'%meanF_altsig.getError()))
      text2.AddText('-sigmaF = '+str('%1.3f'%sigmaF_altsig.getVal())+' #pm '+str('%1.3f'%sigmaF_altsig.getError()))
      text2.AddText('-alphaF = '+str('%1.3f'%alphaF_altsig.getVal())+' #pm '+str('%1.3f'%alphaF_altsig.getError()))
      text2.AddText('-nF = '+str('%1.3f'%nF_altsig.getVal())+' #pm '+str('%1.3f'%nF_altsig.getError()))
      text2.AddText('-sigmaF2 = '+str('%1.3f'%sigmaF2_altsig.getVal())+' #pm '+str('%1.3f'%sigmaF2_altsig.getError()))
      text2.AddText('-sosF = '+str('%1.3f'%sosF_altsig.getVal())+' #pm '+str('%1.3f'%sosF_altsig.getError()))
      text2.AddText('-nSigF = '+str('%1.3f'%DY_fail_total)+' #pm '+str('%1.3f'%DY_fail_error))
      text2.AddText('-acmsF = '+str('%1.3f'%acmsF.getVal())+' #pm '+str('%1.3f'%acmsF.getError()))
      text2.AddText('-betaF = '+str('%1.3f'%betaF.getVal())+' #pm '+str('%1.3f'%betaF.getError()))
      text2.AddText('-gammaF = '+str('%1.3f'%gammaF.getVal())+' #pm '+str('%1.3f'%gammaF.getError()))
      text2.SetTextSize(0.04)


  c1.cd(1)
  text1.Draw()
  text2.Draw()
  if int(ismc):
    c1.SaveAs("mc_"+filename.strip('.root')+".png")
  else:
    c1.SaveAs("data_"+filename.strip('.root')+".png")

  return eff, e_eff

if __name__ == "__main__":

  tdptbin=array('d',[10,20,25,40,50,60,120,500])
  tdptbin_plain=array('d',[1,2,3,4,5,6,7,8])
  tdptbinname=['10~20','20~25','25~40','40~50','50~60','60~120','120~500']
  tdetabin=array('d',[0.0,0.9,1.2,2.1,2.4])

  hist_eff = ROOT.TH2D('MuIDEff', 'MuIDEff', 4, tdetabin, 7, tdptbin_plain)
  hist_eff.Sumw2()
  hist_eff.SetStats(0)
  hist_eff.GetXaxis().SetTitle('Muon #||{#eta}')
  hist_eff.GetYaxis().SetTitle('Muon P_{T} [GeV]')
  hist_eff.SetTitle('')
  for ib in range(1,8):
    hist_eff.GetYaxis().SetBinLabel(ib,tdptbinname[ib-1])

  ptbinnames=['Pt10To20','Pt20To25','Pt25To40','Pt40To50','Pt50To60','Pt60To120','Pt120To500']
  etabinnames=['Etam0p0Top0p9','Etap0p9Top1p2','Etap1p2Top2p1','Etap2p1Top2p4']

  path='./'
  if unc=='LO':
    path='../Inputs/LO/'
  else:
    path='../Inputs/NLO/'

  for files in os.listdir(path):
    if not (files.startswith('Pt') and files.endswith('.root')):continue
    eff,eff_err = Eff(files, ismc, sig_shape, bkg_shape, unc)

    if ptbinnames[0] in files:
      if etabinnames[0] in files:
        hist_eff.SetBinContent(1,1,eff)
        hist_eff.SetBinError(1,1,eff_err)
      if etabinnames[1] in files:
        hist_eff.SetBinContent(2,1,eff)
        hist_eff.SetBinError(2,1,eff_err)
      if etabinnames[2] in files:
        hist_eff.SetBinContent(3,1,eff)
        hist_eff.SetBinError(3,1,eff_err)
      if etabinnames[3] in files:
        hist_eff.SetBinContent(4,1,eff)
        hist_eff.SetBinError(4,1,eff_err)
    if ptbinnames[1] in files:
      if etabinnames[0] in files:
        hist_eff.SetBinContent(1,2,eff)
        hist_eff.SetBinError(1,2,eff_err)
      if etabinnames[1] in files:
        hist_eff.SetBinContent(2,2,eff)
        hist_eff.SetBinError(2,2,eff_err)
      if etabinnames[2] in files:
        hist_eff.SetBinContent(3,2,eff)
        hist_eff.SetBinError(3,2,eff_err)
      if etabinnames[3] in files:
        hist_eff.SetBinContent(4,2,eff)
        hist_eff.SetBinError(4,2,eff_err)
    if ptbinnames[2] in files:
      if etabinnames[0] in files:
        hist_eff.SetBinContent(1,3,eff)
        hist_eff.SetBinError(1,3,eff_err)
      if etabinnames[1] in files:
        hist_eff.SetBinContent(2,3,eff)
        hist_eff.SetBinError(2,3,eff_err)
      if etabinnames[2] in files:
        hist_eff.SetBinContent(3,3,eff)
        hist_eff.SetBinError(3,3,eff_err)
      if etabinnames[3] in files:
        hist_eff.SetBinContent(4,3,eff)
        hist_eff.SetBinError(4,3,eff_err)
    if ptbinnames[3] in files:
      if etabinnames[0] in files:
        hist_eff.SetBinContent(1,4,eff)
        hist_eff.SetBinError(1,4,eff_err)
      if etabinnames[1] in files:
        hist_eff.SetBinContent(2,4,eff)
        hist_eff.SetBinError(2,4,eff_err)
      if etabinnames[2] in files:
        hist_eff.SetBinContent(3,4,eff)
        hist_eff.SetBinError(3,4,eff_err)
      if etabinnames[3] in files:
        hist_eff.SetBinContent(4,4,eff)
        hist_eff.SetBinError(4,4,eff_err)
    if ptbinnames[4] in files:
      if etabinnames[0] in files:
        hist_eff.SetBinContent(1,5,eff)
        hist_eff.SetBinError(1,5,eff_err)
      if etabinnames[1] in files:
        hist_eff.SetBinContent(2,5,eff)
        hist_eff.SetBinError(2,5,eff_err)
      if etabinnames[2] in files:
        hist_eff.SetBinContent(3,5,eff)
        hist_eff.SetBinError(3,5,eff_err)
      if etabinnames[3] in files:
        hist_eff.SetBinContent(4,5,eff)
        hist_eff.SetBinError(4,5,eff_err)
    if ptbinnames[5] in files:
      if etabinnames[0] in files:
        hist_eff.SetBinContent(1,6,eff)
        hist_eff.SetBinError(1,6,eff_err)
      if etabinnames[1] in files:
        hist_eff.SetBinContent(2,6,eff)
        hist_eff.SetBinError(2,6,eff_err)
      if etabinnames[2] in files:
        hist_eff.SetBinContent(3,6,eff)
        hist_eff.SetBinError(3,6,eff_err)
      if etabinnames[3] in files:
        hist_eff.SetBinContent(4,6,eff)
        hist_eff.SetBinError(4,6,eff_err)
    if ptbinnames[6] in files:
      if etabinnames[0] in files:
        hist_eff.SetBinContent(1,7,eff)
        hist_eff.SetBinError(1,7,eff_err)
      if etabinnames[1] in files:
        hist_eff.SetBinContent(2,7,eff)
        hist_eff.SetBinError(2,7,eff_err)
      if etabinnames[2] in files:
        hist_eff.SetBinContent(3,7,eff)
        hist_eff.SetBinError(3,7,eff_err)
      if etabinnames[3] in files:
        hist_eff.SetBinContent(4,7,eff)
        hist_eff.SetBinError(4,7,eff_err)

  c1 = TCanvas()
  pad1 = TPad()
  pad1.Draw()
  pad1.cd()
  hist_eff.Draw('COLZ TEXT E')
  CMSstyle.SetStyle(pad1)
  pad1.SetRightMargin(0.15)
  c1.SetGridx(False);
  c1.SetGridy(False);
  if int(ismc):
    outname='mc_eff_'+sig_shape+'_'+unc
  else:
    outname='data_eff_'+sig_shape+'_'+unc
  c1.SaveAs(outname+'.png')
  c1.SaveAs(outname+'.pdf')
  pad1.Close()

  fout = ROOT.TFile('output.root','recreate')
  fout.cd()
  hist_eff.Write()
  fout.Close()
