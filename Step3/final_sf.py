import ROOT,sys,math
from array import array
from math import sqrt
import CMSTDRStyle
CMSTDRStyle.setTDRStyle().cd()
import CMSstyle

path='../Effi/'

data_f1=ROOT.TFile(path+'data_eff_gauss_nominal/output.root','READ')
print(data_f1)
#data_f2=ROOT.TFile(path+'data_eff_gauss_LO/output.root','READ')
data_f3=ROOT.TFile(path+'data_eff_gauss_puUp/output.root','READ')
data_f4=ROOT.TFile(path+'data_eff_gauss_puDo/output.root','READ')
data_f5=ROOT.TFile(path+'data_eff_gauss_expo_nominal/output.root','READ')
data_f6=ROOT.TFile(path+'data_eff_cbexgauss_nominal/output.root','READ')

data_h1 = data_f1.Get('MuIDEff')
#data_h2 = data_f2.Get('MuIDEff')
data_h3 = data_f3.Get('MuIDEff')
data_h4 = data_f4.Get('MuIDEff')
data_h5 = data_f5.Get('MuIDEff')
data_h6 = data_f6.Get('MuIDEff')

mc_f1=ROOT.TFile(path+'mc_eff_gauss_nominal/output.root','READ')
#mc_f2=ROOT.TFile(path+'mc_eff_gauss_LO/output.root','READ')
mc_f3=ROOT.TFile(path+'mc_eff_gauss_puUp/output.root','READ')
mc_f4=ROOT.TFile(path+'mc_eff_gauss_puDo/output.root','READ')
mc_f5=ROOT.TFile(path+'mc_eff_cbexgauss_nominal/output.root','READ')

mc_h1 = mc_f1.Get('MuIDEff')
#mc_h2 = mc_f2.Get('MuIDEff')
mc_h3 = mc_f3.Get('MuIDEff')
mc_h4 = mc_f4.Get('MuIDEff')
mc_h5 = mc_f5.Get('MuIDEff')

data_h1.Divide(mc_h1) #nominal
#data_h2.Divide(mc_h2) #LO
data_h3.Divide(mc_h3) #nominal + puUp
data_h4.Divide(mc_h4) #nominal + puDo
data_h5.Divide(mc_h1) #nominal + expo
data_h6.Divide(mc_h5) #alt sig

h2_final=data_h1.Clone('plain_SFs')
data_h1.SetNameTitle('NLODY_SF','NLODY_SF')
#data_h2.SetNameTitle('LODY_SF','LODY_SF')
data_h3.SetNameTitle('puUp_SF','puUp_SF')
data_h4.SetNameTitle('puDo_SF','puDo_SF')
data_h5.SetNameTitle('expobkg_SF','expobkg_SF')
data_h6.SetNameTitle('altsig_SF','altsig_SF')

for ix in range(1,data_h1.GetNbinsX()+1):
  for iy in range(1,data_h1.GetNbinsY()+1):
#    err1=abs(data_h1.GetBinContent(ix,iy) - data_h2.GetBinContent(ix,iy))
    err2=abs(data_h1.GetBinContent(ix,iy) - data_h5.GetBinContent(ix,iy))
    err3=0.5*abs(data_h4.GetBinContent(ix,iy) - data_h3.GetBinContent(ix,iy))
    err4=abs(data_h1.GetBinContent(ix,iy) - data_h6.GetBinContent(ix,iy))
    err5=data_h1.GetBinError(ix,iy)
#    h2_final.SetBinError(ix,iy,math.sqrt(err1*err1+err2*err2+err3*err3+err4*err4+err5*err5))
    h2_final.SetBinError(ix,iy,math.sqrt(err2*err2+err3*err3+err4*err4+err5*err5))

tdptbin=array('d',[10,20,25,40,50,60,120,500])
tdptbinname=['10~20','20~25','25~40','40~50','50-60','60~120','120~Inf']
tdetabin=array('d',[0.0, 0.9, 1.2, 2.1, 2.4])

h2_SF = ROOT.TH2D('MuIDSF', 'MuIDSF', 4, tdetabin, 7, tdptbin)
h2_SF.SetStats(0)
h2_SF.GetYaxis().SetTitle('Muon P_{T} [GeV]')
h2_SF.GetXaxis().SetTitle('Muon #||{#eta}')
h2_SF.SetTitle('')
for ib in range(1,8):
  h2_SF.GetYaxis().SetBinLabel(ib,tdptbinname[ib-1])

for ix in range(1,h2_final.GetNbinsX()+1):
  for iy in range(1,h2_final.GetNbinsY()+1):
    h2_SF.SetBinContent(ix,iy,h2_final.GetBinContent(ix,iy))
    h2_SF.SetBinError(ix,iy,h2_final.GetBinError(ix,iy))

c1 = ROOT.TCanvas()
pad1 = ROOT.TPad()
pad1.Draw()
pad1.cd()
h2_final.GetXaxis().SetTitleSize(0.04)
h2_final.GetYaxis().SetTitleSize(0.04)
h2_final.GetXaxis().SetTitleOffset(1.3)
h2_final.GetYaxis().SetTitleOffset(1.9)
h2_final.Draw('COLZ TEXT E')
CMSstyle.SetStyle(pad1)
pad1.SetRightMargin(0.15)
c1.SetGridx(False);
c1.SetGridy(False);
c1.SaveAs('./SF.png')
c1.SaveAs('./SF.pdf')
pad1.Close()

aout=ROOT.TFile.Open('SF.root','recreate')
h2_SF.Write()
h2_final.Write()
data_h1.Write()
#data_h2.Write()
data_h3.Write()
data_h4.Write()
data_h5.Write()
data_h6.Write()
aout.Close()
