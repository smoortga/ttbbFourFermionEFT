import ROOT
from argparse import ArgumentParser
import os


parser = ArgumentParser()

parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttbb_DiLepton_ForValidationPerCouplings_fixed", help='path to the directory of all restricted processes')
parser.add_argument('--run', default = "run_02", help='which run to use')
parser.add_argument('--coupling', default = "cQb1", help='which coupling to use')

args = parser.parse_args()

"""
make sure to do the following before running this script!
cdir=$(pwd)
export LD_LIBRARY_PATH=${cdir}:${cdir}/Delphes:$LD_LIBRARY_PATH
"""

def convertToLatex(c):
    if len(c) == 4: return c[0].upper()+"_{"+c[1:3]+"}^{"+c[3:4]+"}"
    else: 
        print "couldn't convert"
        return c   


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

ROOT.gSystem.Load("libDelphes")
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')

f_SM = ROOT.TFile("%s/%s/%s_0p0/Events/%s/tag_1_delphes_events.root"%(args.InputDir,args.coupling,args.coupling,args.run))
f_cQb1_10 = ROOT.TFile("%s/%s/%s_10p0/Events/%s/tag_1_delphes_events.root"%(args.InputDir,args.coupling,args.coupling,args.run))
f_cQb1_20 = ROOT.TFile("%s/%s/%s_20p0/Events/%s/tag_1_delphes_events.root"%(args.InputDir,args.coupling,args.coupling,args.run))

t_SM = f_SM.Get("Delphes")
t_cQb1_10 = f_cQb1_10.Get("Delphes")
t_cQb1_20 = f_cQb1_20.Get("Delphes")

xmin = 0
xmax = 3
nbins = 30
h_SM = ROOT.TH1D("h_SM",";H_{T} [TeV];a.u.",nbins,xmin,xmax)
h_cQb1_10 = ROOT.TH1D("h_cQb1_10",";H_{T} [TeV];a.u.",nbins,xmin,xmax)
h_cQb1_20 = ROOT.TH1D("h_cQb1_20",";H_{T} [TeV];a.u.",nbins,xmin,xmax)

t_SM.Draw("ScalarHT.HT/1000>>h_SM")
t_cQb1_10.Draw("ScalarHT.HT/1000>>h_cQb1_10")
t_cQb1_20.Draw("ScalarHT.HT/1000>>h_cQb1_20")

c = ROOT.TCanvas("c","c",800,500)
c.SetMargin(0.15,0.1,0.15,0.1)
h_cQb1_20.Draw("hist")
h_cQb1_10.Draw("hist same")
h_SM.Draw("hist same")

#style
h_cQb1_20.SetLineColor(2)
h_cQb1_20.SetLineWidth(2)
h_cQb1_20.GetYaxis().SetRangeUser(0,1500)
h_cQb1_20.GetYaxis().SetTickSize(0)
h_cQb1_20.GetYaxis().SetLabelSize(0)
h_cQb1_20.GetYaxis().SetTitleSize(0.06)
h_cQb1_20.GetYaxis().SetTitleOffset(0.9)
h_cQb1_20.GetXaxis().SetRangeUser(0,2.900)
h_cQb1_20.GetXaxis().SetLabelSize(0.055)
h_cQb1_20.GetXaxis().SetTitleSize(0.06)
h_cQb1_20.GetXaxis().SetTitleOffset(1.2)

h_cQb1_10.SetLineColor(4)
h_cQb1_10.SetLineWidth(2)

h_SM.SetLineColor(1)
h_SM.SetLineWidth(2)

#recast xaxis
#xmin = h->GetXaxis()->GetXmin(); double xmax = h->GetXaxis()->GetXmax(); h->SetLimits(xmin/1000,xmax/1000); 

l = ROOT.TLegend(0.55,0.3,0.85,0.88)
l.SetBorderSize(0)
l.AddEntry(h_SM,"SM only       ","l")
l.AddEntry(h_cQb1_10,"#splitline{SM + EFT     }{%s = 10 TeV^{-2}}"%convertToLatex("cQb1"),"l")
l.AddEntry(h_cQb1_20,"#splitline{SM + EFT     }{%s = 20 TeV^{-2}}"%convertToLatex("cQb1"),"l")

l.Draw("same")

if not os.path.isdir(args.InputDir + "/LimitsVsMcut"): os.mkdir(args.InputDir + "/LimitsVsMcut")
    
c.SaveAs(args.InputDir + "/LimitsVsMcut/Plot_HT.png")
c.SaveAs(args.InputDir + "/LimitsVsMcut/Plot_HT.pdf")
