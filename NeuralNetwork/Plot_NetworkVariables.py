import ROOT
from argparse import ArgumentParser
import os
import sys


parser = ArgumentParser()

parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Training_fixed", help='path to the directory of all restricted processes')
#parser.add_argument('--run', default = "run_02", help='which run to use')
#parser.add_argument('--coupling', default = "cQb1", help='which coupling to use')

args = parser.parse_args()

axis_name_dict = {
    'deltaR_c1c2':['#Delta R(add. b jets)',0,5],
    'pT_l1':['p_{T} first lepton [GeV]',0,400],
    'm_l1l2':['M_{inv}(leptons) [GeV]',0,1000],
    'm_c1c2b1b2':['M_{inv}(4b) [GeV]',0,2000],
    'pT_l2':['p_{T} second lepton [GeV]',0,400],
    'deltaR_l1l2':['#Delta R(leptons)',0,5],
    'm_b2l2':['M_{inv}(b2,l2) [GeV]',0,500],
    'm_c1c2':['M_{inv}(add. b jets) [GeV]',0,1000],
    'pT_b1':['p_{T} first b jet [GeV]',0,500],
    'm_b1l1':['M_{inv}(b1,l1) [GeV]',0,500],
    'deltaR_b2l2':['#Delta R(b2,l2)',0,5],
    'm_b1b2':['M_{inv}(b jets) [GeV]',0,1000],
    'm_c1c2b1b2l1l2':['M_{inv}(leptons and b jets) [GeV]',0,3000],
    'pT_c1':['p_{T} first add. b jet [GeV]',0,400],
    'pT_c2':['p_{T} second add. b jet [GeV]',0,400],
    'pT_b2':['p_{T} second b jet [GeV]',0,400],
    'deltaR_b1l1':['#Delta R(b1,l1)',0,5],
    'deltaR_b1b2':['#Delta R(b jets)',0,5]
}

def convertToLatex(c):
    if len(c) == 4: return c[0].upper()+"_{"+c[1:3]+"}^{"+c[3:4]+"}"
    else: 
        print "couldn't convert"
        return c   


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

files = [i for i in os.listdir(args.InputDir) if ".root" in i]
tL_chain = ROOT.TChain("tree")
tR_chain = ROOT.TChain("tree")
SM_chain = ROOT.TChain("tree")
for f in files:
    if "SM" in f: SM_chain.Add(args.InputDir + "/" + f)
    elif "cQ" in f: tL_chain.Add(args.InputDir + "/" + f)
    else: tR_chain.Add(args.InputDir + "/" + f)

branchnames = [i.GetName() for i in SM_chain.GetListOfBranches()]

nbins=25

for br in branchnames:
    h_SM = ROOT.TH1D("h_SM_"+br,";%s;a.u."%axis_name_dict[br][0],nbins,axis_name_dict[br][1],axis_name_dict[br][2])
    h_tL = ROOT.TH1D("h_tL_"+br,";%s;a.u."%axis_name_dict[br][0],nbins,axis_name_dict[br][1],axis_name_dict[br][2])
    h_tR = ROOT.TH1D("h_tR_"+br,";%s;a.u."%axis_name_dict[br][0],nbins,axis_name_dict[br][1],axis_name_dict[br][2])

    SM_chain.Draw("%s>>h_SM_%s"%(br,br))
    tL_chain.Draw("%s>>h_tL_%s"%(br,br))
    tR_chain.Draw("%s>>h_tR_%s"%(br,br))
    
    h_SM.Scale(1./h_SM.Integral())
    h_tL.Scale(1./h_tL.Integral())
    h_tR.Scale(1./h_tR.Integral())

    c = ROOT.TCanvas("c","c",600,600)
    c.SetMargin(0.10,0.05,0.15,0.05)
    h_SM.Draw("hist")
    h_tL.Draw("hist same")
    h_tR.Draw("hist same")

    #style
    h_SM.SetLineColor(2)
    h_SM.SetLineWidth(2)
    max = h_SM.GetBinContent(h_SM.GetMaximumBin())
    h_SM.GetYaxis().SetRangeUser(0,1.7*max)
    h_SM.GetYaxis().SetTickSize(0)
    h_SM.GetYaxis().SetLabelSize(0)
    h_SM.GetYaxis().SetTitleSize(0.06)
    h_SM.GetYaxis().SetTitleOffset(0.8)
    h_SM.GetXaxis().SetRangeUser(axis_name_dict[br][1],axis_name_dict[br][2])
    h_SM.GetXaxis().SetLabelSize(0.055)
    h_SM.GetXaxis().SetTitleSize(0.06)
    h_SM.GetXaxis().SetTitleOffset(1.1)
    h_SM.GetXaxis().SetNdivisions(505)

    h_tL.SetLineColor(3)
    h_tL.SetLineWidth(2)

    h_tR.SetLineColor(4)
    h_tR.SetLineWidth(2)


    l = ROOT.TLegend(0.12,0.8,0.93,0.93)
    l.SetNColumns(3)
    l.SetBorderSize(0)
    l.AddEntry(h_SM,"SM only","l")
    l.AddEntry(h_tL,"#splitline{EFT with}{left-handed top}","l")
    l.AddEntry(h_tR,"#splitline{EFT with}{right-handed top}","l")

    l.Draw("same")

    if not os.path.isdir(args.InputDir + "/PlotNetworkVariables"): os.mkdir(args.InputDir + "/PlotNetworkVariables")
    
    c.SaveAs("%s/PlotNetworkVariables/Plot_Network_%s.png"%(args.InputDir,br))
    c.SaveAs("%s/PlotNetworkVariables/Plot_Network_%s.pdf"%(args.InputDir,br))
