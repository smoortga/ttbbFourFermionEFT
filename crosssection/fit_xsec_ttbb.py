import os
from argparse import ArgumentParser
import sys
#import matplotlib.pyplot as plt
import numpy as np
from numpy import polyfit, diag, sqrt
from numpy.polynomial import Polynomial as P
from scipy.optimize import curve_fit
import ROOT
from array import array
import math
# plt.rc('text', usetex=True)
# plt.rc('font', family='sans-serif')

ROOT.gROOT.SetBatch(True)

parser = ArgumentParser()

parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttbb_DiLepton_ForValidationPerCouplings_fixed", help='path to the directory of all restricted processes')
#parser.add_argument('--SMxsecDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODEL_ttcc_dilepton_SMRestrictedMassAndCabibbo", help='path to the directory of SM-like xsec measurement')

args = parser.parse_args()




def extract_coupling(coupling,val):
    """ extract the numeric value of a coupling in the directory name"""
    if "minus" in val: return -float(val.replace(coupling+"_minus","").replace("p","."))
    else: return float(val.replace(coupling+"_","").replace("p","."))

#a = 0.0654444216985
def func(x, b, c):
        return 0.0654444216985*(1 + b*x + c*x*x)    
 
def convertToLatex(c):
    if len(c) == 4: return c[0].upper()+"_{"+c[1:3]+"}^{"+c[3:4]+"}"
    else: 
        print "couldn't convert"
        return c   

def makeFitFunction(p1,p2,units="pb"):
    text = ""
    if units == "pb":
        text = text + "#sigma [pb] = "
        text = text + "%.3f ("%0.0654444216985
        text = text + " 1 "
        if p1 > 0: text = text + "+ %.4f C"%math.fabs(p1)
        else: text = text + "- %.4f C"%math.fabs(p1)
        if p2 > 0: text = text + " + %.4f C^{2} )"%math.fabs(p2)
        else: text = text + "- %.4f C^{2} )"%math.fabs(p2)
    
    return text
        



process_dict = {}
subdirs = [d for d in os.listdir(args.InputDir) if os.path.isdir(args.InputDir+"/"+d) and not "fitted" in d and not 'Limits' in d and d[0] == "c"]
for subdir in subdirs:
    process_dict[subdir]= [dd for dd in os.listdir(args.InputDir+"/"+subdir) if os.path.isdir(args.InputDir+"/"+subdir+"/"+dd) and not "fitted_xsec" in dd]

#print process_dict

# coefficients with errors: sigma = p0 ( 1 + p1*C + p2*c^2 )
p0 = []
p0e = []
p1 = []
p1e = []
p2 = []
p2e = []
coupling_names = []

sm_xsec = 0.088 # pb from TOP-16-010
sm_xsec_error = 0.031# pb from TOP-16-010
sm_xsec_frac_error = 0.031/0.088
sm_xsec_frac_error_300fb = 0.1



# output files for limits
if not os.path.isdir(args.InputDir + "/fitted_xsec"): os.mkdir(args.InputDir + "/fitted_xsec")

limits_file = open('%s/fitted_xsec/limits_ttbb_xsec_CMS.txt'%(args.InputDir),"w")
limits_file.write("fractional error on SM measurement: %i %% \n"%int(100*sm_xsec_frac_error))

limits_file_300fb = open('%s/fitted_xsec/limits_ttbb_xsec_CMS_300fb.txt'%(args.InputDir),"w")
limits_file_300fb.write("fractional error on SM measurement: %i %% \n"%int(100*sm_xsec_frac_error_300fb))

for coupling,value in process_dict.iteritems():
    
    coupling_strengths = []
    xsec = []
    xsec_error = []
    nevents = []
    for v in value:
        coupling_strengths.append(extract_coupling(coupling,v))
        f_ = open(args.InputDir + "/" + coupling + "/" + v + "/cross_section.txt", 'r')
        lines = f_.readlines()
        xsec.append(float(lines[-1].split(" ")[2]))
        xsec_error.append(float(lines[-1].split(" ")[3]))
        nevents.append(float(lines[-1].split(" ")[4][:-1])) #need to remove \n character in the end
    
    print np.mean(np.asarray([i for idx, i in enumerate(xsec) if coupling_strengths[idx]==0.0]))
    
    coeff, covmat = curve_fit(func,coupling_strengths,xsec,sigma=xsec_error)#absolute_sigma
    errors = sqrt(diag(covmat))
    
    coupling_names.append(coupling)
    p0.append(coeff[0])
    p0e.append(errors[0])
    p1.append(coeff[1])
    p1e.append(errors[1])
    #p2.append(coeff[2])
    #p2e.append(errors[2])
    
    #**************************
    #
    #   PLOTTING
    #
    #**************************
    
    gr_data = ROOT.TGraphErrors(len(coupling_strengths),array('d',coupling_strengths),array('d',xsec),array('d',np.zeros(len(coupling_strengths))),array('d',xsec_error))
    gr_data.SetLineWidth(1)
    gr_data.SetMarkerStyle(20)
    gr_data.SetMarkerColor(4)
    xmin = -40
    xmax = +40
    nsteps = 200
    x = np.arange(xmin, xmax+ float(xmax-xmin)/float(nsteps), float(xmax-xmin)/float(nsteps))
    y = func(x,*(coeff))
    y_up = func(x,*(coeff + [i for i in errors]))
    y_down = func(x,*(coeff - [i for i in errors]))
    gr_fitted_xsec = ROOT.TGraphErrors(len(x),array('d',x),array('d',y),array('d',np.zeros(len(x))),array('d',[abs(i-j) for i,j in zip(y,y_up)]))
    gr_fitted_xsec.SetMarkerSize(0)
    gr_fitted_xsec.SetFillColor(2)
    gr_fitted_xsec.SetLineColor(2)
    gr_fitted_xsec.SetLineWidth(2)
    
    # CMS result @ 2.3 fb-1
    gr_CMS = ROOT.TGraphErrors(len(x),array('d',x),array('d',[sm_xsec]*len(x)),array('d',np.zeros(len(x))),array('d',[2*sm_xsec_frac_error*sm_xsec]*len(x)))
    gr_CMS.SetFillColorAlpha(41,0.7)
    gr_CMS.SetLineWidth(2)
    # CMS result @ 300 fb-1
    gr_CMS_300fb = ROOT.TGraphErrors(len(x),array('d',x),array('d',[min(y)]*len(x)),array('d',np.zeros(len(x))),array('d',[2*sm_xsec_frac_error_300fb*min(y)]*len(x)))
    gr_CMS_300fb.SetFillColorAlpha(49,0.8)
    gr_CMS_300fb.SetLineWidth(2)
    gr_CMS_300fb.SetLineStyle(2)
    
    mg = ROOT.TMultiGraph("mg",";;cross section [pb]")
    mg.Add(gr_CMS,"3")
    mg.Add(gr_CMS_300fb,"3")
    mg.Add(gr_fitted_xsec,"c3")
    mg.Add(gr_data,"pe")
    
    c = ROOT.TCanvas("c","c",800,700)
    uppad = ROOT.TPad("u","u",0.,0.4,1.,1.)
    downpad = ROOT.TPad("d","d",0.,0.,1.,0.4)
    uppad.Draw()
    downpad.Draw()
    uppad.cd()
    ROOT.gPad.SetMargin(0.15,0.05,0.01,0.1)
    #gr_data.Draw("AP")
    mg.Draw("A")
    ymin = sm_xsec - 5*sm_xsec_frac_error*sm_xsec
    ymax = sm_xsec + 10*sm_xsec_frac_error*sm_xsec
    mg.SetMinimum(ymin)
    mg.SetMaximum(ymax)
    mg.GetYaxis().CenterTitle()
    mg.GetYaxis().SetTitleSize(0.09)
    mg.GetYaxis().SetTitleOffset(0.7)
    mg.GetYaxis().SetLabelSize(0.06)
    mg.GetXaxis().SetTitleSize(0.0)
    mg.GetXaxis().SetTitleOffset(1.05)
    mg.GetXaxis().SetLabelSize(0.0)
    mg.GetXaxis().SetRangeUser(xmin+1,xmax-1)
    
    
    l = ROOT.TLegend(0.35,0.6,0.94,0.88)
    l.SetBorderSize(1)
    l.AddEntry(gr_data,"sample points","pe1")
    l.AddEntry(gr_fitted_xsec,"fitted cross section","l")
    l.AddEntry(gr_CMS,"CMS result 95% CL (2.3 fb^{-1})","fl")
    l.AddEntry(gr_CMS_300fb,"CMS prospect 95% CL (300 fb^{-1})","fl")
    
    l.Draw("same")
    
    t = ROOT.TLatex()
    t.SetTextSize(0.07)
    t.SetTextFont(42)
    t.SetTextColor(2)
    t.DrawLatexNDC(0.15,0.93,makeFitFunction(coeff[0],coeff[1]))
    
    
    downpad.cd()
    ROOT.gPad.SetMargin(0.15,0.05,0.25,0.01)
    y_1 = [(i-sm_xsec)**2/((sm_xsec_frac_error*sm_xsec)**2) for i in y]
    ymin = -1
    ymax = 13
    gr_chi2 = ROOT.TGraph(len(x),array('d',x),array('d',y_1))
    gr_chi2.SetLineWidth(2)
    gr_chi2.SetLineColor(1)
    gr_chi2.SetLineStyle(1)
    # prospects @ 300 fb
    y_2 = [(i-min(y))**2/((sm_xsec_frac_error_300fb*min(y))**2) for i in y]
    gr_chi2_300fb = ROOT.TGraph(len(x),array('d',x),array('d',y_2))
    gr_chi2_300fb.SetLineWidth(2)
    gr_chi2_300fb.SetLineColor(1)
    gr_chi2_300fb.SetLineStyle(2)
    
    
    # get roots
    p_1 = P.fit(x, y_1, 4)
    roots_1 = sorted([i.real for i in (p_1 - 3.84).roots() if i.imag == 0])
    print roots_1
    while len(roots_1)>2: roots_1 = roots_1[1:-1]
    x_limit_1 = np.arange(roots_1[0], roots_1[1]+abs(roots_1[1]-roots_1[0])/200., abs(roots_1[1]-roots_1[0])/200.)
    gr_CMS_limit = ROOT.TGraphErrors(len(x_limit_1),array('d',x_limit_1),array('d',[5]*len(x_limit_1)),array('d',np.zeros(len(x_limit_1))),array('d',[5000]*len(x_limit_1)))
    gr_CMS_limit.SetFillColorAlpha(41,0.7)
    gr_CMS_limit.SetLineWidth(2)
    
    p_2 = P.fit(x, y_2, 4)
    roots_2 = sorted([i.real for i in (p_2 - 3.84).roots() if i.imag == 0])
    print roots_2
    while len(roots_2)>2: roots_2 = roots_2[1:-1]
    x_limit_2 = np.arange(roots_2[0], roots_2[1]+abs(roots_2[1]-roots_2[0])/200., abs(roots_2[1]-roots_2[0])/200.)
    gr_CMS_limit_300fb = ROOT.TGraphErrors(len(x_limit_2),array('d',x_limit_2),array('d',[5]*len(x_limit_2)),array('d',np.zeros(len(x_limit_2))),array('d',[5000]*len(x_limit_2)))
    gr_CMS_limit_300fb.SetFillColorAlpha(49,0.8)
    #gr_CMS_limit_300fb.SetFillStyle(3244)
    gr_CMS_limit_300fb.SetLineWidth(2)
    
    
    mg2 = ROOT.TMultiGraph("mg2",";%s [TeV^{-2}];#chi^{2}"%convertToLatex(coupling))
    mg2.Add(gr_CMS_limit,"3")
    mg2.Add(gr_CMS_limit_300fb,"3")
    mg2.Add(gr_chi2,"l")
    mg2.Add(gr_chi2_300fb,"l")
    
    mg2.Draw("A")
    mg2.SetMinimum(ymin)
    mg2.SetMaximum(ymax)
    mg2.GetYaxis().CenterTitle()
    mg2.GetYaxis().SetTitleSize(0.12)
    mg2.GetYaxis().SetTitleOffset(0.55)
    mg2.GetYaxis().SetLabelSize(0.09)
    mg2.GetYaxis().SetNdivisions(8)
    mg2.GetXaxis().CenterTitle()
    mg2.GetXaxis().SetTitleSize(0.11)
    mg2.GetXaxis().SetTitleOffset(1.)
    mg2.GetXaxis().SetLabelSize(0.09)
    mg2.GetXaxis().SetRangeUser(xmin+1,xmax-1)
    
    # line
    line = ROOT.TLine(xmin+1,3.84,xmax-1,3.84)
    line.SetLineColor(13)
    line.SetLineWidth(1)
    line.SetLineStyle(2)
    line.Draw("same")
    
    
    limits_file.write("%s & [%.2f,%.2f] \n"%(coupling,roots_1[0],roots_1[1]))
    limits_file_300fb.write("%s & [%.2f,%.2f] \n"%(coupling,roots_2[0],roots_2[1]))
    
    c.SaveAs('%s/fitted_xsec/fit_xsec_%s.png'%(args.InputDir,coupling))
    c.SaveAs('%s/fitted_xsec/fit_xsec_%s.pdf'%(args.InputDir,coupling))
    
    print "done, saved in %s/fitted_xsec/fit_xsec_%s.png"%(args.InputDir,coupling)
    
limits_file.close()
limits_file_300fb.close()    

#change order
indices = [7,5,0,4,6,1,3,2]
p0 = np.asarray(p0)[indices]
p0e = np.asarray(p0e)[indices]
p1 = np.asarray(p1)[indices]
p1e = np.asarray(p1e)[indices]
# p2 = np.asarray(p2)[indices]
# p2e = np.asarray(p2e)[indices]
coupling_names = np.asarray(coupling_names)[indices]


c = ROOT.TCanvas("c","c",1100,600)
#uppad = ROOT.TPad("u","u",0.,1.,1.,1.)
midpad = ROOT.TPad("d","d",0.,0.57,1.,1.)
downpad = ROOT.TPad("d","d",0.,0.0,1.,0.57)
#uppad.Draw()
midpad.Draw()
downpad.Draw()

# uppad.cd()
# ROOT.gPad.SetMargin(0.15,0.05,0.0,0.05)
# ROOT.gPad.SetGrid(0,1)
x_ = range(len(coupling_names))
# print x_
# gr_p00 = ROOT.TGraphErrors(len(x_),array('d',x_),array('d',np.asarray([0.0654444216985]*len(x_))),array('d',np.zeros(len(x_))),array('d',np.ones(len(x_))))
# gr_p00.SetMarkerStyle(20)
# gr_p00.SetMarkerSize(1.6)
# gr_p00.SetMarkerColor(4)
# gr_p00.SetLineWidth(2)
# gr_p00.SetLineColor(4)
# mg_p00 = ROOT.TMultiGraph("mg_p00",";;p_{0} [pb]")
# mg_p00.Add(gr_p00)
# mg_p00.Draw("Ap")
# mg_p00.SetTitle("")
# mg_p00.GetXaxis().SetLabelSize(0)
# mg_p00.GetXaxis().SetTickSize(0)
# mg_p00.GetYaxis().CenterTitle(1)
# mg_p00.GetYaxis().SetLabelSize(0.115)
# mg_p00.GetYaxis().SetTickSize(0.02)
# mg_p00.GetYaxis().SetTitle("p_{0} [pb]")
# mg_p00.GetYaxis().SetTitleSize(0.17)
# mg_p00.GetYaxis().SetTitleOffset(0.42)
# mg_p00.SetMinimum(0.0625)
# mg_p00.SetMaximum(0.0675)
# mg_p00.GetYaxis().SetNdivisions(6)

midpad.cd()
ROOT.gPad.SetMargin(0.15,0.05,0.0,0.1)
ROOT.gPad.SetGrid(0,1)
p0 = [i*100 for i in p0]
p0e = [i*100 for i in p0e]
gr_p0 = ROOT.TGraphErrors(len(x_),array('d',x_),array('d',p0),array('d',np.zeros(len(x_))),array('d',p0e))
gr_p0.SetMarkerStyle(20)
gr_p0.SetMarkerSize(1.6)
gr_p0.SetMarkerColor(4)
gr_p0.SetLineWidth(2)
gr_p0.SetLineColor(4)
mg_p0 = ROOT.TMultiGraph("mg_p0",";;p_{1} [pb]")
mg_p0.Add(gr_p0)
mg_p0.Draw("Ap")
mg_p0.SetTitle("")
mg_p0.GetXaxis().SetLabelSize(0)
mg_p0.GetXaxis().SetTickSize(0)
mg_p0.GetYaxis().CenterTitle(1)
mg_p0.GetYaxis().SetLabelSize(0.13)
mg_p0.GetYaxis().SetTickSize(0.02)
mg_p0.GetYaxis().SetTitle("p_{1} [%]")
mg_p0.GetYaxis().SetTitleSize(0.2)
mg_p0.GetYaxis().SetTitleOffset(0.35)
mg_p0.SetMinimum(-0.15)
mg_p0.SetMaximum(0.55)
mg_p0.GetYaxis().SetNdivisions(6)

downpad.cd()
ROOT.gPad.SetMargin(0.15,0.05,0.3,0.0)
ROOT.gPad.SetGrid(0,1)
# Dummy histogram for x labels
# h = ROOT.TH1D("h","",len(x_)+1,-0.5,x[-1]+0.5)
# h.GetXaxis().SetTickSize(0.02)
# h.GetXaxis().SetLabelSize(0.1)
# h.Draw("hist")
p1 = [i*100 for i in p1]
p1e = [i*100 for i in p1e]
gr_p1 = ROOT.TGraphErrors(len(x_),array('d',x_),array('d',p1),array('d',np.zeros(len(x_))),array('d',p1e))
gr_p1.SetMarkerStyle(20)
gr_p1.SetMarkerSize(1.6)
gr_p1.SetMarkerColor(4)
gr_p1.SetLineWidth(2)
gr_p1.SetLineColor(4)
mg_p1 = ROOT.TMultiGraph("mg_p1",";;p_{2} [pb]")
mg_p1.Add(gr_p1)
mg_p1.Draw("Ap same")
mg_p1.SetTitle("")
mg_p1.GetXaxis().SetLabelSize(0.2)
mg_p1.GetXaxis().SetTickSize(0)
mg_p1.GetYaxis().CenterTitle(1)
mg_p1.GetYaxis().SetLabelSize(0.1)
mg_p1.GetYaxis().SetTickSize(0.02)
mg_p1.GetYaxis().SetTitle("p_{2} [%]")
mg_p1.GetYaxis().SetTitleSize(0.15)
mg_p1.GetYaxis().SetTitleOffset(0.45)
mg_p1.SetMinimum(-0.05)
mg_p1.SetMaximum(0.68)
mg_p1.GetYaxis().SetNdivisions(6)
for idx, i in enumerate(coupling_names):
    bin_index = mg_p1.GetXaxis().FindBin(idx)
    mg_p1.GetXaxis().SetBinLabel(bin_index,convertToLatex(i))
mg_p1.GetXaxis().LabelsOption("h")
mg_p1.GetXaxis().SetLabelOffset(0.02)



c.SaveAs('%s/fitted_xsec/fit_xsec_coeff_results.png'%(args.InputDir))
c.SaveAs('%s/fitted_xsec/fit_xsec_coeff_results.pdf'%(args.InputDir))

# f, axarr = plt.subplots(3, sharex=True, figsize=(12,8))
# x_ = range(len(coupling_names))
# latex_coupling_names = [r'$%s$'%convertToLatex(c) for c in coupling_names]
# plt.xticks(x_, latex_coupling_names)
# axarr[0].errorbar(x_, p0, yerr = p0e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
# axarr[0].set_ylabel("p0 [pb]", fontsize = 15)
# axarr[0].grid(True)
# axarr[0].set_xlim(x_[0]-1, x_[-1]+1)
# axarr[0].set_ylim(sm_xsec - 3*sm_xsec_error, sm_xsec + 3*sm_xsec_error)
# textymin, textymax = axarr[0].get_ylim()
# axarr[0].text(x_[0]-1,textymax + (textymax - textymin)/float(10), r'$\sigma=p0 ( 1 + p1 \ C + p2 \ C^2 )$', fontsize=20)
# axarr[0].set_xlim(x_[0]-1, x_[-1]+1)
# axarr[0].fill_between([-1] + x_ + [len(x_)+1], sm_xsec - (sm_xsec_frac_error*sm_xsec), sm_xsec + (sm_xsec_frac_error*sm_xsec), facecolor='red', alpha=0.3, label='SM cross section')
# axarr[0].legend(loc='lower right')
# axarr[1].errorbar(x_, p1, yerr = p1e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
# axarr[1].set_ylabel("p1", fontsize = 15)
# axarr[1].grid(True)
# axarr[2].errorbar(x_, p2, yerr = p2e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
# axarr[2].set_ylabel("p2", fontsize = 15)
# axarr[2].set_xlabel("non-zero EFT coupling", fontsize = 15)
# axarr[2].grid(True)
# f.savefig('%s/fitted_xsec/fit_xsec_coeff_results.png'%(args.InputDir))
# f.savefig('%s/fitted_xsec/fit_xsec_coeff_results.pdf'%(args.InputDir))