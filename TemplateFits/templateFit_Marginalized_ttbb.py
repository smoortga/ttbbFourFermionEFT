import ROOT
import os
from argparse import ArgumentParser
import sys
from math import sqrt
from numpy import polyfit, diag, sqrt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import Polynomial as P
from numpy.random import normal,poisson
import math
from array import array
import pickle


"""
make sure to do the following before running this script!
cdir=$(pwd)
export LD_LIBRARY_PATH=${cdir}:${cdir}/Delphes:$LD_LIBRARY_PATH
"""

def extract_xsec(dir,c,v):
    if os.path.isfile(dir + "/" + c + "/" + c.split("_")[0] + "_" + v.split("_")[0]+ "_" + c.split("_")[1] + "_" + v.split("_")[1]  + "/cross_section.txt"):
        f_ = open(dir + "/" + c + "/" + c.split("_")[0] + "_" + v.split("_")[0]+ "_" + c.split("_")[1] + "_" + v.split("_")[1]  + "/cross_section.txt", 'r')
        lines = f_.readlines()
        xsec = float(lines[-1].split(" ")[2])
        return xsec
    else:
        print "ERROR: could not extract cross section from "+dir + "/" + c + "/" + c.split("_")[0] + "_" + v.split("_")[0]+ "_" + c.split("_")[1] + "_" + v.split("_")[1]  + "/cross_section.txt"
        sys.exit(1)

def get_Original_nevents(dir,c,v):
    events_dir = dir + "/" + c + "/" + c.split("_")[0] + "_" + v.split("_")[0]+ "_" + c.split("_")[1] + "_" + v.split("_")[1]  + "/Events"
    if os.path.isdir(events_dir):
        t_ = ROOT.TChain("Delphes")
        runs = [i for i in os.listdir(events_dir) if "run_" in i]
        for r in runs:
            run_dir = events_dir + "/" + r
            t_.Add(run_dir+ "/tag_1_delphes_events.root")
        return t_.GetEntries()
        
    else:
        print "ERROR: could not retrieve original number of events from "+ dir + "/" + c + "/" + c.split("_")[0] + "_" + v.split("_")[0]+ "_" + c.split("_")[1] + "_" + v.split("_")[1]  + "/Events"
        sys.exit(1)



def convert_coupling_toNumerical(v):
    """ extract the numeric value of a coupling in the directory name"""
    if "minus" in v: return -float(v.replace("minus","").replace("p","."))
    else: return float(v.replace("p","."))
    
def convertToLatex(c):
    if len(c) == 4: return c[0].upper()+"_{"+c[1:3]+"}^{"+c[3:4]+"}"
    else: 
        print "couldn't convert"
        return c   

        
def func2(x, a, b, c):
        return (a + b*x + c*x*x)   

def func4(x, a, b, c, d, e):
        return (a + b*x + c*x*x + d*x*x*x + e*x*x*x*x)   

def func(x, p0, p1a, p1b, p2aa, p2bb, p2ab):
        return p0 + p1a*x[0] + p1b*x[1] + p2aa*x[0]*x[0] + p2bb*x[1]*x[1] + p2ab*x[0]*x[1]   


def fitTemplate(data_hist,SM_templ,tL_templ,tR_templ,coupl_names, coupl_values, savedir, xsec, original_nevents, int_lumi = 300 , xaxis_title = "P(t_{L}) + P(t_{R})", yaxis_title = "#frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",npseudo=200):
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.kError)
    # d = ROOT.Long(0)
#     print d
#     ROOT.gMinuit.mnexcm("SET NOWarnings",array('d',[0]),0,d)
    #ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 2001;")
    #m = ROOT.RooMinuit()
    #m.setPrintLevel(-1)
    
    # scale the data histogram to the 300 fb-1 level
    print data_hist.Integral()
    print xsec, int_lumi, original_nevents, float(xsec*int_lumi)/float(original_nevents)
    data_hist.Scale(float(xsec*int_lumi)/float(original_nevents))
    print data_hist.Integral()
    sm_guess = 2200
    eft_guess = abs(data_hist.Integral() - sm_guess)
    tl_val = convert_coupling_toNumerical(coupl_values.split("_")[0])
    tR_val = convert_coupling_toNumerical(coupl_values.split("_")[1])
    if not (tl_val == 0 and tR_val == 0):
        tl_guess = max(eft_guess*(tl_val**2/(tl_val**2+tR_val**2)),1)
        tR_guess = max(eft_guess*(tR_val**2/(tl_val**2+tR_val**2)),1)
    else: 
        tl_guess=0
        tR_guess=0
    print sm_guess,eft_guess,tl_guess,tR_guess
    
    #########################
    #
    # Add uncertainty to the data
    #
    #########################
    fitted_tL_arr = []
    fitted_tR_arr = []
    fitted_SM_arr = []
    fitted_sigma_tL_arr = []
    fitted_sigma_tR_arr = []
    fitted_sigma_SM_arr = []
    chi2_ndof_arr = []
    for i in range(npseudo):
        pseudo_data_hist = data_hist.Clone()
        #frac_error = SM_unc 
        for binx in range(pseudo_data_hist.GetNbinsX()):
            for biny in range(pseudo_data_hist.GetNbinsY()):
                content = pseudo_data_hist.GetBinContent(binx+1,biny+1,)
                if content != 0:
                    random_var = poisson(int(content))
                    #print content, random_var
                    pseudo_data_hist.SetBinContent(binx+1,biny+1,random_var)
                    if random_var > 0: pseudo_data_hist.SetBinError(binx+1,biny+1,sqrt(random_var))
                    else: pseudo_data_hist.SetBinError(binx+1,biny+1,1)
                    #print content, random_var
                else:
                    pseudo_data_hist.SetBinContent(binx+1,biny+1,0)
                    pseudo_data_hist.SetBinError(binx+1,biny+1,1)
        
    
        #print SM_templ.GetXaxis().GetBinLowEdge(1), SM_templ.GetXaxis().GetBinUpEdge(SM_templ.GetNbinsX())
        varx = ROOT.RooRealVar("vartx","vartx",SM_templ.GetXaxis().GetBinLowEdge(1),SM_templ.GetXaxis().GetBinUpEdge(SM_templ.GetNbinsX()))
        vary = ROOT.RooRealVar("varty","varty",SM_templ.GetYaxis().GetBinLowEdge(1),SM_templ.GetYaxis().GetBinUpEdge(SM_templ.GetNbinsY()))
        SM_Hist = ROOT.RooDataHist("SM_Hist","SM_Hist",ROOT.RooArgList(varx,vary),SM_templ)
        SM_PDF = ROOT.RooHistPdf("SM_PDF","SM_PDF",ROOT.RooArgSet(varx,vary),SM_Hist)
        tL_Hist = ROOT.RooDataHist("tL_Hist","tL_Hist",ROOT.RooArgList(varx,vary),tL_templ)
        tL_PDF = ROOT.RooHistPdf("tL_PDF","tL_PDF",ROOT.RooArgSet(varx,vary),tL_Hist)
        tR_Hist = ROOT.RooDataHist("tR_Hist","tR_Hist",ROOT.RooArgList(varx,vary),tR_templ)
        tR_PDF = ROOT.RooHistPdf("tR_PDF","tR_PDF",ROOT.RooArgSet(varx,vary),tR_Hist)
        data_Hist = ROOT.RooDataHist("pseudo_data_hist","pseudo_data_hist",ROOT.RooArgList(varx,vary),pseudo_data_hist)
        data_PDF = ROOT.RooHistPdf("data_PDF","data_PDF",ROOT.RooArgSet(varx,vary),data_Hist)

    
        N_SM = ROOT.RooRealVar("N_SM","N_SM",int(sm_guess),-(pseudo_data_hist.Integral()/100.),2*(pseudo_data_hist.Integral())) 
        N_tL = ROOT.RooRealVar("N_tL","N_tL",int(tl_guess),-(pseudo_data_hist.Integral()/100.),2*(pseudo_data_hist.Integral()))
        N_tR = ROOT.RooRealVar("N_tR","N_tR",int(tR_guess),-(pseudo_data_hist.Integral()/100.),2*(pseudo_data_hist.Integral()))
        # N_SM = ROOT.RooRealVar("N_SM","N_SM",int(pseudo_data_hist.Integral()/2.),0,2*(pseudo_data_hist.Integral())) 
#         N_tL = ROOT.RooRealVar("N_tL","N_tL",int(pseudo_data_hist.Integral()/4.),0,2*(pseudo_data_hist.Integral()))
#         N_tR = ROOT.RooRealVar("N_tR","N_tR",int(pseudo_data_hist.Integral()/8.),0,2*(pseudo_data_hist.Integral()))

        
        mod = ROOT.RooAddPdf("mod","mod",ROOT.RooArgList(SM_PDF,tL_PDF, tR_PDF),ROOT.RooArgList(N_SM,N_tL,N_tR)); 
        
        fitRes = mod.fitTo(data_Hist,ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1000), ROOT.RooFit.Verbose(False), ROOT.RooFit.Extended(True), ROOT.RooFit.Save())
        
        chi2_var = ROOT.RooChi2Var("chi2_var", "chi2_var", mod, data_Hist)
        chi2 = chi2_var.getVal()
        
        fitted_tL = N_tL.getVal()
        sigma_fitted_tL = N_tL.getError()
        fitted_tR = N_tR.getVal()
        sigma_fitted_tR = N_tR.getError()
        fitted_SM = N_SM.getVal()
        sigma_fitted_SM = N_SM.getError()
        
        # print "chi2/ndof:",chi2/99.
#         print "SM:",fitted_SM,"+-", sigma_fitted_SM
#         print "tL:",fitted_tL,"+-", sigma_fitted_tL
#         print "tR:",fitted_tR,"+-", sigma_fitted_tR
        
        fitted_tL_arr.append(fitted_tL)
        fitted_tR_arr.append(fitted_tR)
        fitted_SM_arr.append(fitted_SM)
        fitted_sigma_tL_arr.append(sigma_fitted_tL)
        fitted_sigma_tR_arr.append(sigma_fitted_tR)
        fitted_sigma_SM_arr.append(sigma_fitted_SM)
        chi2_ndof_arr.append(chi2/(pseudo_data_hist.GetNbinsX()*pseudo_data_hist.GetNbinsY() -1))
        
    
    #print pseudo_data_hist.GetEntries()
    # mean_fitted_tL = np.mean(np.asarray(fitted_tL_arr))
#     err_fitted_tL = np.std(np.asarray(fitted_tL_arr))
#     mean_fitted_tR = np.mean(np.asarray(fitted_tR_arr))
#     err_fitted_tR = np.std(np.asarray(fitted_tR_arr))
#     mean_fitted_SM = np.mean(np.asarray(fitted_SM_arr))
#     err_fitted_SM = np.std(np.asarray(fitted_SM_arr))
    mean_fitted_tL = np.average(np.asarray(fitted_tL_arr), weights = np.asarray([1./i for i in fitted_sigma_tL_arr]))
    #err_fitted_tL = np.std(np.asarray(fitted_tL_arr))#np.average([i-mean_fitted_tL for i in np.asarray(fitted_tL_arr)], weights = np.asarray([1./i for i in fitted_sigma_tL_arr]))
    err_fitted_tL = np.mean(np.asarray(fitted_sigma_tL_arr))
    mean_fitted_tR = np.average(np.asarray(fitted_tR_arr), weights = np.asarray([1./i for i in fitted_sigma_tR_arr]))
    #err_fitted_tR = np.std(np.asarray(fitted_tR_arr))
    err_fitted_tR = np.mean(np.asarray(fitted_sigma_tR_arr))
    mean_fitted_SM = np.average(np.asarray(fitted_SM_arr), weights = np.asarray([1./i for i in fitted_sigma_SM_arr]))
    #err_fitted_SM = np.std(np.asarray(fitted_SM_arr))
    err_fitted_SM = np.mean(np.asarray(fitted_sigma_SM_arr))
    mean_chi2 = np.mean(np.asarray(chi2_ndof_arr))
    std_chi2 = np.std(np.asarray(chi2_ndof_arr))
    print "*************RESULT: %s******************"%coupl_values
    print "SM:",mean_fitted_SM,"+-", err_fitted_SM
    print "tL:", mean_fitted_tL,"+-",err_fitted_tL
    print "tR:", mean_fitted_tR,"+-",err_fitted_tR
    print "chi2/ndof:", mean_chi2,"+-",std_chi2
    print "*********************************************"
    
    
    #########################
    #
    # Plot the postfitted histos/teplates
    #
    #########################
    tL_templ_copy = tL_templ.Clone()
    tR_templ_copy = tR_templ.Clone()
    SM_templ_copy = SM_templ.Clone()
    c_postfit = ROOT.TCanvas("c_postfit","c_postfit",1000,500)
    c_postfit.Divide(2,1)

    tL_templ_copy.Scale(mean_fitted_tL/tL_templ_copy.Integral())
    tR_templ_copy.Scale(mean_fitted_tR/tR_templ_copy.Integral())
    SM_templ_copy.Scale(mean_fitted_SM/SM_templ_copy.Integral())
   
    c_postfit.cd(1)
    ROOT.gPad.SetMargin(0.15,0.02,0.2,0.1)
    tL_templ_copy_ProjectionX = tL_templ_copy.ProjectionX()
    tR_templ_copy_ProjectionX = tR_templ_copy.ProjectionX()
    SM_templ_copy_ProjectionX = SM_templ_copy.ProjectionX()
    tL_templ_copy_ProjectionX.SetFillColor(3)
    tL_templ_copy_ProjectionX.SetLineWidth(0)
    tR_templ_copy_ProjectionX.SetFillColor(ROOT.kBlue-7)
    tR_templ_copy_ProjectionX.SetLineWidth(0)
    SM_templ_copy_ProjectionX.SetFillColor(2)
    SM_templ_copy_ProjectionX.SetLineWidth(0)
    stack_postfit_X = ROOT.THStack("stack_postfit_X","")
    stack_postfit_X.Add(SM_templ_copy_ProjectionX)
    stack_postfit_X.Add(tL_templ_copy_ProjectionX)
    stack_postfit_X.Add(tR_templ_copy_ProjectionX)
    stack_postfit_X.SetMinimum(0)
    stack_postfit_X.SetMaximum(2.5*stack_postfit_X.GetMaximum())
    stack_postfit_X.Draw("hist")
    stack_postfit_X.GetXaxis().SetTitle(xaxis_title)
    stack_postfit_X.GetXaxis().SetTitleOffset(1.55)
    stack_postfit_X.GetXaxis().SetTitleSize(0.05)
    stack_postfit_X.GetXaxis().SetLabelSize(0.045)
    stack_postfit_X.GetYaxis().SetTitle("Entries")
    stack_postfit_X.GetYaxis().SetTitleOffset(1.65)
    stack_postfit_X.GetYaxis().SetTitleSize(0.05)
    stack_postfit_X.GetYaxis().SetLabelSize(0.045)
    error_hist_X = stack_postfit_X.GetHistogram().Clone()
    error_hist_X.SetFillStyle(3354)
    error_hist_X.SetFillColor(1)
    error_hist_X.SetLineWidth(0)
    for binx in range(error_hist_X.GetNbinsX()):
        error_hist_X.SetBinContent(binx+1,SM_templ_copy_ProjectionX.GetBinContent(binx+1)+tL_templ_copy_ProjectionX.GetBinContent(binx+1)+tR_templ_copy_ProjectionX.GetBinContent(binx+1))
        error_hist_X.SetBinError(binx+1,error_hist_X.GetBinContent(binx+1)*(err_fitted_SM+err_fitted_tL+err_fitted_tR)/float(mean_fitted_SM+mean_fitted_tL+mean_fitted_tR))
    error_hist_X.Draw("same E2")
    pseudo_data_hist_X = pseudo_data_hist.ProjectionX()
    pseudo_data_hist_X.SetLineWidth(2)
    pseudo_data_hist_X.SetMarkerStyle(20)
    pseudo_data_hist_X.Draw("same pX1E1")
    l_postfit = ROOT.TLegend(0.2,0.6,0.95,0.89)
    l_postfit.SetBorderSize(0)
    l_postfit.SetTextSize(0.05)
    l_postfit.SetHeader("#bf{postfit}")
    l_postfit.SetNColumns(2)
    l_postfit.AddEntry(tR_templ_copy_ProjectionX,"t_{R} template","f")
    l_postfit.AddEntry(pseudo_data_hist_X,"pseudo-data","pe")
    l_postfit.AddEntry(tL_templ_copy_ProjectionX,"t_{L} template","f")
    l_postfit.AddEntry(error_hist_X,"fit uncertainty","f")
    l_postfit.AddEntry(SM_templ_copy_ProjectionX,"SM template","f")
    
    l_postfit.Draw("same")
    text = ROOT.TLatex()
    text.SetTextSize(0.055)
    text.SetTextFont(42)
    text.DrawLatexNDC(0.15,0.93,"%s = %.0f TeV^{-2}        %s = %.0f TeV^{-2}"%(convertToLatex(coupl_names.split("_")[0]),convert_coupling_toNumerical(coupl_values.split("_")[0]),convertToLatex(coupl_names.split("_")[1]),convert_coupling_toNumerical(coupl_values.split("_")[1])))
    
    c_postfit.cd(2)
    ROOT.gPad.SetMargin(0.15,0.02,0.2,0.1)
    tL_templ_copy_ProjectionY = tL_templ_copy.ProjectionY()
    tR_templ_copy_ProjectionY = tR_templ_copy.ProjectionY()
    SM_templ_copy_ProjectionY = SM_templ_copy.ProjectionY()
    tL_templ_copy_ProjectionY.SetFillColor(3)
    tL_templ_copy_ProjectionY.SetLineWidth(0)
    tR_templ_copy_ProjectionY.SetFillColor(ROOT.kBlue-7)
    tR_templ_copy_ProjectionY.SetLineWidth(0)
    SM_templ_copy_ProjectionY.SetFillColor(2)
    SM_templ_copy_ProjectionY.SetLineWidth(0)
    stack_postfit_Y = ROOT.THStack("stack_postfit_Y","")
    stack_postfit_Y.Add(SM_templ_copy_ProjectionY)
    stack_postfit_Y.Add(tL_templ_copy_ProjectionY)
    stack_postfit_Y.Add(tR_templ_copy_ProjectionY)
    stack_postfit_Y.SetMinimum(0)
    stack_postfit_Y.SetMaximum(2.5*stack_postfit_Y.GetMaximum())
    stack_postfit_Y.Draw("hist")
    stack_postfit_Y.GetXaxis().SetTitle(yaxis_title)
    stack_postfit_Y.GetXaxis().SetTitleOffset(1.55)
    stack_postfit_Y.GetXaxis().SetTitleSize(0.05)
    stack_postfit_Y.GetXaxis().SetLabelSize(0.045)
    stack_postfit_Y.GetYaxis().SetTitle("Entries")
    stack_postfit_Y.GetYaxis().SetTitleOffset(1.65)
    stack_postfit_Y.GetYaxis().SetTitleSize(0.05)
    stack_postfit_Y.GetYaxis().SetLabelSize(0.045)
    error_hist_Y = stack_postfit_Y.GetHistogram().Clone()
    error_hist_Y.SetFillStyle(3354)
    error_hist_Y.SetFillColor(1)
    for binx in range(error_hist_Y.GetNbinsX()):
        error_hist_Y.SetBinContent(binx+1,SM_templ_copy_ProjectionY.GetBinContent(binx+1)+tL_templ_copy_ProjectionY.GetBinContent(binx+1)+tR_templ_copy_ProjectionY.GetBinContent(binx+1))
        error_hist_Y.SetBinError(binx+1,error_hist_Y.GetBinContent(binx+1)*(err_fitted_SM+err_fitted_tL+err_fitted_tR)/float(mean_fitted_SM+mean_fitted_tL+mean_fitted_tR))
    error_hist_Y.Draw("same E2")
    pseudo_data_hist_Y = pseudo_data_hist.ProjectionY()
    pseudo_data_hist_Y.SetLineWidth(2)
    pseudo_data_hist_Y.SetMarkerStyle(20)
    pseudo_data_hist_Y.Draw("same pX1E1")
    l_postfit.Draw("same")
    
    c_postfit.SaveAs("%s/postfit_%s_%s_%s_%s.pdf"%(savedir,coupl_names.split("_")[0], coupl_values.split("_")[0],coupl_names.split("_")[1], coupl_values.split("_")[1]))
    print "Saved fitted output histograms in %s/postfit_%s_%s_%s_%s.pdf"%(savedir,coupl_names.split("_")[0], coupl_values.split("_")[0],coupl_names.split("_")[1], coupl_values.split("_")[1])
    
    return mean_fitted_SM,err_fitted_SM , mean_fitted_tL, err_fitted_tL, mean_fitted_tR, err_fitted_tR
    
    

def main():
    ROOT.gROOT.SetBatch(True)
    
    parser = ArgumentParser()
    parser.add_argument('--TemplateFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Training_fixed/templates/templates.root", help='root file with templates')
    parser.add_argument('--ValidationFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Marginalized/templates/validation_data.root", help='root file with validation data')
    parser.add_argument('--OutputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Marginalized/templates/control_plots/", help='path to the converted delphes training files')
    parser.add_argument('--XsecDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttbb_DiLepton_ForValidationPerCouplings_Marginalized/", help='path to MODELSCAN directory that holds the cross sections')
    parser.add_argument('--ValidationDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Marginalized/", help='path to MODELSCAN directory that holds the cross sections')
    args = parser.parse_args()
    
    if not os.path.isdir(args.OutputDir): os.mkdir(args.OutputDir)
    
    f_templ = ROOT.TFile(args.TemplateFile)
    f_data = ROOT.TFile(args.ValidationFile)
    
    fitting_dict = { #name:class_number
            'cQb1_ctb1':  {}
    }
    
#     fracerr_measured_sm=0.2
#     limits_file = open('%s/limits_ttbb_templateFit.txt'%(args.XsecDir),"w")
#     #limits_file = open('%s/limits_nocut.txt'%(args.InputDir),"w")
#     limits_file.write("fractional error on SM measurement: %i %% \n"%int(100*fracerr_measured_sm))

    
    for c,result_dict in fitting_dict.iteritems():
        print "Processing: %s"%c

        
        # loop over histograms in validation file
        print "*** Retrieving Histograms for pseudodata ***"
        histograms = []
        for key in f_data.GetListOfKeys():
            hist = key.ReadObj()
            if c.split("_")[0] in hist.GetName() and c.split("_")[1] in hist.GetName(): histograms.append(hist)  
        
        
        #extract templates
        print "*** Retrieving Templates ***"
        templates = {}
        for key in f_templ.GetListOfKeys():
            hist = key.ReadObj()
            print hist.GetName() 
            v = hist.GetName().split("_")[1]
            templates[v] = hist
        
        
        xsec = {} #in fb
        original_nevents = {}
        #nevents={"run_01":30000,"run_02":0}
        
        #First extract xsections
        print "*** Retrieving cross sections ***"
        for h in histograms:
            name = h.GetName()
            coupl_values = name.split("_")[2] + "_" + name.split("_")[4]
            xsec[coupl_values] = extract_xsec(args.XsecDir,c,coupl_values)*1000
            original_nevents[coupl_values] = get_Original_nevents(args.XsecDir,c,coupl_values)

        
        # now scan the coupling values
       #  print "*** Starting template fits ***"
#         for h in histograms:
#             name = h.GetName()
#             couplings = name.split("_")[1] + "_" + name.split("_")[3]
#             coupl_values = name.split("_")[2] + "_" + name.split("_")[4]
#             #if not coupl_values == "1p0_1p0": continue
#             print "Processing %s, %s"%(couplings,coupl_values)
#             fitted_SM, sigma_fitted_SM, fitted_tL,sigma_fitted_tL , fitted_tR, sigma_fitted_tR = fitTemplate(h.Clone(),templates["SM"].Clone(),templates["tL"].Clone(),templates["tR"].Clone(),couplings,coupl_values,args.OutputDir, xsec[coupl_values], original_nevents[coupl_values])
#             fitting_dict[c][coupl_values]=[fitted_SM, sigma_fitted_SM, fitted_tL,sigma_fitted_tL , fitted_tR, sigma_fitted_tR,xsec[coupl_values],original_nevents[coupl_values]]
#             print fitting_dict[c][coupl_values]
#         pickle.dump(fitting_dict,open(args.ValidationDir + "/templates/result_dict.pkl",'wb'))
        fitting_dict = pickle.load(open(args.ValidationDir + "/templates/result_dict.pkl",'r'))

        
        #sys.exit(1)   
        
        #int_lumi = 300 #fb-1
        #SM_expected = fitting_dict[c]["0p0_0p0"][2]*int_lumi*fitting_dict[c]["0p0_0p0"][6]/fitting_dict[c]["0p0_0p0"][7]
        # SM_expected = fitting_dict[c]["0p0_0p0"][0]
#         observed = {}
#         observed_unc = {}
#         #chi2 = {}
#         couplings_numerical = {}
#         for v,results in fitting_dict[c].iteritems():
#             #if v == "0p0": continue
#             observed[v] = (fitting_dict[c][v][0]-SM_expected)**2/fitting_dict[c][v][1]**2 + (fitting_dict[c][v][2])**2/fitting_dict[c][v][3]**2 + (fitting_dict[c][v][4])**2/fitting_dict[c][v][5]**2
#             #observed[v] = (fitting_dict[c][v][2]+fitting_dict[c][v][4])/math.sqrt((fitting_dict[c][v][1])**2+(fitting_dict[c][v][3])**2+(fitting_dict[c][v][5])**2)
#             #observed[v] = (fitting_dict[c][v][2]+fitting_dict[c][v][4])/math.sqrt((fitting_dict[c][v][0])+(fitting_dict[c][v][2])+(fitting_dict[c][v][4]))
#             #observed[v] = (fitting_dict[c][v][2]+fitting_dict[c][v][4])
#             #observed[v] = (fitting_dict[c][v][2]+fitting_dict[c][v][4])*int_lumi*fitting_dict[c][v][6]/fitting_dict[c][v][7]
#             observed_unc[v] = (math.sqrt((fitting_dict[c][v][3])**2+(fitting_dict[c][v][5])**2))
#             #observed_unc[v] = (math.sqrt((fitting_dict[c][v][3])**2+(fitting_dict[c][v][5])**2)*int_lumi*fitting_dict[c][v][6]/fitting_dict[c][v][7])
#             #chi2[v] = (observed[v])**2/observed_unc[v]**2
#             #chi2[v] = (observed[v]-SM_expected)**2/observed_unc[v]**2
#             couplings_numerical[v]=[convert_coupling_toNumerical(v.split("_")[0]),convert_coupling_toNumerical(v.split("_")[1])]
#         #print observed
#         #print observed_unc
#         #print couplings_numerical
        
        
        coupling_strengths_x1 = []
        coupling_strengths_x2 = []
        observed_SM_arr = []
        observed_SM_err_arr = []
        observed_tL_arr = []
        observed_tL_err_arr = []
        observed_tR_arr = []
        observed_tR_err_arr = []
        for v,results in fitting_dict[c].iteritems():
            coupling_strengths_x1.append(convert_coupling_toNumerical(v.split("_")[0]))
            coupling_strengths_x2.append(convert_coupling_toNumerical(v.split("_")[1]))
            observed_SM_arr.append(fitting_dict[c][v][0])
            observed_SM_err_arr.append(fitting_dict[c][v][1])
            observed_tL_arr.append(fitting_dict[c][v][2])
            observed_tL_err_arr.append(fitting_dict[c][v][3])
            observed_tR_arr.append(fitting_dict[c][v][4])
            observed_tR_err_arr.append(fitting_dict[c][v][5])
#             observed_err_arr.append(2*observed_unc[v])
            
        # for v, obs in observed.iteritems():
#             coupling_strengths_x1.append(couplings_numerical[v][0])
#             coupling_strengths_x2.append(couplings_numerical[v][1])
#             observed_arr.append(obs)
#             observed_err_arr.append(2*observed_unc[v])
        
        #observed_low_array = [i-j for i,j in zip(observed_arr,observed_err_arr)]
        
        coeff_SM, covmat_SM = curve_fit(func,np.vstack((coupling_strengths_x1, coupling_strengths_x2)),observed_SM_arr,sigma=observed_SM_err_arr)#absolute_sigma
        #errors_SM = sqrt(diag(covmat_SM))
        errors_SM, covmat_errors_SM = curve_fit(func,np.vstack((coupling_strengths_x1, coupling_strengths_x2)),observed_SM_err_arr)
        
        coeff_tL, covmat_tL = curve_fit(func,np.vstack((coupling_strengths_x1, coupling_strengths_x2)),observed_tL_arr,sigma=observed_tL_err_arr)#absolute_sigma
        #errors_tL = sqrt(diag(covmat_tL))
        errors_tL, covmat_errors_tL = curve_fit(func,np.vstack((coupling_strengths_x1, coupling_strengths_x2)),observed_tL_err_arr)
        
        coeff_tR, covmat_tR = curve_fit(func,np.vstack((coupling_strengths_x1, coupling_strengths_x2)),observed_tR_arr,sigma=observed_tR_err_arr)#absolute_sigma
        #errors_tR = sqrt(diag(covmat_tR))
        errors_tR, covmat_errors_tR = curve_fit(func,np.vstack((coupling_strengths_x1, coupling_strengths_x2)),observed_tR_err_arr)
        
        
        
        x1min = -4
        x1max = +4
        x2min = -4
        x2max = +4
        nsteps = 200
        x1 = np.arange(x1min, x1max+ float(x1max-x1min)/float(nsteps), float(x1max-x1min)/float(nsteps))
        x2 = np.arange(x2min, x2max+ float(x2max-x2min)/float(nsteps), float(x2max-x2min)/float(nsteps))
        side_x1 = np.linspace(x1min, x1max+ float(x1max-x1min)/float(nsteps), nsteps)
        side_x2 = np.linspace(x2min, x2max+ float(x2max-x2min)/float(nsteps), nsteps)
        X1, X2 = np.meshgrid(side_x1, side_x2)
        size = X1.shape
        x1_1d = X1.reshape((1, np.prod(size)))
        x2_1d = X2.reshape((1, np.prod(size)))
        xdata = np.vstack((x1_1d, x2_1d))
        
        
        
        y_SM = func(xdata,*(coeff_SM))
        y_SM_down = func(xdata,*(coeff_SM - [i for i in errors_SM]))
        #err_SM = np.asarray([abs(i-j) for i,j in zip(y_SM,y_SM_down)])
        err_SM = func(xdata,*(errors_SM))
        Y_SM = y_SM.reshape(size)
        y_SM_expected = func([0,0],*(coeff_SM))
        
        print y_SM, err_SM, y_SM_expected
        
        y_tL = func(xdata,*(coeff_tL))
        y_tL_down = func(xdata,*(coeff_tL - [i for i in errors_tL]))
        #err_tL = np.asarray([abs(i-j) for i,j in zip(y_tL,y_tL_down)])
        err_tL = func(xdata,*(errors_tL))
        Y_tL = y_tL.reshape(size)
        
        print y_tL, err_tL
        
        y_tR = func(xdata,*(coeff_tR))
        y_tR_down = func(xdata,*(coeff_tR - [i for i in errors_tR]))
        #err_tR = np.asarray([abs(i-j) for i,j in zip(y_tR,y_tR_down)])
        err_tR = func(xdata,*(errors_tR))
        Y_tR = y_tR.reshape(size)
        
        print y_tR, err_tR
        
        #Observation = SM!!
        obs1 = 0.0
        obs2 = 0.0
        coupl_values = "0p0_0p0"
        couplings = "cQb1_ctb1"
        y_SM_obs = fitting_dict[couplings][coupl_values][0]
        err_SM_obs = fitting_dict[couplings][coupl_values][1]
        y_tL_obs = fitting_dict[couplings][coupl_values][2]
        err_tL_obs = fitting_dict[couplings][coupl_values][3]
        y_tR_obs = fitting_dict[couplings][coupl_values][4]
        err_tR_obs = fitting_dict[couplings][coupl_values][5]
        
        
        chi2 = np.asarray([((y_tL_obs-i2)**2/err_tL_obs**2) + ((y_tR_obs-i3)**2/err_tR_obs**2) for i2,i3 in zip(y_tL,y_tR) ])
        #chi2 = np.asarray([((i2)**2/e2**2) + ((i3)**2/e3**2) for i1,e1,i2,e2,i3,e3 in zip(y_SM,err_SM,y_tL,err_tL,y_tR,err_tR) ])
        #chi2 = np.asarray([ ((i2+i3)**2/math.sqrt(e2**2+e3**2)**2) for i1,e1,i2,e2,i3,e3 in zip(y_SM,err_SM,y_tL,err_tL,y_tR,err_tR) ])
        
        #print chi2
        Y_chi2 = chi2.reshape(size)
        
        #print Y_chi2
        
        levels = [5.991]
        CS = plt.contour(X1, X2, Y_chi2,levels,colors=('r'))
        #levels = [1,2,3,4,5,6,7,8,10]
        #CS = plt.contourf(X1, X2, Y_chi2,levels)
        p = CS.collections[0].get_paths()[0]
        v = p.vertices
        px = v[:,0]
        py = v[:,1]
        
        #plt.show()
        
        
        save_dict_gr = {"px":px,"py":py}
        pickle.dump(save_dict_gr,open("%s/contour_Marginalized_templateFit_%s.pkl"%(args.OutputDir,c),"wb"))
        print "saved contour in '%s/contour_Marginalized_templateFit_%s.pkl'"%(args.OutputDir,c)
    
        
        
        gr = ROOT.TGraph(len(px),array('d',px),array('d',py))
        gr.SetMarkerSize(0)
        gr.SetFillColor(2)
        gr.SetLineColor(2)
        gr.SetLineWidth(3)
    
        mg = ROOT.TMultiGraph("mg",";%s [TeV^{-2}];%s [TeV^{-2}]"%(convertToLatex(c.split("_")[0]),convertToLatex(c.split("_")[1])))
        mg.Add(gr,"lc")
    
        ROOT.TGaxis.SetMaxDigits(2)
        canvas = ROOT.TCanvas("canvas","canvas",800,900)
        ROOT.gPad.SetMargin(0.15,0.05,0.15,0.2)
        mg.Draw("AC")

        mg.SetMinimum(x2min)
        mg.SetMaximum(x2max)
        mg.GetYaxis().CenterTitle()
        mg.GetYaxis().SetTitleSize(0.06)
        mg.GetYaxis().SetTitleOffset(1.05)
        mg.GetYaxis().SetLabelSize(0.05)
        mg.GetXaxis().CenterTitle()
        mg.GetXaxis().SetTitleSize(0.06)
        mg.GetXaxis().SetTitleOffset(1.1)
        mg.GetXaxis().SetLabelSize(0.05)
        mg.GetXaxis().SetLimits(x1min,x1max)
    
        l = ROOT.TLegend(0.15,0.8,0.95,0.93)
        l.SetBorderSize(1)
        l.SetNColumns(3)
        l.SetTextSize(0.05)
        l.AddEntry(gr,"template fit","l")
    
        l.Draw("same")
    
    
        
        #c.SaveAs('%s/fitted_xsec_Marginalized/fit_xsec_Marginalized_%s.png'%(args.OutputDir,coupling))
        canvas.SaveAs('%s/fit_xsec_Marginalized_templateFit_%s.pdf'%(args.OutputDir,c))
    
        print "done, saved in '%s/fit_xsec_Marginalized_templateFit_%s.pdf'"%(args.OutputDir,c)
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        # obs_low95CL = []
#         obs_high95CL = []
#         observed_arr = []
#         unc_arr = []
#         c_array= []
#         for v,res in observed.iteritems():
#             observed_arr.append(res)
#             obs_low95CL.append(res-2*observed_unc[v]) 
#             obs_high95CL.append(res+2*observed_unc[v])
#             c_array.append(couplings_numerical[v]) 
#             unc_arr.append(observed_unc[v])
#         
#         plt.errorbar(c_array, observed_arr, yerr=[[abs(observed_arr[i] - obs_low95CL[i]) for i in range(len(observed_arr))], [abs(observed_arr[i] - obs_high95CL[i]) for i in range(len(observed_arr))]], fmt='o', label='data 95% CL')
#         
#         coeff_centr, covmat_centr = curve_fit(func4,c_array,observed_arr,sigma=unc_arr)#,sigma=[0]*len(c_array))#absolute_sigma
#         errors_centr = sqrt(diag(covmat_centr))
#         coeff_up, covmat_up = curve_fit(func4,c_array,obs_high95CL)#,sigma=[0]*len(c_array))#absolute_sigma
#         errors_up = sqrt(diag(covmat_up))
#         coeff_down, covmat_down = curve_fit(func4,c_array,obs_low95CL)#,sigma=[0]*len(c_array))#absolute_sigma
#         errors_down = sqrt(diag(covmat_down))
#         
#         xmin = -15
#         xmax = +15
#         ymax = max(obs_high95CL)
#         ymin = -ymax/4.
#         nsteps = 30
#         x = np.arange(xmin, xmax+ float(xmax-xmin)/float(nsteps), float(xmax-xmin)/float(nsteps))
#         #y = coeff[0] + coeff[1]*x + coeff[2]*x*x
#         #y_centr = coeff_centr[0]+ coeff_centr[1]*x + coeff_centr[2]*x*x
#         y_centr = func4(x,*(coeff_centr))
#         y_up = func4(x,*(coeff_up))
#         y_down = func4(x,*(coeff_down))
#         #y_up = (coeff_centr[0]+errors_centr[0])*(1 + (coeff_centr[1]+errors_centr[1])*x + (coeff_centr[2]+errors_centr[2])*x*x)
#         #y_down = (coeff_centr[0]-errors_centr[0])*(1 + (coeff_centr[1]-errors_centr[1])*x + (coeff_centr[2]-errors_centr[2])*x*x)
#         #y_up = func2(x,*(coeff_centr + [2*i for i in errors_centr]))
#         #y_down = func2(x,*(coeff_centr - [2*i for i in errors_centr]))
#         
#         plt.plot(x,y_centr,c="r",label='fit',linestyle="dotted")
#         plt.fill_between(x, y_down,y_up, facecolor='green', alpha=0.3, label=r'fitted 95% CL')
#         plt.axhline(y=0, linestyle="dashed", linewidth=2, color="navy")
#         p = P.fit(x, y_down, 4)
#         roots = sorted([i.real for i in (p - 0).roots() if i.imag == 0])
#         print roots
#         while len(roots)>2: roots = roots[1:-1]
#         plt.axvline(x=roots[0], linestyle="dashed", linewidth=2, color="navy")
#         plt.axvline(x=roots[1], linestyle="dashed", linewidth=2, color="navy")
#         plt.legend(loc="upper right")
#     
#     
#         # plot points with error bars
#         #plt.errorbar(coupling_strengths, xsec, yerr = xsec_error, fmt= 'o', color = "b")
#         #if roots[0] < xmin: plt.axis([roots[0]-1, roots[1]+1, ymin, ymax])
#         #else: 
#         plt.axis([xmin, xmax, ymin, ymax])
#         plt.grid(True)
#         plt.xlabel(c, fontsize = 15)
#         plt.ylabel('fitted number of EFT events', fontsize = 15)
#     
#         # draw some things on the canvas
#         #plt.text(xmin, ymax+(ymax-ymin)/float(20), r'$\sigma=%.4f ( 1 + (%.4f) \ C + (%.4f) \ C^2 )$'%(coeff[0],coeff[1]/coeff[0],coeff[2]/coeff[0]), fontsize=20, color="r")
#         #plt.text(xmin, ymax+(ymax-ymin)/float(20), r'$\sigma=%.4f ( 1 + (%.4f) \ C + (%.5f) \ C^2 )$'%(coeff[0],coeff[1],coeff[2]), fontsize=12, color="r")
#         #plt.text(xmin + (xmax-xmin)/float(20), ymax-1.5*(ymax-ymin)/float(20), r'%i events'%(sum(nevents.values())), fontsize=17, color="b")
#         plt.text(roots[0] , ymin+(ymax-ymin)/float(20), "%.2f "%roots[0], fontsize=12, color="navy", ha = 'right')
#         plt.text(roots[1] , ymin+(ymax-ymin)/float(20), "%.2f "%roots[1], fontsize=12, color="navy", ha = 'right')
#         
#         plt.savefig('%s/fit_ConfidenceLevelInterval_%s.png'%(args.OutputDir,c))
#         plt.savefig('%s/fit_ConfidenceLevelInterval_%s.pdf'%(args.OutputDir,c))
#         
#         plt.cla()
#         
#         limits_file.write("%s & [%.2f,%.2f] \n"%(c,roots[0],roots[1]))  
#        #  print ""
# #         print SM_expected
# #         print observed
# #         print observed_unc
# #         print chi2
# #         print ""
#         
#     limits_file.close()    
#         
#     #print fitting_dict
# 
#     
#     #fitTemplate(inputdir+"/CONVERTED_DELPHES_TEST_TOP_ETA/C1Qu1p0_Order0_run_01_tag_1_delphes_events.root", inputdir+"/CONVERTED_DELPHES_TEST_TOP_ETA/C1Qu1p0_Order2_run_01_tag_1_delphes_events.root",inputdir+"/CONVERTED_DELPHES_TEST_TOP_ETA/C1Qu1p0_Order0_run_01_tag_1_delphes_events.root",template_hist_name="hist_top_eta")
#     
#     # f_data = ROOT.TFile("/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ValidationPerCoupling_WithInterference/templates/validation_data.root")
# #     data_hist = ROOT.TH1D()
# #     data_hist = f_data.Get("h_C1tq_10p0_outLLRR")
# #     f_EFT_templ = ROOT.TFile("/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_Training_WithInterference/templates/templates.root")
# #     EFT_templ = ROOT.TH1D()
# #     EFT_templ = f_EFT_templ.Get("h_C1tq_outLLRR")
# #     f_SM_templ = ROOT.TFile("/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_Training_WithInterference/templates/templates.root")
# #     SM_templ = ROOT.TH1D()
# #     SM_templ = f_SM_templ.Get("h_SM_outLLRR")
# #     #SM_templ.Scale(10.)
# #     
# #     outdir="/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ValidationPerCoupling_WithInterference/templates/"
# #     fitTemplate(data_hist,EFT_templ,SM_templ,"C1tq","10p0",outdir)


if __name__ == "__main__":
    main()