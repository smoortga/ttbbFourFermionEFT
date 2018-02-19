import os
from argparse import ArgumentParser
import sys
#import matplotlib.pyplot as plt
import numpy as np
from numpy import polyfit, diag, sqrt
from scipy.optimize import curve_fit
from numpy.polynomial import Polynomial as P
import numpy as np
# import root_numpy as rootnp
# from keras.models import load_model
import pickle
import sys
import ROOT
import math
from array import array
# plt.rc('text', usetex=True)
# plt.rc('font', family='sans-serif')

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

parser = ArgumentParser()

parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttbb_DiLepton_ForValidationPerCouplings_fixed", help='path to the directory of all restricted processes')
parser.add_argument('--ValidationDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_ForValidationPerCouplings_fixed", help='path to the directory of all restricted processes')

args = parser.parse_args()

classes_dict = { #name:class_number
            'SM': 0,
            #LLLL
            "cQQ1": 1,
            "cQQ8": 1,
            #LLRR
            "cQt1": 3,
            "cQb1": 3,
            "cQt8": 3,
            "cQb8": 3,
            #RRRR
            "ctb1": 2,
            "ctb8": 2
        }

def get_Original_nevents(filename,nevents_dict):
    for run,n in nevents_dict.iteritems():
        if run in filename: return n
     
    print "ERROR: couldn't retrieve original number of events from file %s"%filename
    print "Returning -1"
    return -1

def convertToLatex(c):
    if len(c) == 4: return c[0].upper()+"_{"+c[1:3]+"}^{"+c[3:4]+"}"
    else: 
        print "couldn't convert"
        return c  

def GetEventsPassingCut(filelist,variable_name,xmin,xmax,nbins,cut,original_n_events):
    chain = ROOT.TChain("tree")
    for f in filelist:
        chain.Add(args.ValidationDir+"/"+f)
    #print chain.GetEntries()
    tmp_hist = ROOT.TH1D("tmp_hist","",nbins,xmin,xmax)
    chain.Draw(variable_name + " >> tmp_hist")
    #print tmp_hist.GetEntries()
    selected_n_events = tmp_hist.Integral(tmp_hist.FindBin(cut),tmp_hist.GetNbinsX()+1)
    result = float(selected_n_events)/float(original_n_events)
    print "fracion selected: %.2f%%"%(result*100)
    return result

def extract_coupling(coupling,val):
    """ extract the numeric value of a coupling in the directory name"""
    if "minus" in val: return -float(val.replace(coupling+"_minus","").replace("p","."))
    else: return float(val.replace(coupling+"_","").replace("p","."))
    
def extract_coupling_string(coupling,val):
    """ extract the string value of a coupling in the directory name"""
    return val.replace(coupling+"_","")
    
def func2(x, a, b, c):
        return (a + b*x + c*x*x)   

def func4(x, a, b, c, d, e):
        return (a + b*x + c*x*x + d*x*x*x + e*x*x*x*x)    

def func(x, a, b, c):
        return a*(1 + b*x + c*x*x) 

def makeFitFunction(p0,p1,p2,units="pb"):
    text = ""
    if units == "pb":
        text = text + "#sigma [pb] = "
        text = text + "%.4f ("%p0
        text = text + " 1 "
        if p1 > 0: text = text + "+ %.4f C"%math.fabs(p1)
        else: text = text + "- %.4f C"%math.fabs(p1)
        if p2 > 0: text = text + " + %.4f C^{2} )"%math.fabs(p2)
        else: text = text + "- %.4f C^{2} )"%math.fabs(p2)
    
    return text
        

    

process_dict = {}
subdirs = [d for d in os.listdir(args.InputDir) if os.path.isdir(args.InputDir+"/"+d) and not "fitted" in d and not 'Limits' in d and not "Discriminator" in d]
for subdir in subdirs:
    process_dict[subdir]= [dd for dd in os.listdir(args.InputDir+"/"+subdir) if os.path.isdir(args.InputDir+"/"+subdir+"/"+dd) and not "fitted_xsec" in dd]

#print process_dict

validation_files = {}
files = [f for f in os.listdir(args.ValidationDir) if ".root" in f]
couplings = set([f.split("_")[0] for f in os.listdir(args.ValidationDir) if ".root" in f])
coupl_strengths = set([f.split("_")[1] for f in os.listdir(args.ValidationDir) if ".root" in f])
for c in couplings:
    validation_files[c]={}
    for s in coupl_strengths:
        validation_files[c][s]=[]
for f in files:
    coupling_name = f.split("_")[0]
    strength_name = f.split("_")[1]
    validation_files[coupling_name][strength_name].append(f)


# coefficients with errors: sigma = p0 ( 1 + p1*C + p2*c^2 )
p0 = []
p0e = []
p1 = []
p1e = []
p2 = []
p2e = []
coupling_names = []


nevents={"run_01":30000,"run_02":0}

result_dict = {}


for coupling,orders in process_dict.iteritems():
    result_dict[coupling] = {}
    for o in sorted(orders):
        print "Processing %s %s"%(coupling,extract_coupling(coupling,o))
        ##################
        #
        # Get cross section
        #
        ##################
        if os.path.isfile(args.InputDir + "/" + coupling + "/" + o + "/cross_section.txt"):
            f_ = open(args.InputDir + "/" + coupling + "/" + o + "/cross_section.txt", 'r')
            lines = f_.readlines()
            xsec = float(lines[-1].split(" ")[2])
            xsec_error = float(lines[-1].split(" ")[3])
        else: continue
        #xsec = 0.01
        #xsec_error = 0.001
        
        ##################
        #
        # NN output
        #
        ##################
        
        filelist = validation_files[coupling][extract_coupling_string(coupling,o)]
        #filelist = validation_files[coupling][o.split("_")[1]]
        if len(filelist) == 0: continue
        sum_original_n_events = 0
        for f in filelist:
            sum_original_n_events += get_Original_nevents(f,nevents)
        if "8" in coupling:
            frac_passing_evts = GetEventsPassingCut(filelist,variable_name = "m_c1c2b1b2",xmin=0,xmax=3000, nbins=1000,cut=1100,original_n_events = sum_original_n_events)
        elif "1" in coupling:
            frac_passing_evts = GetEventsPassingCut(filelist,variable_name = "m_c1c2b1b2",xmin=0,xmax=3000, nbins=1000,cut=1100,original_n_events = sum_original_n_events)
        result_dict[coupling][extract_coupling_string(coupling,o)] = [xsec,xsec_error,frac_passing_evts, sum_original_n_events]
       
if not os.path.isdir(args.ValidationDir + "/validation_output"): os.mkdir(args.ValidationDir + "/validation_output")
pickle.dump(result_dict,open(args.ValidationDir + "/validation_output/SensitiveVariable_result_dict.pkl",'wb'))
#result_dict = pickle.load(open(args.ValidationDir + "/validation_output/SensitiveVariable_result_dict.pkl",'r'))

#sys.exit(1)


# find out the mean of the SM contributions as an average value
# sm_values = []
# for coupling,orders in result_dict.iteritems():
#     for o,results in orders.iteritems():
#         if o == "0p0": sm_values.append(results[0]*results[2])
# sm_xsec = np.mean(np.asarray(sm_values))
# sm_xsec_error = np.std(np.asarray(sm_values))
    
# sm_xsec = 0.085*0.11*0.05 # pb from TOP-16-010
# sm_xsec_error = 0.041*0.11*0.05# pb from TOP-16-010
sm_xsec_frac_error = 0.1
sm_xsec_frac_error_300fb = 0.1

fracerr_measured_sm = sm_xsec_frac_error

# define a file for the limits to be stored: 
limits_file = open('%s/limits_ttbb_cutonSensitiveVariable.txt'%(args.InputDir),"w")
#limits_file = open('%s/limits_ttbb_NN_Binary_cut0p75.txt'%(args.InputDir),"w")
#limits_file = open('%s/limits_ttbb_nocut.txt'%(args.InputDir),"w")
limits_file.write("fractional error on SM measurement: %i %% \n"%int(100*fracerr_measured_sm))


for coupling,values in result_dict.iteritems():
    print coupling
    #if coupling == "cQQ1": continue
    coupling_strengths = []
    xsec = []
    xsec_error = []
    #nevents = []
    for v,result in values.iteritems():
        if v == "0p0": 
            sm_xsec = result[0]*result[2]
            sm_xsec_error = math.sqrt( (result[1]*result[2])**2 + (result[0]*math.sqrt(result[3]*result[2])/result[3])**2 )
        coupling_strengths.append(extract_coupling(coupling,coupling+"_"+v))
        xsec.append(result[0]*result[2])
        xsec_error.append( math.sqrt( (result[1]*result[2])**2 + (result[0]*math.sqrt(result[3]*result[2])/result[3])**2 ) )
        #nevents.append(float(lines[-1].split(" ")[4][:-1])) #need to remove \n character in the end
    
    #print xsec,xsec_error
    coeff, covmat = curve_fit(func,coupling_strengths,xsec,sigma=xsec_error)#absolute_sigma
    errors = sqrt(diag(covmat))
    
    coupling_names.append(coupling)
    p0.append(coeff[0])
    p0e.append(errors[0])
    p1.append(coeff[1])
    p1e.append(errors[1])
    p2.append(coeff[2])
    p2e.append(errors[2])
    
    
    #**************************
    #
    #   PLOTTING
    #
    #**************************
    
    total_frac_error = math.sqrt((sm_xsec_frac_error)**2 + (math.sqrt(sm_xsec*1000*300)/(sm_xsec*1000*300))**2)
    
    
    gr_data = ROOT.TGraphErrors(len(coupling_strengths),array('d',coupling_strengths),array('d',xsec),array('d',np.zeros(len(coupling_strengths))),array('d',xsec_error))
    gr_data.SetLineWidth(1)
    gr_data.SetMarkerStyle(20)
    gr_data.SetMarkerColor(4)
    xmin = -10
    xmax = +10
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
    # gr_CMS = ROOT.TGraphErrors(len(x),array('d',x),array('d',[sm_xsec]*len(x)),array('d',np.zeros(len(x))),array('d',[2*sm_xsec_frac_error*sm_xsec]*len(x)))
#     gr_CMS.SetFillColorAlpha(3,0.5)
#     gr_CMS.SetLineWidth(2)
    # CMS result @ 300 fb-1
    gr_CMS_300fb = ROOT.TGraphErrors(len(x),array('d',x),array('d',[min(y)]*len(x)),array('d',np.zeros(len(x))),array('d',[2*total_frac_error*min(y)]*len(x)))
    gr_CMS_300fb.SetFillColorAlpha(49,0.8)
    gr_CMS_300fb.SetLineWidth(2)
    gr_CMS_300fb.SetLineStyle(2)
    
    mg = ROOT.TMultiGraph("mg",";;cross section [pb]")
    #mg.Add(gr_CMS,"3")
    mg.Add(gr_CMS_300fb,"3")
    mg.Add(gr_fitted_xsec,"c3")
    mg.Add(gr_data,"pe")
    
    ROOT.TGaxis.SetMaxDigits(2)
    c = ROOT.TCanvas("c","c",800,700)
    uppad = ROOT.TPad("u","u",0.,0.4,1.,1.)
    downpad = ROOT.TPad("d","d",0.,0.,1.,0.4)
    uppad.Draw()
    downpad.Draw()
    uppad.cd()
    ROOT.gPad.SetMargin(0.15,0.05,0.01,0.1)
    #gr_data.Draw("AP")
    mg.Draw("A")
    if "8" in coupling:
        ymin = sm_xsec - 10*total_frac_error*sm_xsec
        ymax = sm_xsec + 20*total_frac_error*sm_xsec
    elif "1" in coupling:
        ymin = sm_xsec - 10*total_frac_error*sm_xsec
        ymax = sm_xsec + 20*total_frac_error*sm_xsec
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
    l.AddEntry(gr_data,"sample points (M_{4b} > 1.1 TeV)","pe1")
    l.AddEntry(gr_fitted_xsec,"fitted cross section","l")
    #l.AddEntry(gr_CMS,"CMS result 95% CL (2.3 fb^{-1})","fl")
    l.AddEntry(gr_CMS_300fb,"CMS prospect 95% CL (300 fb^{-1})","fl")
    
    l.Draw("same")
    
    t = ROOT.TLatex()
    t.SetTextSize(0.07)
    t.SetTextFont(42)
    t.SetTextColor(2)
    t.DrawLatexNDC(0.25,0.93,makeFitFunction(coeff[0],coeff[1],coeff[2]))
    
    
    # text = ROOT.TLatex(0.7,0.9,"M_{4b} < 600 GeV")
#     text.SetTextSize(0.06)
#     text.SetTextFont(42)
#     text.DrawLatexNDC(0.7,0.9,"M_{4b} < 600 GeV")
    
    downpad.cd()
    ROOT.gPad.SetMargin(0.15,0.05,0.25,0.01)
    y_2 = [(i-min(y_down))**2/((total_frac_error*sm_xsec)**2) for i in y_down]
    ymin = -1
    ymax = 13
    # gr_chi2 = ROOT.TGraph(len(x),array('d',x),array('d',y_1))
#     gr_chi2.SetLineWidth(2)
#     gr_chi2.SetLineColor(1)
#     gr_chi2.SetLineStyle(1)
    # prospects @ 300 fb
    #y_2 = [(i-min(y))**2/((sm_xsec_frac_error_300fb*min(y))**2) for i in y_down]
    gr_chi2_300fb = ROOT.TGraph(len(x),array('d',x),array('d',y_2))
    gr_chi2_300fb.SetLineWidth(2)
    gr_chi2_300fb.SetLineColor(1)
    gr_chi2_300fb.SetLineStyle(2)
    
    
    # get roots
    # p_1 = P.fit(x, y_1, 4)
#     roots_1 = sorted([i.real for i in (p_1 - 3.84).roots() if i.imag == 0])
#     print roots_1
#     while len(roots_1)>2: roots_1 = roots_1[1:-1]
#     x_limit_1 = np.arange(roots_1[0], roots_1[1]+abs(roots_1[1]-roots_1[0])/200., abs(roots_1[1]-roots_1[0])/200.)
#     gr_CMS_limit = ROOT.TGraphErrors(len(x_limit_1),array('d',x_limit_1),array('d',[5]*len(x_limit_1)),array('d',np.zeros(len(x_limit_1))),array('d',[5000]*len(x_limit_1)))
#     gr_CMS_limit.SetFillColorAlpha(3,0.5)
#     gr_CMS_limit.SetLineWidth(2)
    
    p_2 = P.fit(x, y_2, 4)
    roots_2 = sorted([i.real for i in (p_2 - 3.84).roots() if i.imag == 0])
    print roots_2
    while len(roots_2)>2: roots_2 = roots_2[1:-1]
    x_limit_2 = np.arange(roots_2[0], roots_2[1]+abs(roots_2[1]-roots_2[0])/200., abs(roots_2[1]-roots_2[0])/200.)
    gr_CMS_limit_300fb = ROOT.TGraphErrors(len(x_limit_2),array('d',x_limit_2),array('d',[5]*len(x_limit_2)),array('d',np.zeros(len(x_limit_2))),array('d',[5000]*len(x_limit_2)))
    gr_CMS_limit_300fb.SetFillColorAlpha(49,0.8)
    gr_CMS_limit_300fb.SetLineWidth(2)
    
    
    mg2 = ROOT.TMultiGraph("mg2",";%s [TeV^{-2}];#chi^{2}"%convertToLatex(coupling))
    #mg2.Add(gr_CMS_limit,"3")
    mg2.Add(gr_CMS_limit_300fb,"3")
    #mg2.Add(gr_chi2,"l")
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
    
    
    limits_file.write("%s & [%.2f,%.2f] \n"%(coupling,roots_2[0],roots_2[1]))
    #limits_file_300fb.write("%s & [%.2f,%.2f] \n"%(coupling,roots_2[0],roots_2[1]))
    
    if not os.path.isdir(args.InputDir + "/fitted_xsec_AfterCutOnSensitiveVariable"): os.mkdir(args.InputDir + "/fitted_xsec_AfterCutOnSensitiveVariable")

    c.SaveAs('%s/fitted_xsec_AfterCutOnSensitiveVariable/fit_xsec_AfterCutOnSensitiveVariable_%s.png'%(args.InputDir,coupling))
    c.SaveAs('%s/fitted_xsec_AfterCutOnSensitiveVariable/fit_xsec_AfterCutOnSensitiveVariable_%s.pdf'%(args.InputDir,coupling))
    
    print "done, saved in %s/fitted_xsec_AfterCutOnSensitiveVariable/fit_xsec_AfterCutOnSensitiveVariable_%s.pdf"%(args.InputDir,coupling)
    
limits_file.close()
#limits_file_300fb.close() 
    
 
    

