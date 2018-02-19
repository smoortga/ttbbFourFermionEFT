import os
from argparse import ArgumentParser
import sys
import numpy as np
import root_numpy as rootnp
import ROOT
from keras.models import load_model
import pickle
from rootpy.plotting import Hist
#from rootpy.root2array import fill_hist_with_ndarray


def pad(A, length):
    if len(A) > length: return np.asarray(A)
    arr = np.zeros(length)
    arr[:len(A)] = A
    return np.asarray(arr)

def retrieve_discr(input_dir,filename,model,scaler):
    if len(filename) == 0:
        print "ERROR: no files found"
        sys.exit()
    X = rootnp.root2array(input_dir + "/" + filename[0],"tree")
    X = rootnp.rec2array(X)
    for i in range(len(filename)):
        if i == 0: continue
        X_ = rootnp.root2array(input_dir + "/" + filename[i],"tree")
        X_ = rootnp.rec2array(X_)
        X = np.concatenate((X,X_))
    
    X = scaler.transform(X)
    
    discrSM = model.predict(X)[:,0]
    discrtL = model.predict(X)[:,1]
    discrtR = model.predict(X)[:,2]
    
    #print len(discrSM)
    
    discr_dict = {"SMvsEFT":np.asarray([]),"tLvstR":np.asarray([])}
    discr_dict["SMvsEFT"] = np.asarray([(discrtL[jdx]+discrtR[jdx]) for jdx in range(len(discrSM))])
    discr_dict["tLvstR"] = np.asarray([(discrtL[jdx])/(discrtL[jdx]+discrtR[jdx]) for jdx in range(len(discrSM))])
    
    return discr_dict

def main():

    parser = ArgumentParser()
    parser.add_argument('--TrainingFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Training_fixed/training_output/model_checkpoint_save.hdf5", help='path to training')
    parser.add_argument('--ScalerFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Training_fixed/training_output/scaler.pkl", help='Scaler')
    parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Training_fixed", help='path to the converted delphes training files')
    parser.add_argument('--InputDirForValidation', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Marginalized", help='path to the converted delphes training files')
    args = parser.parse_args()
    
    if not os.path.isdir(args.InputDir+"/templates/"): os.mkdir(args.InputDir+"/templates/")
    f_templates = ROOT.TFile(args.InputDir+"/templates/templates.root","RECREATE")
    
    #get training and scaler
    model_ = load_model(args.TrainingFile)
    scaler_ = pickle.load(open(args.ScalerFile,'r'))
    
    # make files_dict
    files_dict = {}
    #SM part
    files_dict["SM"] = []
    files_dict["tL"] = []
    files_dict["tR"] = []
    fileslist = [f for f in os.listdir(args.InputDir) if ".root" in f]
    for f in fileslist:
        if "SM" in f: files_dict["SM"].append(f)
        elif "t" in f.split("_")[0]: files_dict["tR"].append(f)
        else: files_dict["tL"].append(f)
    

    
    # loop over couplings
    for c,arr in files_dict.iteritems():
        print "Processing %s"%c
        discriminators = retrieve_discr(args.InputDir,arr,model_,scaler_)
        hist = ROOT.TH2D("h_%s_out"%(c),";discriminator P(t_{L}) + P(t_{R});discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",10,0,1,10,0,1)
        for smvseft,tlvstr in zip(discriminators["SMvsEFT"],discriminators["tLvstR"]):
            hist.Fill(smvseft,tlvstr)
        f_templates.cd()
        hist.Write()

    f_templates.Close()
    
    print "Output saved in",args.InputDir+"/templates/templates.root"
    
    
    
    # NOW DO THE EXACT SAME THING BUT FOR THE VALIDATION FILES (Marginalized)
    
    if not os.path.isdir(args.InputDirForValidation+"/templates/"): os.mkdir(args.InputDirForValidation+"/templates/")
    f_templates_val = ROOT.TFile(args.InputDirForValidation+"/templates/validation_data.root","RECREATE")
    
    # make files_dict
    files_dict = {}

    files = [f for f in os.listdir(args.InputDirForValidation) if ".root" in f]
    couplings = set([f.split("_")[0] + "_" + f.split("_")[2] for f in os.listdir(args.InputDirForValidation) if ".root" in f])
    coupl_strengths = set([f.split("_")[1] + "_" + f.split("_")[3] for f in os.listdir(args.InputDirForValidation) if ".root" in f])
    for c in couplings:
        files_dict[c]={}
        for s in coupl_strengths:
            files_dict[c][s]=[]

    for f in files:
        coupling_name = f.split("_")[0] + "_" + f.split("_")[2]
        strength_name = f.split("_")[1] + "_" + f.split("_")[3]
        files_dict[coupling_name][strength_name].append(f)
    
        
    # loop over couplings
    for c,dict in files_dict.iteritems():
        print "Processing %s"%c
        for val,arr in dict.iteritems():
            print "    %s"%val
            discriminators = retrieve_discr(args.InputDirForValidation,arr,model_,scaler_)
#            for output,disc in discriminators.iteritems():
            hist = ROOT.TH2D("h_%s_%s_%s_%s_out"%(c.split("_")[0],val.split("_")[0],c.split("_")[1],val.split("_")[1]),";discriminator P(t_{L}) + P(t_{R});discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",10,0,1,10,0,1)
            for smvseft,tlvstr in zip(discriminators["SMvsEFT"],discriminators["tLvstR"]):
                hist.Fill(smvseft,tlvstr)
            f_templates_val.cd()
            hist.Write()

    f_templates_val.Close()
    
    print "Output saved in",args.InputDirForValidation+"/templates/validation_data.root"
    
    
    

if __name__ == "__main__":
    main()