import os, sys
import uproot 
#import pyxrootd
import pandas as pd 
import numpy as np
from glob import glob

# di= "/pnfs/desy.de/cms/tier2//store/mc/RunIIFall17NanoAODv5/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/"
# flist17 = glob(di+"*/*.root")

infile = sys.argv[1]

branches = [
"Jet_btagDeepFlavC",
"Jet_btagDeepFlavB",
#"Jet_btagDeepFlavUDS",
"Jet_btagDeepC",
"Jet_btagDeepB",
"Jet_pt",
"Jet_btagCMVA",
"Jet_btagCSVV2",
#"GenJet_hadronFlavour",
"Jet_qgl",
"Jet_hadronFlavour",
"Jet_partonFlavour",
"Jet_jetId",
"Jet_puId",
"Jet_eta",
"Jet_phi"]

outdir = "/nfs/dust/cms/user/spmondal/ROOTtoPKL/ctagWP"
os.system("mkdir -p "+outdir)

outf = "%s/%s.pkl"%(outdir,infile.split("/")[-1].rstrip('.root'))

sumdf = pd.DataFrame()
for df in uproot.iterate(infile, 'Events', 
    branches,
    entrysteps=300000,
    flatten=True,
    outputtype=pd.DataFrame
                        ):
    
    print(df.shape)
    sumdf = pd.concat([sumdf, df])

print "Total:", sumdf.shape
sumdf.to_pickle(outf) 
