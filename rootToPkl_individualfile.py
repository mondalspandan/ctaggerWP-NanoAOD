import os, sys
import uproot3 as uproot
#import pyxrootd
import pandas as pd, awkward as ak
import numpy as np, pickle
from glob import glob

# di= "/pnfs/desy.de/cms/tier2//store/mc/RunIIFall17NanoAODv5/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/"
# flist17 = glob(di+"*/*.root")

infile = sys.argv[1]
maxjets = 3000000
taggers = ["DeepFlav","PNet","RobustParTAK4"]
outdir = "/eos/user/s/spmondal/ctagWP/2022"
os.system("mkdir -p "+outdir)

branches = [
"Jet_pt",
"Jet_hadronFlavour",
"Jet_partonFlavour",
"Jet_jetId",
# "Jet_puId",   # Run 2/CHS only
"Jet_eta",
"Jet_phi"]

for tag in taggers:
    for suff in ["B","CvL","CvB","C"]:                 # Remove "C" for stock Nano
        if tag=="PNet" and suff=="C": suff = "ProbC"
        branches.append("Jet_btag%s%s"%(tag,suff))
print (branches)

outf = "%s/%s.pkl"%(outdir,infile.split("/")[-1].rstrip('.root').replace('*','_'))

sumdf = pd.DataFrame()
for df in uproot.iterate(infile,'Events', 
    branches,
    entrysteps=300000,
    outputtype=pd.DataFrame,
    flatten=True):   

    sumdf = pd.concat([sumdf, df])    
    print("Added %d, total %d"%(df.shape[0],sumdf.shape[0]))
    if sumdf.shape[0] > maxjets: break

print ("Total:", sumdf.shape[0])
pickle.dump(sumdf,open(outf,"wb"))
