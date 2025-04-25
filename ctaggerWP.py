#  Coded in Python 3
#  c-tagger WP measurement using NanoAOD.
#  Author: Spandan Mondal,
#  RWTH Aachen University

import pickle, glob, sys, pandas as pd, argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import mplhep as hep, numpy as np, pickle
hep.style.use(hep.style.CMS)

parser = argparse.ArgumentParser("Measurement of c-tagger WPs")
parser.add_argument('-i','--inputdir',type=str,default="/mnt/c/Work/ctaggerWP/",help="Path to input dir containing pkls (if multiple).")
parser.add_argument('-f','--inputfile',type=str,default="",help="Input file name.")
parser.add_argument('-s','--skipseeds',action="store_true",default=False,help="Skip 1D scans that identify starting points for 2D scans.")
parser.add_argument('-p','--plotoldresults',action="store_true",default=False,help="Only plot stored results, skip everything else.")
parser.add_argument('-o','--outname',type=str,default="2024v15")
parser.add_argument('-n','--njets',type=str,default="-1",help="Maximum number of jets to process.")
args = parser.parse_args()
print (args)

if args.inputfile == "":
    indir = args.inputdir
    filelist = glob.glob(indir+'/*.pkl')
else:
    filelist = [args.inputfile]

df = pd.DataFrame()
taggers = ["UParTAK4"] #"PNet"] #"DeepFlav","PNet","RobustParTAK4"]
tagname = {
#    "DeepFlav" : "DeepJet",
#    "PNet" : "ParticleNet",
#    "RobustParTAK4" : "Robust ParT",
    "UParTAK4" : "UParT"
}
xax = "CvL"
yax = "CvB"
xname = "CvsL"
yname = "CvsB"

for fl in filelist:
    thisdf = pickle.load(open(fl,'rb'))
    df = pd.concat([df,thisdf])


def cleandf(tdf):
    cdf = tdf
    print(tdf.shape)
    # cdf = tdf[(tdf['Jet_btagDeepFlavB'] < 1) & (tdf['Jet_btagDeepFlavB'] > 0)]
    if "Jet_jetId" not in cdf.keys():
        cdf["Jet_jetId"] = ((cdf["Jet_neHEF"] < 0.99) & (cdf["Jet_neEmEF"] < 0.9) & (cdf["Jet_chMultiplicity"]+cdf["Jet_neMultiplicity"] > 1) & (cdf["Jet_chHEF"] > 0.01) & (cdf["Jet_chMultiplicity"] > 0)) * 5
    cdf = cdf[ (cdf['Jet_pt'] > 20)  & (abs(cdf['Jet_eta']) < 2.5) & (cdf['Jet_jetId'] >= 5)]
    njets = int(args.njets)
    print(cdf.shape)
    if njets > 0: cdf = cdf[:njets]    
    return cdf

def getTaggers(df):
    df['truthb'] = (df['Jet_hadronFlavour'] == 5).astype(int)
    df['truthc'] = (df['Jet_hadronFlavour'] == 4).astype(int)
    df['truthudsg'] = (df['Jet_hadronFlavour'] < 4).astype(int)

 #   df['Jet_btagPNetC'] = df['Jet_btagPNetProbC']
 #   for tag in taggers:
 #       df["Jet_btag"+tag+"CvsUDSGL"] = df["Jet_btag"+tag+"C"]/(1-df["Jet_btag"+tag+"B"])
    
    # df['dCvL'] = df['Jet_btagDeepC']/(1 - df['Jet_btagDeepB'])
    # df['dCvB'] = df['Jet_btagDeepC']/(df['Jet_btagDeepC'] + df['Jet_btagDeepB'])
    # df['dBvL'] = df['Jet_btagDeepB']/(1 - df['Jet_btagDeepC'])
    # df['dBvC'] = df['Jet_btagDeepB']/(df['Jet_btagDeepC'] + df['Jet_btagDeepB'])

    # df['dfCvL'] = df['Jet_btagDeepFlavC']/(1 - df['Jet_btagDeepFlavB'])
    # df['dfCvB'] = df['Jet_btagDeepFlavC']/(df['Jet_btagDeepFlavC'] + df['Jet_btagDeepFlavB'])
    # df['dfBvL'] = df['Jet_btagDeepFlavB']/(1 - df['Jet_btagDeepFlavC'])
    # df['dfBvC'] = df['Jet_btagDeepFlavB']/(df['Jet_btagDeepFlavC'] + df['Jet_btagDeepFlavB'])

    return df

df = cleandf(df)
df = getTaggers(df)
print (df)

def getEff(flav,pref,cvsl=0.,cvsb=0.):
    alljets = df[df['truth'+flav]==1]
    njets = alljets.shape[0]
    passjets = alljets[(alljets[pref+xax]>=cvsl) & (alljets[pref+yax]>=cvsb)]
    return float(passjets.shape[0])/njets


#Set these, everything else *should* be automatic
mistags = { 
    'Loose':  {'b':.35,'udsg':.9},
    'Medium': {'b':.25,'udsg':.25},
    'Tight':  {'b':.20,'udsg':.03},
    'XT':     {'b':.10,'udsg':.01},
}

# Manual seeds, e.g., from a previous campaign. Overwritten by next section
# seeds ={'Loose':
#             {'d': {'cvl' : 0.072, 'cvb' : 0.415},
#             'df': {'cvl' : 0.035, 'cvb' : 0.568}},
#         'Medium':
#             {'d': {'cvl' : 0.26, 'cvb' : 0.07},
#             'df': {'cvl' : 0.27, 'cvb' : 0.33}},
#         'Tight':
#             {'d': {'cvl' : 0.23, 'cvb' : 0.418},
#             'df': {'cvl' : 0.27, 'cvb' : 0.58}},
#         }
# seeds = {'Loose': {'d': {'cvl': 0.1, 'cvb': 0.22}, 'df': {'cvl': 0.05, 'cvb': 0.32}}, 'Medium': {'d': {'cvl': 0.19, 'cvb': 0.35000000000000003}, 'df': {'cvl': 0.11, 'cvb': 0.49}}, 'Tight': {'d': {'cvl': 0.43, 'cvb': 0.43}, 'df': {'cvl': 0.29, 'cvb': 0.58}}}

#2022EE
#seeds = {'Loose': {'Jet_btagDeepFlav': {'cvl': 0.041999999999999996, 'cvb': 0.207}, 'Jet_btagPNet': {'cvl': 0.05399999999999999, 'cvb': 0.183}, 'Jet_btagRobustParTAK4': {'cvl': 0.03874999999999999, 'cvb': 0.06799999999999999}}, 'Medium': {'Jet_btagDeepFlav': {'cvl': 0.10749999999999998, 'cvb': 0.29999999999999993}, 'Jet_btagPNet': {'cvl': 0.16, 'cvb': 0.3069999999999999}, 'Jet_btagRobustParTAK4': {'cvl': 0.11649999999999999, 'cvb': 0.12999999999999992}}, 'Tight': {'Jet_btagDeepFlav': {'cvl': 0.3, 'cvb': 0.24399999999999977}, 'Jet_btagPNet': {'cvl': 0.485, 'cvb': 0.2619999999999998}, 'Jet_btagRobustParTAK4': {'cvl': 0.3535, 'cvb': 0.09799999999999981}}}

#2022
#seeds = {'Loose': {'Jet_btagDeepFlav': {'cvl': 0.041999999999999996, 'cvb': 0.207}, 'Jet_btagPNet': {'cvl': 0.05099999999999999, 'cvb': 0.181}, 'Jet_btagRobustParTAK4': {'cvl': 0.03874999999999999, 'cvb': 0.06799999999999999}}, 'Medium': {'Jet_btagDeepFlav': {'cvl': 0.10749999999999998, 'cvb': 0.29999999999999993}, 'Jet_btagPNet': {'cvl': 0.1525, 'cvb': 0.3029999999999999}, 'Jet_btagRobustParTAK4': {'cvl': 0.11649999999999999, 'cvb': 0.12999999999999992}}, 'Tight': {'Jet_btagDeepFlav': {'cvl': 0.2935, 'cvb': 0.30999999999999983}, 'Jet_btagPNet': {'cvl': 0.45899999999999996, 'cvb': 0.3599999999999999}, 'Jet_btagRobustParTAK4': {'cvl': 0.352, 'cvb': 0.12999999999999984}}}

#NanoAODv11 2022
#seeds = {'Loose': {'Jet_btagDeepFlav': {'cvl': 0.03668749999999999, 'cvb': 0.27681249999999996}}, 'Medium': {'Jet_btagDeepFlav': {'cvl': 0.09838867560029031, 'cvb': 0.36874999999999997}}, 'Tight': {'Jet_btagDeepFlav': {'cvl': 0.33549999999999996, 'cvb': 0.31199999999999983}}}

#NanoAODv11 2022EE
#seeds = {'Loose': {'Jet_btagDeepFlav': {'cvl': 0.03665161319077015, 'cvb': 0.27731249999999996}}, 'Medium': {'Jet_btagDeepFlav': {'cvl': 0.09812499999999982, 'cvb': 0.37024999999999997}}, 'Tight': {'Jet_btagDeepFlav': {'cvl': 0.33549999999999996, 'cvb': 0.31249999999999983}}}

#2023
# seeds = {'Loose': {'Jet_btagDeepFlav': {'cvl': 0.04183593750000002, 'cvb': 0.23450000000000001}, 'Jet_btagPNet': {'cvl': 0.052089843749999996, 'cvb': 0.22050000000000003}, 'Jet_btagRobustParTAK4': {'cvl': 0.0380859375, 'cvb': 0.08650000000000001}}, 'Medium': {'Jet_btagDeepFlav': {'cvl': 0.10192187499999998, 'cvb': 0.32199999999999995}, 'Jet_btagPNet': {'cvl': 0.1480703125000003, 'cvb': 0.3529999999999999}, 'Jet_btagRobustParTAK4': {'cvl': 0.10924999999999999, 'cvb': 0.15274999999999994}}, 'Tight': {'Jet_btagDeepFlav': {'cvl': 0.24999999999999994, 'cvb': 0.2619999999999998}, 'Jet_btagPNet': {'cvl': 0.43399999999999994, 'cvb': 0.2999999999999989}, 'Jet_btagRobustParTAK4': {'cvl': 0.30849999999999994, 'cvb': 0.11274999999999982}}, 'XT': {'Jet_btagDeepFlav': {'cvl': 0.37099999999999994, 'cvb': 0.4405000000000001}, 'Jet_btagPNet': {'cvl': 0.634499999999997, 'cvb': 0.5490000000000002}, 'Jet_btagRobustParTAK4': {'cvl': 0.469, 'cvb': 0.2750000000000001}}}

#2023BPix
seeds = {'Loose': {'Jet_btagDeepFlav': {'cvl': 0.04186718750000002, 'cvb': 0.24200000000000002}, 'Jet_btagPNet': {'cvl': 0.052214843749999997, 'cvb': 0.22800000000000004}, 'Jet_btagRobustParTAK4': {'cvl': 0.038101562500000005, 'cvb': 0.09125000000000001}}, 'Medium': {'Jet_btagDeepFlav': {'cvl': 0.10204687499999998, 'cvb': 0.32799999999999996}, 'Jet_btagPNet': {'cvl': 0.14869531250000032, 'cvb': 0.3582499999999999}, 'Jet_btagRobustParTAK4': {'cvl': 0.1095581091940403, 'cvb': 0.15749999999999995}}, 'Tight': {'Jet_btagDeepFlav': {'cvl': 0.25049999999999994, 'cvb': 0.26674999999999977}, 'Jet_btagPNet': {'cvl': 0.43599999999999994, 'cvb': 0.3029999999999989}, 'Jet_btagRobustParTAK4': {'cvl': 0.30849999999999994, 'cvb': 0.11574999999999983}}, 'XT': {'Jet_btagDeepFlav': {'cvl': 0.37099999999999994, 'cvb': 0.4440000000000001}, 'Jet_btagPNet': {'cvl': 0.634499999999997, 'cvb': 0.5520000000000002}, 'Jet_btagRobustParTAK4': {'cvl': 0.469, 'cvb': 0.2805000000000001}}}

#1D scans
if not args.skipseeds and not args.plotoldresults:
    seeds = {}
    for wp in mistags:
        seeds[wp] = {}
        for pref in taggers:
            tag = "Jet_btag"+pref
            seeds[wp][tag] = {}
            lspace = np.linspace(0,1,101)
            
            leffs = [1.1]
            target = mistags[wp]['udsg']
            for cvsl in lspace:
                eff = getEff('udsg',tag,cvsl,0)
                leffs.append(eff)
                if eff < target: break

            idx = np.abs(np.array(leffs)-target).argmin()
            seeds[wp][tag]['cvl'] = lspace[idx]

            beffs = [1.1]
            target = mistags[wp]['b']
            for cvsb in lspace:
                eff = getEff('b',tag,0,cvsb)
                beffs.append(eff)
                if eff < target: break

            idx = np.abs(np.array(beffs)-mistags[wp]['b']).argmin()
            seeds[wp][tag]['cvb'] = lspace[idx]

print ("Using seeds:",seeds)


#2D scan
accuracy = 1e-4
step = 1e-3
def findcuts(tag,wp):
    cvsl = seeds[wp][tag]['cvl']
    cvsb = seeds[wp][tag]['cvb']
    targetl = mistags[wp]['udsg']
    targetb = mistags[wp]['b']

    stepcvl = step
    stepcvb = step
    leff = getEff('udsg',tag,cvsl,cvsb)
    beff = getEff('b',tag,cvsl,cvsb)
    if leff < targetl: leffwasless = True
    else: leffwasless = False
    if beff < targetb: beffwasless = True
    else: beffwasless = False

    count = 0
    lastl = 0
    lastb = 0
    while abs(leff-targetl)>accuracy or abs(beff-targetb)>accuracy :
        leff = getEff('udsg',tag,cvsl,cvsb)
        beff = getEff('b',tag,cvsl,cvsb)
        if abs(leff-targetl)>accuracy:
            if leff < targetl:
                if not leffwasless: stepcvl /= 2
                cvsl -= stepcvl
                leffwasless = True
            else:
                if leffwasless: stepcvl /= 2
                cvsl += stepcvl
                leffwasless = False

        if abs(beff-targetb)>accuracy:
            if beff < targetb:
                if not beffwasless: stepcvb /= 2
                cvsb -= stepcvb
                beffwasless = True
            else:
                if beffwasless: stepcvb /= 2
                cvsb += stepcvb
                beffwasless = False

        if cvsl< 0 or cvsb < 0:
            print ("Did not converge! Try lowering mistag rates!")
            sys.exit(1)
        
        count+=1
        if count%50 == 0:
            print("    - Completed %d iterations. Diffs: %f, %f."%(count,abs(leff-targetl),abs(beff-targetb)))
            if abs(leff-targetl) == lastl and abs(beff-targetb)==lastb:
                print("    - No change. Quitting.")
                break
            lastl = abs(leff-targetl)
            lastb = abs(beff-targetb)
        if count >=500:
            print("Could not reach target in 500 iterations. Quitting.")
            break
        
        # print ("CvsL =",cvsl, "CvsB =",cvsb, "udsg eff = ", leff, "b eff =", beff)
    ceff = getEff('c',tag,cvsl,cvsb)
    return cvsl, cvsb, ceff, beff, leff

# === Just an example result to plot ===
#UL16 PostVFP
# results = {'d': {'Loose': [0.08794335937499997, 0.20418749999999997, 0.8951978206046597], 'Medium': [0.18000195312500208, 0.2211875000000011, 0.5296963919076899], 'Tight': [0.4070312499999934, 0.13628125000000385, 0.2944899964582163]}, 'df': {'Loose': [0.038759765625000005, 0.30453125, 0.9284062594591812], 'Medium': [0.09871484374999896, 0.3525468749999977, 0.588260984189797], 'Tight': [0.26931249999999646, 0.24746874999999965, 0.3373544562119247]}}

# #UL16 PreVFP
# results = {'d': {'Loose': [0.08785449218749995, 0.21440625, 0.8917158872364161], 'Medium': [0.18107812500000198, 0.227687500000001, 0.5225992567294728], 'Tight': [0.41687499999999295, 0.13765625000000636, 0.28483995044863153]}, 'df': {'Loose': [0.038513671875, 0.32731249999999995, 0.9242750374715114], 'Medium': [0.0977929687499992, 0.3700781249999983, 0.5792017712560998], 'Tight': [0.2699999999999961, 0.2559062499999996, 0.3258926432643675]}}
# --------------------------------------

if not args.plotoldresults:
    results = {}
    for pref in taggers:
        tag = "Jet_btag"+pref
        results[tag] = {}
        for wp in mistags:
            cvsl, cvsb, ceff, beff, leff = findcuts(tag,wp)
            print ("Tagger",tag,"WP:",wp, xname, "=",cvsl, yname, "=",cvsb, "c eff =", ceff, "b eff =", beff, "udsg eff =", leff)
            results[tag][wp] = [cvsl,cvsb,ceff]

print("results =",results)

#Plot
ncjets = df[df['truthc']==1].shape[0]
nbjets = df[df['truthb']==1].shape[0]
nljets = df[df['truthudsg']==1].shape[0]
print ("Total c jets:",ncjets)
print ("Total b jets:",nbjets)
print ("Total udsg jets:",nljets)

newseeds = {}

for tag_ in  taggers:
    pref = "Jet_btag"+tag_
    patches = []
    for flav,color in [('b','red'),('udsg','blue'),('c','green')]:
        alljets = df[df['truth'+flav]==1][:ncjets]
        alljets.loc[alljets[pref+xax]>1.] = 1.
        alljets.loc[alljets[pref+yax]>1.] = 1.

        plt.scatter(alljets[pref+xax], alljets[pref+yax],color=color,s=0.1)
        patches.append(mpatches.Patch(color=color, label=flav+' jets'))
    tag = tagname[tag_]
    plt.xlabel(tag + " " + xname)
    plt.ylabel(tag + " " + yname)
    plt.ylim(0.,1.35)
    plt.xlim(0.,1.)
    patches.append(mpatches.Patch(color="white", label='', alpha=0))
    
    hep.cms.label(rlabel="%s (13.6 TeV)"%args.outname)

    ax = plt.gca()
    for wp,c in zip(results[pref],['cyan','#F9FBB2','magenta','#7C6354','#0e3b43']):
        cvsl,cvsb,ceff = results[pref][wp]
        rect = mpatches.Rectangle((cvsl,cvsb), 1-cvsl, 1-cvsb, linewidth=3, linestyle='--', edgecolor=c, facecolor='none')
        if wp == "Medium": t = plt.text(cvsl+0.005,cvsb+0.01,"c eff = %.1f%%"%(ceff*100),color='black',ha='left',va='bottom',size=15)
        else: t = plt.text(cvsl+0.005,cvsb-0.01,"c eff = %.1f%%"%(ceff*100),color='black',ha='left',va='top',size=15)
        t.set_bbox(dict(facecolor=c, alpha=0.8, edgecolor = c))
        ax.add_patch(rect)
        patches.append(Line2D([0], [0], color=c, linewidth=3, linestyle='--',label=wp))

        if wp not in newseeds:
            newseeds[wp] = {}
        newseeds[wp][pref] = {'cvl' : cvsl, 'cvb' : cvsb}
    
    

    plt.legend(handles=patches,loc="upper right",ncol=2)
    plt.text(0.06,1.3,r"t$\bar{t}$ jets"+"\n"+r"$\mathit{p}_{T}>20$ GeV"+"\n"+r"$|\eta|<$2.5",va='top',ha='left',size=20)
    plt.savefig('ctaggerWP_%s_%s.png'%(args.outname,tag_))
    plt.clf()

print()
print("Seeds to hardcode for future:")
print(newseeds)



# Write yaml
def yamlname(st):
    yamldict = {
        "Loose" : "L",
        "Medium": "M",
        "Tight" : "T",
        "Jet_btagDeepFlav" :    "deepJet",
        "Jet_btagPNet"      :   "particleNet",
        "Jet_btagRobustParTAK4":"robustParticleTransformer"
    }
    if st in yamldict: st = yamldict[st]
    return st

yamlfile = f"wp_database_ctagging_{args.outname}.yml"
outfl = open(yamlfile,'w')
for tag in results:
    outfl.write(yamlname(tag)+":\n")
    for wp in results[tag]:
        outfl.write("    "+yamlname(wp)+":\n")
        lst = results[tag][wp]
        outfl.write(f"        eff: {lst[2]*100:.1f}\n")
        outfl.write(f"        CvL: {lst[0]:.3f}\n")
        outfl.write(f"        CvB: {lst[1]:.3f}\n")

outfl.close()

print()
print("Wrote yaml file",yamlfile)
