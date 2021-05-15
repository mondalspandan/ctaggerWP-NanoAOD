#  Coded in Python 3
#  c-tagger WP measurement using NanoAOD.
#  Author: Spandan Mondal,
#  RWTH Aachen University

import pickle, glob, sys, pandas as pd, argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import mplhep as hep, numpy as np
hep.set_style(hep.style.CMS)

parser = argparse.ArgumentParser("Measurement of c-tagger WPs")
parser.add_argument('-i','--inputdir',type=str,default="/mnt/c/Work/ctaggerWP/",help="Path to input dir containing pkls (if multiple).")
parser.add_argument('-f','--inputfile',type=str,default="/mnt/c/Work/ctaggerWP/UL16NanoAODAPVv2_TTToHadronic.pkl",help="Input file name.")
parser.add_argument('-s','--skipseeds',action="store_true",default=False,help="Skip 1D scans that identify starting points for 2D scans.")
parser.add_argument('-p','--plotoldresults',action="store_true",default=False,help="Only plot stored results, skip everything else.")
parser.add_argument('-o','--outname',type=str,default="preVFP")
parser.add_argument('-n','--njets',type=str,default="-1",help="Maximum number of jets to process.")
args = parser.parse_args()
print (args)

if args.inputfile == "":
    indir = args.inputdir
    filelist = glob.glob(indir+'/*.pkl')
else:
    filelist = [args.inputfile]

df = pd.DataFrame()

for fl in filelist:
    thisdf = pd.read_pickle(fl)
    df = pd.concat([df,thisdf])

def cleandf(tdf):
    cdf = tdf[(tdf['Jet_btagDeepB'] < 1) & (tdf['Jet_btagDeepC'] < 1)
             & (tdf['Jet_btagDeepB'] > 0) & (tdf['Jet_btagDeepC'] > 0) ]
    cdf = cdf[ ((cdf['Jet_puId'] > 0) | (cdf['Jet_pt'] > 50)  ) & (cdf['Jet_pt'] > 20) & (cdf['Jet_jetId'] >= 5) & (abs(cdf['Jet_eta']) < 2.4)]
    njets = int(args.njets)
    if njets > 0: cdf = cdf[:njets]    
    return cdf

def getTaggers(df):
    df['truthb'] = (df['Jet_hadronFlavour'] == 5).astype(int)
    df['truthc'] = (df['Jet_hadronFlavour'] == 4).astype(int)
    df['truthudsg'] = (df['Jet_hadronFlavour'] < 4).astype(int)
    
    df['dCvL'] = df['Jet_btagDeepC']/(1 - df['Jet_btagDeepB'])
    df['dCvB'] = df['Jet_btagDeepC']/(df['Jet_btagDeepC'] + df['Jet_btagDeepB'])
    df['dBvL'] = df['Jet_btagDeepB']/(1 - df['Jet_btagDeepC'])
    df['dBvC'] = df['Jet_btagDeepB']/(df['Jet_btagDeepC'] + df['Jet_btagDeepB'])

    df['dfCvL'] = df['Jet_btagDeepFlavC']/(1 - df['Jet_btagDeepFlavB'])
    df['dfCvB'] = df['Jet_btagDeepFlavC']/(df['Jet_btagDeepFlavC'] + df['Jet_btagDeepFlavB'])
    df['dfBvL'] = df['Jet_btagDeepFlavB']/(1 - df['Jet_btagDeepFlavC'])
    df['dfBvC'] = df['Jet_btagDeepFlavB']/(df['Jet_btagDeepFlavC'] + df['Jet_btagDeepFlavB'])

    return df

df = cleandf(df)
df = getTaggers(df)

def getEff(flav,pref,cvsl=0.,cvsb=0.):
    alljets = df[df['truth'+flav]==1]
    njets = alljets.shape[0]
    passjets = alljets[(alljets[pref+'CvL']>=cvsl) & (alljets[pref+'CvB']>=cvsb)]
    return float(passjets.shape[0])/njets


#Set these, everything else *should* be automatic
mistags = { 
    'Loose':  {'b':.35,'udsg':.9},
    'Medium': {'b':.25,'udsg':.25},
    'Tight':  {'b':.20,'udsg':.03},
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
seeds = {'Loose': {'d': {'cvl': 0.1, 'cvb': 0.22}, 'df': {'cvl': 0.05, 'cvb': 0.32}}, 'Medium': {'d': {'cvl': 0.19, 'cvb': 0.35000000000000003}, 'df': {'cvl': 0.11, 'cvb': 0.49}}, 'Tight': {'d': {'cvl': 0.43, 'cvb': 0.43}, 'df': {'cvl': 0.29, 'cvb': 0.58}}}


#1D scans
if not args.skipseeds or not args.plotoldresults:
    seeds = {}
    for wp in mistags:
        seeds[wp] = {}
        for tag in ['d','df']:
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
accuracy = 1e-5
def findcuts(tag,wp):
    cvsl = seeds[wp][tag]['cvl']
    cvsb = seeds[wp][tag]['cvb']
    targetl = mistags[wp]['udsg']
    targetb = mistags[wp]['b']

    stepcvl = 1e-3
    stepcvb = 1e-3
    leff = getEff('udsg',tag,cvsl,cvsb)
    beff = getEff('b',tag,cvsl,cvsb)
    if leff < targetl: leffwasless = True
    else: leffwasless = False
    if abs(beff-targetb)>accuracy: beffwasless = True
    else: beffwasless = False

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
        
        # print ("CvsL =",cvsl, "CvsB =",cvsb, "udsg eff = ", leff, "b eff =", beff)
    ceff = getEff('c',tag,cvsl,cvsb)
    return cvsl, cvsb, ceff, beff, leff

# Just an example result to plot

#UL16 PostVFP
results = {'d': {'Loose': [0.08794335937499997, 0.20418749999999997, 0.8951978206046597], 'Medium': [0.18000195312500208, 0.2211875000000011, 0.5296963919076899], 'Tight': [0.4070312499999934, 0.13628125000000385, 0.2944899964582163]}, 'df': {'Loose': [0.038759765625000005, 0.30453125, 0.9284062594591812], 'Medium': [0.09871484374999896, 0.3525468749999977, 0.588260984189797], 'Tight': [0.26931249999999646, 0.24746874999999965, 0.3373544562119247]}}

#UL16 PreVFP
results = {'d': {'Loose': [0.08785449218749995, 0.21440625, 0.8917158872364161], 'Medium': [0.18107812500000198, 0.227687500000001, 0.5225992567294728], 'Tight': [0.41687499999999295, 0.13765625000000636, 0.28483995044863153]}, 'df': {'Loose': [0.038513671875, 0.32731249999999995, 0.9242750374715114], 'Medium': [0.0977929687499992, 0.3700781249999983, 0.5792017712560998], 'Tight': [0.2699999999999961, 0.2559062499999996, 0.3258926432643675]}}

if not args.plotoldresults:
    results = {}
    for tag in ['d','df']:
        results[tag] = {}
        for wp in mistags:
            cvsl, cvsb, ceff, beff, leff = findcuts(tag,wp)
            print ("Tagger",tag,"WP:",wp, "CvsL =",cvsl, "CvsB =",cvsb, "c eff =", ceff, "b eff =", beff, "udsg eff =", leff)
            results[tag][wp] = [cvsl,cvsb,ceff]

print("results =",results)

#Plot
ncjets = df[df['truthc']==1].shape[0]
print ("Total c jets:",ncjets)

for pref in ['d','df']:
    patches = []
    for flav,color in [('b','red'),('udsg','blue'),('c','green')]:
        alljets = df[df['truth'+flav]==1][:ncjets]
        alljets.loc[alljets[pref+'CvL']>1.] = 1.
        alljets.loc[alljets[pref+'CvB']>1.] = 1.

        plt.scatter(alljets[pref+'CvL'], alljets[pref+'CvB'],color=color,s=0.1)
        patches.append(mpatches.Patch(color=color, label=flav+' jets'))
    tag = "DeepCSV"
    if pref == "df": tag = "DeepJet"
    plt.xlabel(tag + " CvsL")
    plt.ylabel(tag + " CvsB")
    plt.ylim(0.,1.35)
    plt.xlim(0.,1.)
    
    hep.cms.label(year="2016UL")

    ax = plt.gca()
    for wp,c in zip(results[pref],['cyan','yellow','magenta']):
        cvsl,cvsb,ceff = results[pref][wp]
        rect = mpatches.Rectangle((cvsl,cvsb), 1-cvsl, 1-cvsb, linewidth=3, linestyle='--', edgecolor=c, facecolor='none')
        if wp == "Medium": t = plt.text(cvsl+0.005,cvsb+0.01,"c eff = %.1f%%"%(ceff*100),color='black',ha='left',va='bottom',size=15)
        else: t = plt.text(cvsl+0.005,cvsb-0.01,"c eff = %.1f%%"%(ceff*100),color='black',ha='left',va='top',size=15)
        t.set_bbox(dict(facecolor=c, alpha=0.8, edgecolor = c))
        ax.add_patch(rect)
        patches.append(Line2D([0], [0], color=c, linewidth=3, linestyle='--',label=wp))

    plt.legend(handles=patches,loc="upper right",ncol=2)
    plt.text(0.06,1.3,r"t$\bar{t}$ jets"+"\n"+r"$\mathit{p}_{T}>20$ GeV"+"\n"+r"$|\eta|<$2.4",va='top',ha='left',size=20)
    plt.savefig('ctaggerWP_%s_%s.png'%(args.outname,pref))
    plt.clf()