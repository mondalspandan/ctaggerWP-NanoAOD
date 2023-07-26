# ctaggerWP-NanoAOD

Dependencies: Python 3, pandas, matplotlib, numpy, mplhep.

## Step 1:
1. Edit `rootToPkl_individualfile.py` to replace "outdir".
2. Run `python rootToPkl_individualfile.py "/path/to/nanoaod/file.root"`. XRootD paths accepted. Put the file path in quotes if using wildcards.

## Step 2:
Run `python ctaggerWP.py -f /path/to/output/pkl/from/step1.pkl`.