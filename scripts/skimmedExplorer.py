import root_numpy as rtnp
from collections import Counter

fileLoc = '/storage/af/group/rdst_analysis/BPhysics/data/cmsMC/CP_BdToDstDs_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/ntuples_B2DstMu/skimmed'
arr = rtnp.root2array(fileLoc+'/Low_bare.root', branches=['MC_CharmedDstSisPdgId'])

print Counter(arr['MC_CharmedDstSisPdgId'])
