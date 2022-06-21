import uproot
import random
import awkward as ak
#from uproot_methods import TLorentzVectorArray
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys, os
from scipy import constants
import math
from math import pi
import pickle

from awkward.layout import ListOffsetArray64


def convertNBIBToFrac(x):
    return x/2992.

def smear(arr,sigma):
    #Convert it to a 1D numpy array and perform smearing
    numpy_arr = np.asarray(arr.layout.content)
    smeared_arr = np.random.normal(numpy_arr, sigma)
    #Convert it back to awkward form
    return ak.Array(ListOffsetArray64(arr.layout.offsets, ak.Array(smeared_arr).layout))

SpeedOfLight = constants.c/1e6 # mm/ns

mpl.rcParams['patch.force_edgecolor'] = True
mpl.rcParams['patch.linewidth']       = 0.5
mpl.rcParams['savefig.dpi']           = 300
mpl.rcParams['agg.path.chunksize'] = 10000

colors  = ["#A4036F","#5213ba","#048ba8","#16db93","#16f421","#efea5a","#f29e4c","#f96e21","#ff1010","#964b00","#000000","#ffffff"]
markers = ["o","s","D","^","v","<",">","*","X","p"]


# hard scatter only (0 BIB)
noBIBFile = uproot.open(f"DiHiggs_full.root")
noBIBTree = noBIBFile["LCTupleDefault"]

hsMoX = noBIBTree["mcmox"].array()
hsMoY = noBIBTree["mcmoy"].array()
hsMoZ = noBIBTree["mcmoz"].array()
hsMCE = noBIBTree["mcene"].array()
hsStat= noBIBTree["mcgst"].array()

hsMoX_Squared = np.square(hsMoX)
hsMoY_Squared = np.square(hsMoY)
hsMoR = np.sqrt(hsMoX_Squared + hsMoY_Squared)
hsTheta = np.arctan2(hsMoR,hsMoZ)
hsPhi = np.arctan2(hsMoY,hsMoX)

hsfilter = hsStat == 1
hsMoZ = hsMoZ[hsfilter]
hsMoR = hsMoR[hsfilter]
hsMCE = hsMCE[hsfilter]
hsTheta = hsTheta[hsfilter]
hsStat = hsStat[hsfilter]
hsPhi = hsPhi[hsfilter]

hsfilter = hsMCE > 1
hsMoZ = hsMoZ[hsfilter]
hsMoR = hsMoR[hsfilter]
hsMCE = hsMCE[hsfilter]
hsTheta = hsTheta[hsfilter]
hsStat = hsStat[hsfilter]
hsPhi = hsPhi[hsfilter]

hsZComplete = [y for x in hsMoZ for y in x]
hsRComplete = [y for x in hsMoR for y in x]
hsEComplete = [y for x in hsMCE for y in x]
hsThComplete= [y for x in hsTheta for y in x]
hsStComplete= [y for x in hsStat for y in x]
hsPhiComplete=[y for x in hsPhi for y in x]

SignalZ = []
SignalR = []
SignalE = []
SignalTh= []
SignalSt= []
SignalPhi=[]

for i in range(len(hsZComplete)):
    SignalZ.append(hsZComplete[i])
    SignalR.append(hsRComplete[i])
    SignalE.append(hsEComplete[i])
    SignalTh.append(hsThComplete[i])
    SignalSt.append(hsStComplete[i])
    SignalPhi.append(hsPhiComplete[i])
            
pickleoutput = [SignalZ,SignalR,SignalE,SignalTh,SignalPhi,SignalSt]
picklefilename = f"DiHiggs_ZREThPhSt_ECut.pickle"
with open(picklefilename, 'wb') as f:
    pickle.dump(pickleoutput,f)
