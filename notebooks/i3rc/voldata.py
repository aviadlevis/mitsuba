"""
Define the volumetric cloud fields for the experiments. 
Create a dict(): 'voldata' with the following keys
   1. 'case1': 1D step cloud
   2. 'case2':
   3. 'case3':
   4. 'case4':
"""
__all__ = ['voldata']

import numpy as np
from mtspywrapper import GridVolume

voldata = {
    'case1' : GridVolume(),
    'case2' : GridVolume()
}

def initCase1(voldata):
    tauField = np.hstack((np.full(16, 2.00), np.full(16, 18.00)))
    boundingBox = [0, 0, 0, 0.5, 0.5, 0.25]  # [xmin, ymin, zmin, xmax, ymax, zmax] in km units
    geometricalThickness = boundingBox[5] - boundingBox[2]
    betaField = tauField/geometricalThickness
    voldata['case1'].setData(betaField, boundingBox)

def initCase2(voldata):
    betaField = np.loadtxt(fname='i3rc/case2/mmcr_tau_32km_020898') 
    boundingBox = [0, 0, 0, 0.64, 0.5, 0.25]  # [xmin, ymin, zmin, xmax, ymax, zmax] in km units

def initCase3(voldata):
    pass

def initCase4(voldata):
    pass

initCase1(voldata)
initCase2(voldata)
initCase3(voldata)
initCase4(voldata)