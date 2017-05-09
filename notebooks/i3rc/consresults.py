"""
Load the consensus results for all test cases. 
Each test case results is in the directory:
       MTSPATH/notebooks/i3rc/consensus_result_case#    where #=1,2,3,4
A readme.txt file is in every folder describing the format of the resutls.
"""

__all__ =  ['consresults']


import os, sys
import numpy as np

consresults = {
    'case1' : dict()
}

MTSPATH = next(path for path in sys.path if path.endswith('mitsuba'))

for case in consresults.iterkeys():
    
    RESULTPATH = os.path.join(MTSPATH, 
                              'notebooks', 
                              'i3rc', 
                              'consensus_result_' + case)
    resultFileNames = os.listdir(RESULTPATH)
    resultFileNames.remove('readme.txt')
    
    for fname in resultFileNames:
        fpath = os.path.join(RESULTPATH, fname)
        key = fname[fname.find('exp'):-4] 
        consresults[case][key] = np.loadtxt(fpath)