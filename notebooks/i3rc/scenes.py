"""
Define a dictionary conataining mitsuba scenes for each experiment. 
Experiments are as follows:
Exp.#    
        1.    SZA = 0,  w0 = 1
        2.    SZA = 60, w0 = 1
        3.    SZA = 0,  w0 = 0.99
        4.    SZA = 60, w0 = 0.99

"RQ" is the radiative quantity. RQ takes the following values:
        - RQ=R     (reflectance)
        - RQ=T     (transmittance)
        - RQ=A     (absorptance)
        - RQ=H     (net horizontal flux)
        - RQ=Iu    (nadir reflectivity)
        - RQ=I601  (reflectivity at 60 view, 0 azimuth)
        - RQ=I602  (reflectivity at 60 view, 180 azimuth)
        - RQ=Id    (zenith transmissivity) 
"""

__all__ = ['scenes']

experiments = [
    'exp1_Iu',
    'exp2_Iu',
    'exp3_Iu',
    'exp4_Iu'
]

from mitsuba.core import *
from mitsuba.render import Scene
import numpy as np
from voldata import *

scenes = dict()
for experiment in experiments:
    expNum = int(experiment[3])

    # Set single-scattering albedo according to the experiment number
    albedo = Spectrum(1.0)
    if (expNum == 3) or (expNum == 4):
        albedo = Spectrum(0.99)

    # Set solar zenith angle according to the experiment number
    solarZenith = 0.0
    if (expNum == 2) or (expNum == 4):
        solarZenith = np.deg2rad(60.0)
    solarDirection = Vector(0.0, 0.0, -np.cos(solarZenith))

    scene = Scene()
    pmgr = PluginManager.getInstance()
    
    # Create a sensor, film & sample generator
    scene.addChild(pmgr.create({
        'type' : 'radiancemeter',
        'film' : {
            'type' : 'mfilm',
            'fileFormat' : 'numpy'
            },
        'sampler' : {
            'type' : 'ldsampler',
            'sampleCount' : 100
        }
    }))
    
    # Set the integrator
    scene.addChild(pmgr.create({
        'type' : 'volpath_simple',
    }))
    
    # Add heterogeneous medium
    scene.addChild(pmgr.create({
        'type' : 'heterogeneous',
        'method' : 'simpson',
        'density' : voldata['case1'].getMitsubaParams(),
        'albedo' : {
            'type' : 'constvolume',
            'value' : albedo
            },
        'phase' : {
            'type' : 'hg',
            'g' : 0.85
            },
        'scale' : voldata['case1'].scale
    }))
    
    # Create medium bounding box
    scene.addChild(pmgr.create({
        'type' : 'cube',
        'xPeriodic' : True,
        'yPeriodic' : True, 
        'toWorld' : voldata['case1'].getWorldTransform(), 
        'interior' : scene.getMedia()[0]
    }))
    
    scene.addChild(pmgr.create({
        'type' : 'directional',
        'id' : 'Solar',
        'direction' : solarDirection
    }))
    scene.configure()

    scenes[experiment] = scene