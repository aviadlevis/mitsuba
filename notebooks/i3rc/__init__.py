"""
I3RC benchmark dataset
Intercomparison of 3D Radiation Codes (I3RC): https://i3rc.gsfc.nasa.gov/
-------------------------------------------------------------------------------------
Description
I3RC is an ongoing project initiated in the late 1990s. Its goals include: 
- Comparing methods available for 3D atmospheric radiative transfer calculations.
- Providing benchmark results for testing and debugging 3D radiative transfer codes.
- Publishing an open source toolkit (community 3D Monte Carlo code).
- Helping atmospheric science education by creating an archive of illustrative images 
  and other resources on 3D radiative transfer.
-------------------------------------------------------------------------------------
"""

import os, sys
import numpy as np
from mtspywrapper import GridVolume
from mitsuba.core import PluginManager, Vector, Point, Transform, Spectrum
from mitsuba.render import Scene

MTSPATH = next(path for path in sys.path if path.endswith('mitsuba'))
BASEPATH = os.path.join(MTSPATH, 'notebooks','i3rc')

class Case(object):
    """
    Case class is the basic class for importing I3RC results and scenes stored per case. 
    Each case inherits from this base class adding its own parameters and attributes.

    Attributes:
        singleScatteringAlbedo: A Float in the range [0,1].
        solarZenithAngle: A Float in the range [0,180].
        solarAzimuthAngle: A Float in the range [0, 360].
        experiment: A String: exp#_RQ. See inherited methods for valid experiment strings.
        resultFolder: A String.
        volumetricData: GridVolume type. See mtspywrapper.GridVolume.
        xBoundary, yBoundary: 'periodic'/'open'
    """

    def __init__(self):
        self.singleScatteringAlbedo = None
        self.solarZenithAngle = None
        self.solarAzimuthAngle = None
        self.viewZenithAngle = None
        self.viewAzimuthAngle = None
        self.experiment = None
        self.resultFolder = None
        self.volumetricData = GridVolume() 
        self.xBoundary = 'periodic'
        self.yBoundary = 'periodic'
    
    
    def getResults(self):
        """Import consensus results from resultFolder."""
        
        assert(self.resultFolder is not None), "resultFolder is None: choose a case."
        assert(self.experiment is not None),   "Experiment is None: choose an Experiment."
        
        resultFileNames = os.listdir(self.resultFolder)
        fname = [filename for filename in resultFileNames if self.experiment in filename]
        assert(len(fname) == 1), "Two or more experiment files with the same name: {}, {}.".format(fname[0], fname[1])
        
        fpath = os.path.join(self.resultFolder, fname[0])
        results = np.loadtxt(fpath)
        return results
        
        
    def setExperiment(self, experiment):
        """
        Set the experiment parameter.         
        Set sensor parameters according to RQ (radiative quantity). RQ takes the following values:
            - RQ=R     (reflectance)
            - RQ=T     (transmittance)
            - RQ=A     (absorptance)
            - RQ=H     (net horizontal flux)
            - RQ=Iu    (nadir reflectivity)
            - RQ=I601  (reflectivity at 60 view, 0 azimuth)
            - RQ=I602  (reflectivity at 60 view, 180 azimuth)
            - RQ=Id    (zenith transmissivity)
        """
        self.experiment = experiment
        
        RQ = experiment[5:]
        if (RQ == 'Iu'):
            self.viewZenithAngle, self.viewAzimuthAngle = 0.0, 0.0
        elif (RQ == 'I601'):
            self.viewZenithAngle, self.viewAzimuthAngle = 60.0, 0.0
        elif (RQ == 'I602'):
            self.viewZenithAngle, self.viewAzimuthAngle = 60.0, 180.0
        
              
    def getScene(self, nSamples):
        """Import the Mitsuba scene according to the
           Case and Experiment. nSamples are the number of "photons" to use """
        
        assert(self.experiment is not None), "Experiment is not set: Use setExperiment(EXP)."
        volData = self.volumetricData
        scene = Scene()
        pmgr = PluginManager.getInstance()  
        
        width = volData.shape[0]
        try:
            height = volData.shape[1]
        except IndexError:
            height = 1
            
        # Create a sensor, film & sample generator
        scene.addChild(pmgr.create({
            'type' : 'orthographic',
            'toWorld' : self.sensorTransformation,
            'film' : {
                'type' : 'mfilm',
                'fileFormat' : 'numpy',
                'width' : width,
                'height' : height
             },
            'sampler' : {
                'type' : 'ldsampler',
                'sampleCount' : nSamples
            }
        }))
        # Set the integrator
        scene.addChild(pmgr.create({
            'type' : 'volpath_simple',
            'rrDepth' : 30
        }))
        # Add heterogeneous medium
        scene.addChild(pmgr.create({
            'type' : 'heterogeneous',
            'method' : 'simpson',
            'xBoundary' : self.xBoundary,
            'yBoundary' : self.yBoundary,
            'density' : volData.getMitsubaParams(),
            'albedo' : {
                'type' : 'constvolume',
                'value' : self.singleScatteringAlbedo
             },
            'phase' : {
                'type' : 'hg',
                'g' : 0.85
             },
            'scale' : volData.scale
        }))
        # Create medium bounding box
        scene.addChild(pmgr.create({
            'type' : 'cube',
            'toWorld' : volData.getWorldTransform(), 
            'interior' : scene.getMedia()[0]
        }))
        # Direction Solar source 
        scene.addChild(pmgr.create({
            'type' : 'directional',
            'id' : 'Solar',
            'direction' : self.solarPhotonDirection
        }))
        scene.configure()
        return scene

    
    @property
    def solarPhotonDirection(self):
        """ Changing Solar direction to photon flow direction """
        zenith  = np.deg2rad(self.solarZenithAngle + 180.0)
        azimuth = np.deg2rad(self.solarAzimuthAngle + 180.0)
        solarDirection = Vector(
            np.sin(zenith) * np.cos(azimuth), 
            np.sin(zenith) * np.sin(azimuth), 
            np.cos(zenith)
        )
        return solarDirection
    
    
    @property
    def sensorTransformation(self):
        """ 
        Return sensor transformation using the viewZenithAngle[deg], viewAzimuthAngle[deg] and 
        domain size: volumetricData.boundingBox. 
        Transofrmation is expressed as a translation (T1) scale (T2) and rotation (T3) around an axis (V) 
        """
        [xmin, ymin, zmin, xmax, ymax, zmax] = self.volumetricData.boundingBox
        
        # T1*T2: Move to the top of the domain (twice the height). Origin is (0,0) and scale to cover the entire domain
        cosAzimuth, sinAzimuth = np.cos(np.deg2rad(self.viewAzimuthAngle)), np.sin(np.deg2rad(self.viewAzimuthAngle))
        cosZenith, sinZenith = np.cos(np.deg2rad(self.viewZenithAngle)), np.sin(np.deg2rad(self.viewZenithAngle))
        
        assert(cosZenith > 0.0), "Error: cosZenith < 0.0 ===> Zenith > 90.0"
        scaleVector = Vector(
            (xmax/2) * (cosZenith + (1 - cosZenith) * np.abs(sinAzimuth)), 
            (ymax/2) * (cosZenith + (1 - cosZenith) * np.abs(cosAzimuth)), 
            1.0
        )
        translationVector = Vector(
            cosAzimuth * sinZenith, 
            sinAzimuth * sinZenith, 
            cosZenith
        )        
        target = Point(xmax/2, ymax/2, zmax)
        origin = target + (xmax**2 + ymax**2) * translationVector 
        up = Vector(0,1,0)
        
        T = Transform.lookAt(origin, target, up) * Transform.scale(scaleVector)
        return T 
    
    
    @property
    def radiance2BRF(self):
        return np.pi/np.cos(np.deg2rad(self.solarZenithAngle))
        
        
#----------------------------------------------------------------------------------------------------#
#                                Case1: 1D Step Cloud                                                #
#----------------------------------------------------------------------------------------------------#

class Case1(Case):
    """ 
    Case1: 1D step cloud (optical path goes from tau=2 to tau=18)
        - Project page: https://i3rc.gsfc.nasa.gov/input/step_cloud/index.html
        - Readme: http://i3rc.gsfc.nasa.gov/input/step_cloud/README.txt
    """
    
    def __init__(self):
        super(Case1, self).__init__()
        self.resultFolder = os.path.join(BASEPATH, 'case1','consensus_results')
        
        # Create the volumetric data for the step cloud case
        tauField = np.hstack((np.full(16, 2.00), np.full(16, 18.00)))
        boundingBox = [0, 0, 0, 0.5, 0.5, 0.25]   # [xmin, ymin, zmin, xmax, ymax, zmax] in km units
        geometricalThickness = boundingBox[5] - boundingBox[2]
        betaField = tauField/geometricalThickness
        self.volumetricData.setData(betaField, boundingBox)
        
        
    def setExperiment(self, experiment):
        """
        SZA - Solar Zenith Angle
        w0 - Single Scattering Albedo
        
        Experiments are as follows:
        Exp.#    
            1.    SZA = 0,  w0 = 1
            2.    SZA = 60, w0 = 1
            3.    SZA = 0,  w0 = 0.99
            4.    SZA = 60, w0 = 0.99
        """
        super(Case1, self).setExperiment(experiment)
       
        self.solarAzimuthAngle = 0.0
        exp = self.experiment[:4]
        if (exp == 'exp1'):
            self.solarZenithAngle, self.singleScatteringAlbedo = 0.0, Spectrum(1.0)
        elif (exp == 'exp2'):
            self.solarZenithAngle, self.singleScatteringAlbedo = 60.0, Spectrum(1.0)
        elif (exp == 'exp3'):
            self.solarZenithAngle, self.singleScatteringAlbedo = 0.0, Spectrum(0.99)
        elif (exp == 'exp4'):
            self.solarZenithAngle, self.singleScatteringAlbedo = 60.0, Spectrum(0.99)
        else:
            raise "Error: {} is not a valid experiment for Case1.".format(exp)
        
        
#----------------------------------------------------------------------------------------------------#
#                         Case2: 2D field derived from the ARM cloud radar                           #
#----------------------------------------------------------------------------------------------------#

class Case2(Case):
    """ 
    Case2: 2D field derived from the ARM cloud radar 
        - Project page: https://i3rc.gsfc.nasa.gov/input/MMCR/high_res/020898/index.html
        - Readme: https://i3rc.gsfc.nasa.gov/input/MMCR/high_res/020898/README.txt
    """
    
    def __init__(self):
        super(Case2, self).__init__()
        self.resultFolder = os.path.join(BASEPATH, 'case2','consensus_results')
        
        # Load the volumetric data 
        fileName = 'mmcr_tau_32km_020898'
        betaField = np.loadtxt(os.path.join(BASEPATH, 'case2', fileName)).T.reshape((640,1,54))
        boundingBox = [0, 0, 0, 32, 0.5, 54]   # [xmin, ymin, zmin, xmax, ymax, zmax] in km units
        self.volumetricData.setData(betaField, boundingBox)
        
        
    def setExperiment(self, experiment):
        """
        SZA - Solar Zenith Angle
        w0 - Single Scattering Albedo
        As - Albedo Surface
        PF - Phase Function 
        
        Experiments are as follows:
        Exp.#    
            1.    SZA = 0,  w0 = 1
            2.    SZA = 60, w0 = 1
            3.    SZA = 0,  w0 = 0.99
            4.    SZA = 60, w0 = 0.99
            5.    SZA = 60, w0 = 1, As = 0.4, PF = HG
            6:    SZA = 0,  w0 = 1, As = 0, PF = C1
            7:    SZA = 60, w0 = 1, As = 0, PF = C1
            8:    SZA = 60, w0 = 1, As = 0.4, PF = C1
        """
        super(Case2, self).setExperiment(experiment)
       
        self.solarAzimuthAngle = 0.0
        exp = self.experiment[:4]
        if (exp == 'exp1'):
            self.solarZenithAngle, self.singleScatteringAlbedo = 0.0, Spectrum(1.0)
        elif (exp == 'exp2'):
            self.solarZenithAngle, self.singleScatteringAlbedo = 60.0, Spectrum(1.0)
        elif (exp == 'exp3'):
            self.solarZenithAngle, self.singleScatteringAlbedo = 0.0, Spectrum(0.99)
        elif (exp == 'exp4'):
            self.solarZenithAngle, self.singleScatteringAlbedo = 60.0, Spectrum(0.99)
        else:
            raise "Error: {} is not a valid experiment for Case1.".format(exp)
        
        
#----------------------------------------------------------------------------------------------------#
#                                   Case3: Landsat cloud case                                        #
#----------------------------------------------------------------------------------------------------#

class Case3(Case):
    """ 
    Case3: Landsat cloud case 
        - Project page: https://i3rc.gsfc.nasa.gov/input/Landsat/index.html
        - Readme: https://i3rc.gsfc.nasa.gov/input/Landsat/README.txt
    """
    
    def __init__(self):
        super(Case3, self).__init__()
        self.resultFolder = os.path.join(BASEPATH, 'case3','consensus_results')
        
        # Load the volumetric data 
        geometricalThickness = np.loadtxt(os.path.join(BASEPATH, 'case3', 'scene43.dz.128x128'))
        tauField = np.loadtxt(os.path.join(BASEPATH, 'case3', 'scene43.tau.128x128'))
        
        nx, ny, nz = 128, 128, 128
        dz = geometricalThickness.max()/nz
        nCellsOccupied = np.round(geometricalThickness/dz).astype(np.int)
        roundedGeometricalThickness = nCellsOccupied*dz
        betaField2D = np.divide(tauField,
                                roundedGeometricalThickness, 
                                out=np.zeros_like(tauField), 
                                where=roundedGeometricalThickness!=0)
        betaField3D = np.zeros(shape=(nx, ny, nz), dtype=np.float32)
        for x in range(nx):
            for y in range(ny):
                betaField3D[x,y,:nCellsOccupied[x,y]] = betaField2D[x,y]

        boundingBox = [0, 0, 0, 3.84, 3.84, geometricalThickness.max()]   # [xmin, ymin, zmin, xmax, ymax, zmax] in km units
        self.volumetricData.setData(betaField, boundingBox)
        
        
    def setExperiment(self, experiment):
        """
        SZA - Solar Zenith Angle
        w0 - Single Scattering Albedo
  
        Experiments are as follows:
        Exp.#    
            1.    SZA = 0,  w0 = 1
            2.    SZA = 60, w0 = 1
            3.    SZA = 0,  w0 = 0.99
            4.    SZA = 60, w0 = 0.99
        """
        super(Case2, self).setExperiment(experiment)
       
        self.solarAzimuthAngle = 0.0
        exp = self.experiment[:4]
        if (exp == 'exp1'):
            self.solarZenithAngle, self.singleScatteringAlbedo = 0.0, Spectrum(1.0)
        elif (exp == 'exp2'):
            self.solarZenithAngle, self.singleScatteringAlbedo = 60.0, Spectrum(1.0)
        elif (exp == 'exp3'):
            self.solarZenithAngle, self.singleScatteringAlbedo = 0.0, Spectrum(0.99)
        elif (exp == 'exp4'):
            self.solarZenithAngle, self.singleScatteringAlbedo = 60.0, Spectrum(0.99)
        else:
            raise "Error: {} is not a valid experiment for Case1.".format(exp)
        
        
         