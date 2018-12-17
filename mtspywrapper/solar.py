import numpy as np
from mitsuba.render import Emitter
from mitsuba.core import PluginManager, Vector

class pySolarEmitter(object):
    
    def __init__(self, zenith=0.0, azimuth=0.0):
        """"""
        self.__pmgr = PluginManager.getInstance() 
        self.set_solar_angles(zenith, azimuth)
        
    def set_solar_angles(self, zenith, azimuth):
        self.__zenith = zenith
        self.__azimuth = azimuth
        self.__direction = self.angles_to_photon_direction(zenith, azimuth)
    
    def angles_to_photon_direction(self, zenith, azimuth):
        """ 
        Changing Solar angles to photon flow direction
        Inputs
            zenith [deg]
            azimuth [deg]
            
        Output
            photon_direction [Vector]
        """
        zenith  = np.deg2rad(zenith + 180.0)
        azimuth = np.deg2rad(azimuth + 180.0)
        photon_direction = Vector(
            np.sin(zenith) * np.cos(azimuth), 
            np.sin(zenith) * np.sin(azimuth), 
            np.cos(zenith)
        )
        return photon_direction 
    
    def get_mitsuba_emitter(self):
        emitter = self.__pmgr.create({
            'type' : 'directional',
            'id' : 'Solar',
            'direction' : self.__direction
        })
        return emitter
    
    @property
    def brf_factor(self):
        return np.pi/np.cos(np.deg2rad(self.__zenith))
    
    @property
    def zenith(self):
        return self.__zenith
    
    @property
    def azimuth(self):
        return self.__azimuth    
    
    