import numpy as np
from mtspywrapper import pyMedium
from mitsuba.core import PluginManager, Vector, Point, Transform, Spectrum
from mitsuba.render import Sensor


class pySensor(object):
    
    def __init__(self):
        """"""
        self.__pmgr = PluginManager.getInstance() 
        self.set_sensor_type(None)
        self.set_film(None, None)
        self.set_sampler(None)
    
    def set_sensor_type(self, sensor_type):
        self.__sensor_type = sensor_type
        
    def get_sensor_type(self):
        return self.__sensor_type
    
    def set_film(self, width, height):
        self.__film = {
            'type' : 'mfilm',
            'fileFormat' : 'numpy',
            'width' : width,
            'height' : height
         }
    
    def get_film(self):
        return self.__film
    
    def set_sampler(self, num_samples, sampler_type='ldsampler'):
        self.__sampler = {
            'type' : 'ldsampler',
            'sampleCount' : num_samples
        }
        
    def get_sampler(self):
        return self.__sampler
    
    def get_mitsuba_sensor(self):
        sensor = self.__pmgr.create({
            'type' : self.get_sensor_type(),
            'toWorld' : self.get_world_transform(),
            'film' : self.get_film(),
            'sampler' : self.get_sampler()
        })
        return sensor
    
    
class pyParallelRaySensor(pySensor):
    
    def __init__(self, medium):
        """"""
        super(pyParallelRaySensor, self).__init__()
        self.__medium_bounding_box = medium.bounding_box
        width = medium.shape[0]
        try:
            height = medium.shape[1]
        except IndexError:
            height = 1
        self.set_sensor_type('orthographic')
        self.set_world_transform(0.0, 0.0)
        self.set_film(width, height)

    def set_world_transform(self, view_zenith, view_azimuth):
        """ 
        Set a parallal ray sensor world transform 
        Input: 
          view_zenith [deg]
          view_azimuth [deg]  
        """
        [xmin, ymin, zmin, xmax, ymax, zmax] = self.__medium_bounding_box

        cos_azimuth, sin_azimuth = np.cos(np.deg2rad(view_azimuth)), np.sin(np.deg2rad(view_azimuth))
        cos_zenith, sin_zenith = np.cos(np.deg2rad(view_zenith)), np.sin(np.deg2rad(view_zenith))

        assert(cos_zenith > 0.0), "Error: cos(zenith) < 0.0 ===> zenith > 90.0"
        scale_vector = Vector(
            (xmax/2) * (cos_zenith + (1 - cos_zenith) * np.abs(sin_azimuth)), 
            (ymax/2) * (cos_zenith + (1 - cos_zenith) * np.abs(cos_azimuth)), 
            1.0
        )
        translation_vector = Vector(
            cos_azimuth * sin_zenith, 
            sin_azimuth * sin_zenith, 
            cos_zenith
        )        
        target = Point(xmax/2, ymax/2, zmax)
        origin = target + (xmax**2 + ymax**2) * translation_vector 
        up = Vector(0,1,0)

        self.__world_transform = Transform.lookAt(origin, target, up) * Transform.scale(scale_vector)
        
    def get_world_transform(self):
        return self.__world_transform        