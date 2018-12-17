from struct import pack, unpack
from os import path, unlink
from tempfile import NamedTemporaryFile
from mitsuba.core import PluginManager, Transform, Vector, Spectrum
import numpy as np
from mitsuba.render import Medium


class pyMedium(object):
    
    def __init__(self):
        """
        This class wraps the c++ gridvolume implemented in mitsuba (src/volume/gridvolue.cpp).
        """
        self.__pmgr = PluginManager.getInstance() 
        self.__filename = None
        self.__mitsuba_density = None
        self.set_single_scattering_albedo(1.0)
        self.__phase = None
        self.__scale = None
        self.__bounding_box = None
        self.__ndim = None
        self.__shape = None
        self.set_boundary()
    
    def get_mitsuba_medium(self):
        medium = self.__pmgr.create({
            'type' : 'heterogeneous',
            'method' : 'simpson',
            'xBoundary' : self.__x_boundary,
            'yBoundary' : self.__y_boundary,
            'density' :  self.get_mitsuba_density(),
            'albedo' : self.get_mitsuba_albedo(),
            'phase' : self.get_mitsuba_phase(),
            'scale' : self.scale
        })
        return medium
    
    def get_mitsuba_bounding_box(self):
        bounding_box = self.__pmgr.create({
            'type' : 'cube',
            'toWorld' : self.get_world_transform(), 
            'interior' : self.get_mitsuba_medium()
        })
        return bounding_box
        
        
    def set_boundary(self, x='open', y='open'):
        """ 
        Set boundary conditions to open / periodic
        Input
          x: 'open' / 'periodic' 
          y: 'open' / 'periodic'
        """
        self.__x_boundary = x
        self.__y_boundary = y
    
    def get_boundary(self):
        """ 
        Output boundary conditions of (x, y)
        """
        return (self.__x_boundary, self.__y_boundary)
        
    def set_single_scattering_albedo(self, ssa):
        self.__single_scattering_albedo = Spectrum(ssa)
        
    def get_mitsuba_albedo(self):
        albedo = {'type' : 'constvolume',
                  'value' : self.single_scattering_albedo}
        return albedo
    
    def set_hg_phase(self, g=0.85):
        self.__phase = {'type' : 'hg',
                        'g' : g}
    
    def get_mitsuba_phase(self):
        return self.__phase
    
    def set_density(self, volume, bounding_box):
        """
        Generates a binary file (.vol) from a 3d matrix (numpy array)
        Input 
          volume: 3D matrix of float representing the voxels values of the object
          bounding_box: bounding box of the object [xmin, ymin, zmin, xmax, ymax, zmax]
        """
        # Scale the data to the interval [0,1]: needed for woodcock integration
        # this scale is later used as a parameters in mitsuba
        self.__scale = volume.max()   
        self.__bounding_box = bounding_box
        self.__ndim = 3
        self.__shape = volume.shape
        
        volume = volume/self.scale
        
        fid = NamedTemporaryFile(delete=False)
        self.set_mitsuba_density(fid.name)
        
        fid.write('VOL')                # Bytes 1-3 ASCII Bytes 'V', 'O', and 'L'
        fid.write(pack('B',3))   # Byte 4 File format version number (currently 3)
        fid.write(pack('I',1))   # Bytes 5-8 Encoding identifier (32-bit integer).The following choices are available:
                                        #       1. Dense float32-based representation
                                        #       2. Dense float16-based representation (currently not supported by this implementation)       
                                        #       3. Dense uint8-based representation (The range 0..255 will be mapped to 0..1)
                                        #       4. Dense quantized directions. The directions are stored in spherical coordinates with a total storage cost of 16 bit per entry.
        
                                        
        # Add dimensions to reach a 4D structure (fourth dimention for multi-spectral data)
        for i in range(volume.ndim, 4):
            volume = volume[...,np.newaxis]
            
        # Duplicate dimensions with 1 cell (currently mitsuba accepts only >2 grid points per dimension)
        shape = volume.shape
        dup = [1,1,1,1]
        for i in range(3):
            # Singelton on that dimension - requieres duplication 
            if (shape[i] == 1): 
                self.__ndim -= 1   
                dup[i] = 2
                
        volume = np.tile(volume, dup)
        shape = volume.shape
        ncells = shape[0]*shape[1]*shape[2]
        
        fid.write(pack(4*'I',*shape))        # Bytes 9-24 Number of cells along the X,Y,Z axes (32 bit integer); Bytes 21-24 Number of channels (32 bit integer, supported values: 1 or 3)           
        fid.write(pack(6*'f',*bounding_box)) # Bytes 25-48 Axis-aligned bounding box of the data stored in single precision order: (xmin, ymin, zmin, xmax, ymax, zmax)
    
        # Write the data: Bytes 49-*
        # Binary data of the volume stored in the specified encoding. The data are ordered so that the following C-style indexing operation makes sense
        # after the file has been mapped into memory: data[((zpos*yres + ypos)*xres + xpos)*channels + chan]
        # where (xpos, ypos, zpos, chan) denotes the lookup location.
        fid.write(pack('f'*ncells, *volume.ravel(order='F')));        
        fid.close()
        
    def get_density(self):
        """
        Generates 3D matrix (ndarray) from a binary of .vol type
        Output
          volume: 3D matrix of float representing the voxels values of the object
          bounding_box: bounding box of the object [xmin, ymin, zmin, xmax, ymax, zmax]
        """     
        
        fid = open(self.filename)
        
        # Reading first 48 bytes of volFileName as header , count begins from zero  
        header = fid.read(48)  
        
        # Converting header bytes 8-21 to volume size [xsize,ysize,zsize] , type = I : 32 bit integer
        size = unpack(3*'I', bytearray(header[8:20]))
    
        # Converting header bytes 24-47 to bounding box [xmin,ymin,zmin],[xmax,ymax,zmax] type = f : 32 bit float
        # bounding_box = unpack(6*'f', bytearray(header[24:48]))
    
        # Converting data bytes 49-* to a 3D matrix size of [xsize,ysize,zsize], 
        # type = f : 32 bit float   
        binary_data = fid.read()
        nCells = size[0]*size[1]*size[2]
        volume = np.array(unpack(nCells*'f', bytearray(binary_data)))
        volume = volume.reshape(size, order='F')
        fid.close()

        for ax in range(3):
            u_volume, counts =  np.unique(volume, axis=ax, return_counts=True)
            if np.all(counts==2):
                volume = u_volume
            
        volume *= self.scale
        return volume
    
    def set_mitsuba_density(self, filename):
        self.__filename = filename
        self.__mitsuba_density = {
            'type' : 'gridvolume',
            'filename' : filename
        }
    
    def get_world_transform(self):
        bb = self.bounding_box
        bottom_left = Vector(bb[0], bb[1], bb[2])
        top_right = Vector(bb[3], bb[4], bb[5])
        scale_vector = (top_right - bottom_left)/2.0
        translate_vector = bottom_left + scale_vector
        transform = Transform.translate(translate_vector) * \
            Transform.scale(scale_vector)
        return transform     
        
    def get_mitsuba_density(self):
        return self.__mitsuba_density
    
    @property 
    def filename(self):
        return self.__filename 
    
    @property
    def scale(self):
        return self.__scale
    
    @property
    def single_scattering_albedo(self):
        return self.__single_scattering_albedo  
    
    @property
    def x_boundary(self):
        return self.__x_boundary

    @property
    def y_boundary(self):
        return self.__y_boundary
      
    @property
    def ndim(self):
        return '{}D'.format(self.__ndim)
    
    @property
    def shape(self):
        return self.__shape
    
    @property 
    def bounding_box(self):
        if (self.__bounding_box is None)&(self.__mitsuba_density is not None):
            fid = open(self.filename)
            header = fid.read(48)  
            # Converting header bytes 24-47 to bounding box [xmin,ymin,zmin],[xmax,ymax,zmax] type = f : 32 bit float
            bounding_box = unpack(6*'f', bytearray(header[24:48]))
            self.__bounding_box = bounding_box
        return self.__bounding_box
            
            
    def __del__(self):
        unlink(self.__filename)
        assert not path.exists(self.filename)    