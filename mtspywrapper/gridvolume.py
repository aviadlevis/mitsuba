import struct, os
import numpy as np
from tempfile import NamedTemporaryFile
from mitsuba.core import Transform, Vector

class GridVolume(object):
    
    def __init__(self):
        """
        This class wraps the gridvolume c++ class implemented in mitsuba (src/volume/gridvolue.cpp). 
        """
        self.__filename = None
        self.__mitsubaparams = None
        self.__scale = None
        self.__boundingBox = None
        self.__ndim = None
        self.__shape = None
            
    def setData(self, volData, boundingBox):
        """
        Generates a binary file (.vol) from a 3d matrix (numpy array)  
        Input: volData - 3D matrix of float representing the voxels values of the object
               boundingBox - bounding box of the object [xmin,ymin,zmin],[xmax,ymax,zmax]
               filename - .vol file name
        """        
        # Scale the data to the interval [0,1]: needed for woodcock integration
        # this scale is later used as a parameters in mitsuba
        self.__scale = volData.max()   
        self.__boundingBox = boundingBox
        self.__ndim = 3
        self.__shape = volData.shape
        
        volData = volData/self.scale
        
        fid = NamedTemporaryFile(delete=False)
        self.setMitsubaParams(fid.name)
        
        fid.write('VOL')                # Bytes 1-3 ASCII Bytes 'V', 'O', and 'L'
        fid.write(struct.pack('B',3))   # Byte 4 File format version number (currently 3)
        fid.write(struct.pack('I',1))   # Bytes 5-8 Encoding identifier (32-bit integer).The following choices are available:
                                        #       1. Dense float32-based representation
                                        #       2. Dense float16-based representation (currently not supported by this implementation)       
                                        #       3. Dense uint8-based representation (The range 0..255 will be mapped to 0..1)
                                        #       4. Dense quantized directions. The directions are stored in spherical coordinates with a total storage cost of 16 bit per entry.
        
                                        
        # Add dimensions to reach a 4D structure (fourth dimention for multi-spectral data)
        for i in range(volData.ndim, 4):
            volData = volData[...,np.newaxis]
            
        # Duplicate dimensions with 1 cell (currently mitsuba accepts only >2 grid points per dimension)
        shape = volData.shape
        dup = [1,1,1,1]
        for i in range(3):
            # Singelton on that dimension - requieres duplication 
            if (shape[i] == 1): 
                self.__ndim -= 1   
                dup[i] = 2
                
        volData = np.tile(volData, dup)
        shape = volData.shape
        ncells = shape[0]*shape[1]*shape[2]
        
        fid.write(struct.pack(4*'I',*shape))        # Bytes 9-24 Number of cells along the X,Y,Z axes (32 bit integer); Bytes 21-24 Number of channels (32 bit integer, supported values: 1 or 3)           
        fid.write(struct.pack(6*'f',*boundingBox))  # Bytes 25-48 Axis-aligned bounding box of the data stored in single precision order: (xmin, ymin, zmin, xmax, ymax, zmax)
    
        # Write the data: Bytes 49-*
        # Binary data of the volume stored in the specified encoding. The data are ordered so that the following C-style indexing operation makes sense
        # after the file has been mapped into memory: data[((zpos*yres + ypos)*xres + xpos)*channels + chan]
        # where (xpos, ypos, zpos, chan) denotes the lookup location.
        fid.write(struct.pack('f'*ncells, *volData.ravel(order='F')));        
        fid.close()
        
    def getData(self):
        """
        Generates 3d matrix (ndarray) from a binary of .vol type
        Output: volData - 3D matrix of float representing the voxels values of the object
                boundingBox - bounding box of the object [xmin,ymin,zmin],[xmax,ymax,zmax]
        """     
        
        fid = open(self.filename)
        
        # Reading first 48 bytes of volFileName as header , count begins from zero  
        header = fid.read(48)  
        
        # Converting header bytes 8-21 to volume size [xsize,ysize,zsize] , type = I : 32 bit integer
        size = struct.unpack(3*'I', bytearray(header[8:20]))
    
        # Converting header bytes 24-47 to bounding box [xmin,ymin,zmin],[xmax,ymax,zmax] type = f : 32 bit float
        boundingBox = struct.unpack(6*'f', bytearray(header[24:48]))
    
        # Converting data bytes 49-* to a 3D matrix size of [xsize,ysize,zsize], 
        # type = f : 32 bit float   
        bindata = fid.read()
        nCells = size[0]*size[1]*size[2]
        volData = np.array(struct.unpack(nCells*'f', bytearray(bindata)))
        volData = volData.reshape(size, order='F')
        fid.close()

        for ax in range(3):
            u_volData, counts =  np.unique(volData, axis=ax, return_counts=True)
            if np.all(counts==2):
                volData = u_volData
            
        volData *= self.scale
        return volData, boundingBox
    
    def setMitsubaParams(self, filename):
        self.__filename = filename
        self.__mitsubaparams = {
            'type' : 'gridvolume',
            'filename' : filename
        }
       
    def getMitsubaParams(self):
        return self.__mitsubaparams
    
    def getWorldTransform(self):
        bb = self.boundingBox
        buttomLeft = Vector(bb[0], bb[1], bb[2])
        topRight = Vector(bb[3], bb[4], bb[5])
        scaleVector = (topRight - buttomLeft)/2.0
        translateVector = buttomLeft + scaleVector
        tform = Transform.translate(translateVector) * \
            Transform.scale(scaleVector)
        return tform     
        
    @property 
    def filename(self):
        return self.__filename 
    
    @property
    def scale(self):
        return self.__scale
    
    @property
    def ndim(self):
        return '{}D'.format(self.__ndim)
    
    @property
    def shape(self):
        return self.__shape
    
    @property 
    def boundingBox(self):
        if (self.__boundingBox is None)&(self.getMitsubaParams() is not None):
            fid = open(self.filename)
            header = fid.read(48)  
            # Converting header bytes 24-47 to bounding box [xmin,ymin,zmin],[xmax,ymax,zmax] type = f : 32 bit float
            boundingBox = struct.unpack(6*'f', bytearray(header[24:48]))
            self.__boundingBox = boundingBox
        return self.__boundingBox
            
            
    def __del__(self):
        os.unlink(self.__filename)
        assert not os.path.exists(self.filename)    