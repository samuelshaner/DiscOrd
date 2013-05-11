import numpy as np

class Cell(object):
    
    def __init__(self, num_angles, id):
        self.ang_flux = np.zeros(num_angles)
        self.flux     = 0.0
        self.old_flux = 0.0
        self.id       = id

    def setMaterial(self, material):
    
        self.material = material
