import numpy as np


class Cell(object):

    def __init__(self, num_angles):

		self.ang_flux     = np.zeros(num_angles)
		self.flux         = 0.0
		self.old_flux     = 0.0
		self.source       = np.zeros(num_angles)

    def setMaterial(self, material):
    
        self.material = material
