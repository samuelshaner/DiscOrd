import numpy as np
import matplotlib.pyplot as plt
import quadrature
from math import *
from cell import *
from material import *
import plotter as pttr
import os
import sys
import getopt
import time

class Mesh(object):

	def __init__(self, mesh_size, order, iterations, tolerance, pitch=1.26, fuel_diameter=0.8):
		# create lists for storing values
		self.width  = int(pitch/mesh_size)
		self.quad          = quadrature.LevelSymmetricQuadrature().getQuadrature(order)
		self.fuel_diameter = fuel_diameter
		self.mesh_size     = mesh_size
		self.order         = order
		self.tol			   = tolerance
		self.iterations    = iterations
		self.num_angles    = self.quad['num_angles']
		self.pitch         = pitch
		self.cells         = np.empty(self.width**2, dtype=object)	
		self.ang_flux      = np.zeros(4*self.width*self.num_angles/2))
		self.makeCells()

	def setFuel(self, fuel):

		self.fuel = fuel

	def setModerator(self, moderator):
	
		self.moderator = moderator

	def makeMeshMaterials(self):

		cw = self.width

		width_to_fuel = cw / 2 - int(self.fuel_diameter / self.mesh_size / 2)

		for y in range(cw):
			for x in range(cw):
				if x < width_to_fuel or x > cw - width_to_fuel:
					self.cells[y*cw+x].setMaterial(self.moderator)
				else:
					if y < width_to_fuel or y > self.width - width_to_fuel:
						self.cells[y*cw+x].setMaterial(self.moderator)
					else:
						self.cells[y*cw+x].setMaterial(self.fuel)

	def initializeSource(self):

		self.source = np.zeros((self.width**2, self.num_angles))
		for y in range(self.width):
			for x in range(self.width):
				if self.cells[y*self.width+x].material.mat_type == 'fuel':
					for i in range(self.num_angles):
						self.cells[y*self.width+x].source[i] = 1.0


	def makeCells(self):

		for y in range(self.width):
			for x in range(self.width):
				self.cells[y*self.width+x] = Cell(self.num_angles)

	def solveSn(self):

		cw = self.width
		na = self.num_angles
		eps = 1.0
	
		# loop over iterations
		for iteration in range(self.iterations):

			print 'Sn iteration ' + str(iteration) + ' eps ' + str(eps)

			# sweep in the positive mu and eta direction (bottom left corner)
			for y in range(cw):
				for x in range(cw):

					cell = self.cells[y*cw+x]

					for angle in range(na/4):

						# compute cell centered flux
						cell.ang_flux[angle] = (cell.source[angle] + 2 * self.quad['mu'][angle]/self.mesh_size * self.ang_flux[y*na/2 + angle] + 2 * self.quad['eta'][angle]/self.mesh_size * self.ang_flux[1*cw*na/2 + x*na/2 + angle]) / (cell.material.sigma_t + 2 * self.quad['mu'][angle]/self.mesh_size + 2 * self.quad['eta'][angle]/self.mesh_size)

						# sweep from left to right
						self.ang_flux[y*na/2 + angle] = 2 * cell.ang_flux[angle] - self.ang_flux[y*na/2 + angle]

						# sweep from bottom to top
						self.ang_flux[1*cw*na/2 + x*na/2 + angle] = 2 * cell.ang_flux[angle] - self.ang_flux[1*cw*na/2 + x*na/2 + angle]


			# reflect on right boundary
			for y in range(cw):			
				for angle in range(na/4):
					self.ang_flux[2*cw*na/2 + y*na/2 + angle] = self.ang_flux[y*na/2 + angle]
				
			# reflect flux on top boundary
			for x in range(cw):			
				for angle in range(na/4):
					self.ang_flux[3*cw*na/2 + x*na/2 + angle] = self.ang_flux[1*cw*na/2 + x*na/2 + angle]


			# sweep in the negative mu and positive eta direction (bottom right corner)
			for y in range(cw):
				for x in reversed(range(cw)):

					cell = self.cells[y*cw+x]

					for angle in range(na/4):	

						# compute cell centered flux
						cell.ang_flux[na/4 + angle] = (cell.source[angle] + 2 * self.quad['mu'][angle]/self.mesh_size * self.ang_flux[2, y, angle] + 2 * self.quad['eta'][angle]/self.mesh_size * self.ang_flux[1, x, angle + na/4]) / (cell.material.sigma_t + 2 * self.quad['mu'][angle]/self.mesh_size + 2 * self.quad['eta'][angle]/self.mesh_size)

						# sweep from right to left
						self.ang_flux[2, y, angle] = 2 * cell.ang_flux[na/4 + angle] - self.ang_flux[2, y, angle]

						# sweep from bottom to top
						self.ang_flux[1, x, angle + na/4] = 2 * cell.ang_flux[na/4 + angle] - self.ang_flux[1, x, angle + na/4]


			# reflect on left boundary
			for y in range(cw):			
				for angle in range(na/4):
					self.ang_flux[0,y,angle] = self.ang_flux[2,y,angle]
				
			# reflect flux on top boundary
			for x in range(cw):			
				for angle in range(na/4):
					self.ang_flux[3,x,angle + na/4] = self.ang_flux[1,x,angle + na/4]


			# sweep in the negative mu and eta direction (top right corner)
			for y in reversed(range(cw)):
				for x in reversed(range(cw)):

					cell = self.cells[y*cw+x]

					for angle in range(na/4):	

						# compute cell centered flux
						cell.ang_flux[na/2 + angle] = (cell.source[angle] + 2 * self.quad['mu'][angle]/self.mesh_size * self.ang_flux[2, y, angle + na/4] + 2 * self.quad['eta'][angle]/self.mesh_size * self.ang_flux[3, x, angle + na/4]) / (cell.material.sigma_t + 2 * self.quad['mu'][angle]/self.mesh_size + 2 * self.quad['eta'][angle]/self.mesh_size)

						# sweep from right to left
						self.ang_flux[2, y, angle + na/4] = 2 * cell.ang_flux[na/2 + angle] - self.ang_flux[2, y, angle + na/4]

						# sweep from top to bottom
						self.ang_flux[3, x, angle + na/4] = 2 * cell.ang_flux[na/2 + angle] - self.ang_flux[3, x, angle + na/4]


			# reflect on left boundary
			for y in range(cw):			
				for angle in range(na/4):
					self.ang_flux[0,y,angle + na/4] = self.ang_flux[2,y,angle + na/4]
				
			# reflect flux on bottom boundary
			for x in range(cw):			
				for angle in range(na/4):
					self.ang_flux[1,x,angle + na/4] = self.ang_flux[3,x,angle + na/4]


			# sweep in the positive mu and negative eta direction (top left corner)
			for y in reversed(range(cw)):
				for x in range(cw):

					cell = self.cells[y*cw+x]

					for angle in range(na/4):										
						
						# compute cell centered flux
						cell.ang_flux[3*na/4+ angle] = (cell.source[angle] + 2 * self.quad['mu'][angle]/self.mesh_size * self.ang_flux[0, y, angle + na/4] + 2 * self.quad['eta'][angle]/self.mesh_size * self.ang_flux[3, x, angle]) / (cell.material.sigma_t + 2 * self.quad['mu'][angle]/self.mesh_size + 2 * self.quad['eta'][angle]/self.mesh_size)

						# sweep from left to right
						self.ang_flux[0, y, angle + na/4] = 2 * cell.ang_flux[3*na/4+ angle] - self.ang_flux[0, y, angle + na/4]

						# sweep from bottom to top
						self.ang_flux[3, x, angle] = 2 * cell.ang_flux[3*na/4+ angle] - self.ang_flux[3, x, angle]

			# reflect on right boundary
			for y in range(cw):			
				for angle in range(na/4):
					self.ang_flux[2,y,angle + na/4] = self.ang_flux[0,y,angle + na/4]
				
			# reflect flux on bottom boundary
			for x in range(cw):			
				for angle in range(na/4):
					self.ang_flux[1,x,angle] = self.ang_flux[3,x,angle]


			self.computeScalarFlux()	
			pttr.plotScalarFlux(self, self.order, iteration)
			eps = self.computeL2Norm()
			
			if (eps < self.tol):
				break



	def computeL2Norm(self):

		eps = 0.0

		for y in range(self.width):
			for x in range(self.width):
				eps += (self.cells[y*self.width+x].flux - self.cells[y*self.width+x].old_flux)**2

		sqrt(eps)

		return eps
	

	def computeScalarFlux(self):
		
		# set old flux
		for y in range(self.width):
			for x in range(self.width):
				self.cells[y*self.width+x].old_flux = self.cells[y*self.width+x].flux
				self.cells[y*self.width+x].flux = 0.0
		
		# compute scalar flux in each cell
		for y in range(self.width):
			for x in range(self.width):
				cell = self.cells[y*self.width+x]

				for angle in range(self.num_angles/4):
					cell.flux += self.quad['weight'][angle] * (cell.ang_flux[angle] + cell.ang_flux[angle + self.num_angles/4] + cell.ang_flux[angle + self.num_angles/2] + cell.ang_flux[angle + 3*self.num_angles/4])


		# compute average scalar flux
		avg_flux = 0.0	
		for y in range(self.width):
			for x in range(self.width):
				avg_flux += self.cells[y*self.width+x].flux

		avg_flux = avg_flux / (self.width**2)


		# normalize flux
		for y in range(self.width):
			for x in range(self.width):
				self.cells[y*self.width+x].flux = self.cells[y*self.width+x].flux / avg_flux		





def main():

	print 'parsing command line input...'
	
	# parse commandl line options
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:s:i:t:", ["order", "size", "num_iter", "tolerance"])
	except getopt.GetoptError, err:
		print str(err)
		usage()
		sys.exit(2)

	sn_order = 4
	cell_size = .01
	tol = .001
	num_iter = 50

	for o, a in opts:
		if o in ("-o", "--order"):
			sn_order = int(a)
		elif o in ("-s", "--size"):
			cell_size = float(a)
		elif o in ("-i", "--iterations"):
			num_iter = int(a)
		elif o in ("-t", "--tolerance"):
			tol = float(a)
		else:
			assert False, "unhandled option"


	mesh = Mesh(order=sn_order, mesh_size=cell_size, iterations=num_iter, tolerance=tol)
	
	fuel = Material('fuel', 100.0)
	moderator = Material('moderator', 0.25)

	mesh.setFuel(fuel)
	mesh.setModerator(moderator)

	mesh.makeMeshMaterials()
	
	pttr.plotMaterial(mesh)

	mesh.initializeSource()

	start = time.time()

	mesh.solveSn()

	stop = time.time()

	print 'Ran Sn solver with ' + str(sn_order*(sn_order+2)/2) + ' angles in ' + str(stop-start) + ' seconds'



if __name__ == '__main__':

	main()












