
import Image
import ImageDraw
import ImageFont
from math import *

def plotMaterial(mesh):

	# create image
	bit_length = round(500.0 / mesh.width)
	size = int(bit_length * mesh.width)
	img = Image.new('RGB', (size,size), 'white')
	draw = ImageDraw.Draw(img)
	

	# draw cell interior
	for y in range(mesh.width):
		for x in range(mesh.width):

			# fuel red; moderator blue
			if mesh.cells[y*mesh.width+x].material.mat_type == 'fuel':
				draw.rectangle([x*bit_length, y*bit_length, (x+1)*bit_length, (y+1)*bit_length], (255,0,0))
			else:
				draw.rectangle([x*bit_length, y*bit_length, (x+1)*bit_length, (y+1)*bit_length], (0,0,255))


	# draw horizontal grid lines
	for y in range(1,mesh.width):
		draw.line((0, y*bit_length, size,y*bit_length), fill=(0,0,0))

	# draw vertical grid lines
	for x in range(1,mesh.width):
		draw.line((x*bit_length, 0, x*bit_length, size), fill=(0,0,0))

	# save image
	img.save('material.png')
    

def plotScalarFlux(mesh, order, iteration):

	# create image
	bit_length = round(500.0 / mesh.width)
	size = int(bit_length * mesh.width)
	img = Image.new('RGB', (size,size), 'white')
	draw = ImageDraw.Draw(img)

	# find max and min flux
	max_flux = mesh.cells[0].flux
	min_flux = mesh.cells[0].flux
	for y in range(mesh.width):
		for x in range(mesh.width):
			max_flux = max(mesh.cells[y*mesh.width+x].flux, max_flux)
			min_flux = min(mesh.cells[y*mesh.width+x].flux, min_flux)

	flux_range = max_flux - min_flux

	# draw flux map
	for y in range(mesh.width):
		for x in range(mesh.width):
			cell = mesh.cells[y*mesh.width+x]


			# get color
			if ((cell.flux-min_flux) / flux_range <= 1.0/3.0):
				red = 0.0
				green = 3.0 * (cell.flux-min_flux) / flux_range
				blue = 1.0
			elif ((cell.flux-min_flux) / flux_range <= 2.0/3.0):
				red = 3.0 * (cell.flux-min_flux) / flux_range - 1.0
				green = 1.0
				blue = -3.0 * (cell.flux-min_flux) / flux_range + 2.0
			else:
				red = 1.0
				green = -3.0 * (cell.flux-min_flux) / flux_range + 3.0
				blue = 0.0

			# convert color to RGB triplet
			red = int(255*red)
			green = int(255*green)
			blue = int(255*blue)

			# draw pin and pin power
			draw.rectangle([x*bit_length, y*bit_length, (x+1)*bit_length, (y+1)*bit_length], (red,green,blue))

	# save image
	img.save('flux_' + str(int(floor(order/10))) + str(order % 10) + '_' + str(int(floor(iteration/10))) + str(iteration % 10) + '.png')
touc









