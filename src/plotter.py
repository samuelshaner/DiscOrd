
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
from math import *
import matplotlib.pyplot as plt
import numpy as np

def plotMaterial(mesh, spacing, plot_cells):

    bit_size = 500.0

    # create image
    bit_length = round(bit_size / mesh.width)
    size = int(bit_length * mesh.width)
    img = Image.new('RGB', (size,size), 'white')
    draw = ImageDraw.Draw(img)

    # draw cell interior
    for y in range(mesh.width):
        for x in range(mesh.width):

            # fuel red; moderator blue
            if mesh.cells[y*mesh.width+x].material.mat_type == 'fuel':
                draw.rectangle([x*bit_length, size - y*bit_length - bit_length, (x+1)*bit_length, size - y*bit_length], (255,0,0))
            else:
                draw.rectangle([x*bit_length, size - y*bit_length - bit_length, (x+1)*bit_length, size - y*bit_length], (0,0,255))

            if mesh.cells[y*mesh.width+x].id in plot_cells:
                draw.rectangle([x*bit_length, size - y*bit_length - bit_length, (x+1)*bit_length, size - y*bit_length], (255,255,255))

    # draw horizontal grid lines
    for y in range(1,mesh.width):
        draw.line((0, y*bit_length, size,y*bit_length), fill=(0,0,0), width=1)

    # draw vertical grid lines
    for x in range(1,mesh.width):
        draw.line((x*bit_length, 0, x*bit_length, size), fill=(0,0,0), width=1)

    # save image
    img.save('material_' + str(spacing)[2:] + '.png')

def plotScalarFlux(mesh, order, spacing, iteration):

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
            draw.rectangle([x*bit_length, size - y*bit_length - bit_length, (x+1)*bit_length, size - y*bit_length], (red,green,blue))

    # save image
    img.save('flux_' + str(spacing)[2:] + '_' + str(int(floor(order/10))) + str(order % 10) + '_' + str(int(floor(iteration/10))) + str(iteration % 10) + '.png')


def plotAngularFlux(cell, quad):

    # create image
    img = Image.new('RGB', (500,500), 'white')
    draw = ImageDraw.Draw(img)

    # get mus, eta's, and index for smallest polar angle (xi)
    xi_min = min(quad['xi'])
    min_xi_indices = []
    for i,a in enumerate(quad['xi']):
        if abs(a - xi_min) < .0001:
            min_xi_indices.append(i)

    # get the angular fluxes for the selected angles and make theta array
    ang_fluxes = []
    theta = []

    for i in min_xi_indices:
        ang_fluxes.append(cell.ang_flux[i])
        theta.append(acos(quad['mu'][i]))

    for i in min_xi_indices:
        ang_fluxes.append(cell.ang_flux[quad['num_angles_per_octant'] + i])
        theta.append(pi - acos(quad['mu'][i]))

    for i in min_xi_indices:
        ang_fluxes.append(cell.ang_flux[2*quad['num_angles_per_octant'] + i])
        theta.append(pi + acos(quad['mu'][i]))

    for i in min_xi_indices:
        ang_fluxes.append(cell.ang_flux[3*quad['num_angles_per_octant'] + i])
        theta.append(2*pi - acos(quad['mu'][i]))


    # rarrange theta array from min to max
    theta_min = min(theta)
    theta_sort = []
    ang_fluxes_sort = []
    for j in range(len(theta)):
        for i,a in enumerate(theta):
            if a == theta_min:
                theta_sort.append(a)
                ang_fluxes_sort.append(ang_fluxes[i])
                theta.remove(theta_min)
                ang_fluxes.remove(ang_fluxes[i])
                if len(theta) != 0:
                    theta_min = min(theta)

    # append 0 and 2 pi values
    ang_flux_max = max(ang_fluxes_sort)
    ang_flux_0 = ang_fluxes_sort[0]
    theta_0 = - theta_sort[0]
    ang_fluxes_sort.insert(0, ang_flux_0)
    theta_sort.insert(0, theta_0)
    ang_flux_2pi = ang_fluxes_sort[-1]
    theta_2pi = 2*pi + (2*pi - theta_sort[-1])
    ang_fluxes_sort.append(ang_flux_2pi)
    theta_sort.append(theta_2pi)

    # draw slices for each theta
    for i, angle in enumerate(theta_sort[1:-1]):
        i = i+1

        radius = int(ang_fluxes_sort[i]/ang_flux_max*200)
        ang_start = int((theta_sort[i]+theta_sort[i-1])/2*180/pi)
        ang_end   = int((theta_sort[i+1]+theta_sort[i])/2*180/pi)

        #               print 'plotting slice radius ' + str(radius) + ' theta_min ' + str(ang_start) + ' theta_max ' + str(ang_end)

        ang_start = 360 - ang_start
        ang_end = 360 - ang_end

        draw.pieslice((250 - radius, 250 - radius, 250 + radius, 250 + radius), ang_end, ang_start, fill=(0,0,0), outline=(255,255,255))

    fig = plt.figure()
    plt.plot(theta_sort, ang_fluxes_sort)

    fig.savefig('ang_flux.png')
    img.save('ang_flux_circle_' + str(cell.id) + '.png')
