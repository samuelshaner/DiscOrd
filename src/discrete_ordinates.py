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

    def __init__(self, mesh_size, order, tolerance, pitch=1.26, fuel_diameter=0.70, geometry='square'):

        # initialize variables and create lists for storing values
        self.width         = int(pitch//mesh_size)
        self.quad          = quadrature.LevelSymmetricQuadrature().getQuadrature(order)
        self.fuel_diameter = fuel_diameter
        self.mesh_size     = mesh_size
        self.order         = order
        self.tol           = tolerance
        self.num_angles    = self.quad['num_angles']
        self.pitch         = pitch
        self.cells         = np.empty(self.width**2, dtype=object)
        self.ang_flux      = np.zeros((4,self.width,self.num_angles//2))
        self.geometry      = geometry

        # create the mesh cells
        for y in range(self.width):
            for x in range(self.width):
                self.cells[y*self.width+x] = Cell(self.num_angles, y*self.width+x)


    # set pointer to fuel material
    def setFuel(self, fuel):

        self.fuel = fuel

    # set pointer to moderator material
    def setModerator(self, moderator):

        self.moderator = moderator

    # give each cell a pointer to a material
    def makeMeshMaterials(self):

        cw = self.width

        width_to_fuel = cw / 2 - int(self.fuel_diameter / self.mesh_size / 2)

        if self.geometry == 'square':
            for y in range(cw):
                for x in range(cw):
                    if x < width_to_fuel or x >= cw - width_to_fuel:
                        self.cells[y*cw+x].setMaterial(self.moderator)
                    else:
                        if y < width_to_fuel or y >= self.width - width_to_fuel:
                            self.cells[y*cw+x].setMaterial(self.moderator)
                        else:
                            self.cells[y*cw+x].setMaterial(self.fuel)

        elif self.geometry == 'circle':
            for y in range(cw):
                for x in range(cw):
                    radius = sqrt((self.pitch/2.0-(y+0.5)*self.mesh_size)**2 + (self.pitch/2.0-(x+0.5)*self.mesh_size)**2)

                    if radius <= self.fuel_diameter/2.0:
                        self.cells[y*cw+x].setMaterial(self.fuel)
                    else:
                        self.cells[y*cw+x].setMaterial(self.moderator)

    # solve the Sn problem
    def solveSn(self, update, num_iter):

        # create shortened variables for num_angles and width and initialize eps
        cw = self.width
        na = self.num_angles
        eps = 1.0

        # loop over iterations
        for iteration in range(num_iter):

            print('Sn iteration ' + str(iteration) + ' eps ' + str(eps))

            ###################################################################
            # sweep in the positive mu and eta direction (bottom left corner)
            ###################################################################
            for y in range(cw):
                for x in range(cw):

                    cell = self.cells[y*cw+x]
                    source = cell.material.source

                    for angle in range(na//4):

                        # compute cell centered flux
                        cell.ang_flux[angle] = (source + 2 * self.quad['mu'][angle]/self.mesh_size * self.ang_flux[0, y, angle] +
                                                2 * self.quad['eta'][angle]/self.mesh_size * self.ang_flux[1, x, angle]) / \
                                                (cell.material.sigma_t + 2 * self.quad['mu'][angle]/self.mesh_size + 2 * self.quad['eta'][angle]/self.mesh_size)

                        # sweep from left to right
                        self.ang_flux[0, y, angle] = 2 * cell.ang_flux[angle] - self.ang_flux[0,y, angle]

                        # sweep from bottom to top
                        self.ang_flux[1, x, angle] = 2 * cell.ang_flux[angle] - self.ang_flux[1, x, angle]

            # update the boundary angular fluxes
            if update:
                # reflect on right boundary
                for y in range(cw):
                    for angle in range(na//4):
                        self.ang_flux[2, y, angle] = self.ang_flux[0, y, angle]

                # reflect flux on top boundary
                for x in range(cw):
                    for angle in range(na//4):
                        self.ang_flux[3, x, angle] = self.ang_flux[1, x, angle]


            ###################################################################
            # sweep in the negative mu and eta direction (top right corner)
            ###################################################################
            for y in reversed(range(cw)):
                for x in reversed(range(cw)):

                    cell = self.cells[y*cw+x]
                    source = cell.material.source

                    for angle in range(na//4):

                        # compute cell centered flux
                        cell.ang_flux[na//2 + angle] = (source + 2 * self.quad['mu'][angle]/self.mesh_size * self.ang_flux[2, y, angle + na//4] +
                                                        2 * self.quad['eta'][angle]/self.mesh_size * self.ang_flux[3, x, angle + na//4]) / \
                                                        (cell.material.sigma_t + 2 * self.quad['mu'][angle]/self.mesh_size + 2 * self.quad['eta'][angle]/self.mesh_size)

                        # sweep from right to left
                        self.ang_flux[2, y, angle + na//4] = 2 * cell.ang_flux[na//2 + angle] - self.ang_flux[2, y, angle + na//4]

                        # sweep from top to bottom
                        self.ang_flux[3, x, angle + na//4] = 2 * cell.ang_flux[na//2 + angle] - self.ang_flux[3, x, angle + na//4]

            # update the boundary angular fluxes
            if update:
                # reflect on left boundary
                for y in range(cw):
                    for angle in range(na//4):
                        self.ang_flux[0, y, angle + na//4] = self.ang_flux[2, y, angle + na//4]

                # reflect flux on bottom boundary
                for x in range(cw):
                    for angle in range(na//4):
                        self.ang_flux[1, x, angle + na//4] = self.ang_flux[3, x, angle + na//4]



            # sweep in the negative mu and positive eta direction (bottom right corner)
            for y in range(cw):
                for x in reversed(range(cw)):

                    cell = self.cells[y*cw+x]
                    source = cell.material.source

                    for angle in range(na//4):

                        # compute cell centered flux
                        cell.ang_flux[na//4 + angle] = (source + 2 * self.quad['mu'][angle]/self.mesh_size * self.ang_flux[2, y, angle] +
                                                        2 * self.quad['eta'][angle]/self.mesh_size * self.ang_flux[1, x, angle + na//4]) / \
                                                        (cell.material.sigma_t + 2 * self.quad['mu'][angle]/self.mesh_size + 2 * self.quad['eta'][angle]/self.mesh_size)

                        # sweep from right to left
                        self.ang_flux[2, y, angle] = 2 * cell.ang_flux[na//4 + angle] - self.ang_flux[2, y, angle]

                        # sweep from bottom to top
                        self.ang_flux[1, x, angle + na//4] = 2 * cell.ang_flux[na//4 + angle] - self.ang_flux[1, x, angle + na//4]

            # update the boundary angular fluxes
            if update:
                # reflect on left boundary
                for y in range(cw):
                    for angle in range(na//4):
                        self.ang_flux[0, y, angle] = self.ang_flux[2, y, angle]

                # reflect flux on top boundary
                for x in range(cw):
                    for angle in range(na//4):
                        self.ang_flux[3, x, angle + na//4] = self.ang_flux[1, x, angle + na//4]


            # sweep in the positive mu and negative eta direction (top left corner)
            for y in reversed(range(cw)):
                for x in range(cw):

                    cell = self.cells[y*cw+x]
                    source = cell.material.source

                    for angle in range(na//4):

                        # compute cell centered flux
                        cell.ang_flux[3*na//4 + angle] = (source + 2 * self.quad['mu'][angle]/self.mesh_size * self.ang_flux[0, y, angle + na//4] +
                                                          2 * self.quad['eta'][angle]/self.mesh_size * self.ang_flux[3, x, angle]) / \
                                                          (cell.material.sigma_t + 2 * self.quad['mu'][angle]/self.mesh_size + 2 * self.quad['eta'][angle]/self.mesh_size)

                        # sweep from left to right
                        self.ang_flux[0, y, angle + na//4] = 2 * cell.ang_flux[3*na//4+ angle] - self.ang_flux[0, y, angle + na//4]

                        # sweep from bottom to top
                        self.ang_flux[3, x, angle] = 2 * cell.ang_flux[3*na//4+ angle] - self.ang_flux[3, x, angle]

            # update the boundary angular fluxes
            if update:
                # reflect on right boundary
                for y in range(cw):
                    for angle in range(na//4):
                        self.ang_flux[2, y, angle + na//4] = self.ang_flux[0, y, angle + na//4]

                # reflect flux on bottom boundary
                for x in range(cw):
                    for angle in range(na//4):
                        self.ang_flux[1, x, angle] = self.ang_flux[3, x, angle]

            # if vacuum case, plot scalar flux and compute fuel rxn rate
            if update == False:
                self.RR_isolated = self.computeRRFuel()

                # plot the scalar flux for vacuum boundary case
                self.computeScalarFlux()
                pttr.plotScalarFlux(self, self.order, self.mesh_size, iteration+100)

                # zero out angular flux
                for i in range(4):
                    for x in range(cw):
                        for angle in range(na//2):
                            self.ang_flux[i, x, angle] = 0.0

            # plot the scalar flux for reflective boundary case and compute eps
            if update:
                self.computeScalarFlux()
                pttr.plotScalarFlux(self, self.order, self.mesh_size, iteration)
                eps = self.computeL2Norm()

            # check for convergence; if convgerged comput dancoff factor and flux ratio
            if (update and eps < self.tol):
                self.RR_lattice = self.computeRRFuel()
                self.RR_total = self.computeRRTotal()
                self.flux_ratio = self.computeFluxRatio()
                self.dancoff = 1 - (1.0 - self.RR_lattice / self.RR_total) / (1.0 - self.RR_isolated / self.RR_total)
                self.dancoff2 = 1 - (1.0 / self.fuel.sigma_t)
                print('EPS converged ' + str(eps))
                print('RR isolated ' + str(self.RR_isolated))
                print('RR lattice ' + str(self.RR_lattice))
                print('flux ratio ' + str(self.flux_ratio))
                print('Dancoff factor ' + str(self.dancoff))
                break

    # compute the scalar fuel to coolant ratio
    def computeFluxRatio(self):

        # zero flux
        flux_fuel = 0.0
        flux_coolant = 0.0
        fuel_cells = 0
        coolant_cells = 0

        # compute scalar flux in each cell
        for y in range(self.width):
            for x in range(self.width):
                cell = self.cells[y*self.width+x]

                if cell.material.mat_type == 'fuel':
                    flux_fuel += cell.flux
                    fuel_cells += 1
                else:
                    flux_coolant += cell.flux
                    coolant_cells += 1

        ratio = (flux_fuel/fuel_cells) / (flux_coolant/coolant_cells)

        return ratio


    # compute the rxn rate in the fuel
    def computeRRFuel(self):

        # set old flux
        for y in range(self.width):
            for x in range(self.width):
                self.cells[y*self.width+x].flux = 0.0

        # compute scalar flux in each cell
        for y in range(self.width):
            for x in range(self.width):
                cell = self.cells[y*self.width+x]

                for angle in range(self.num_angles//4):
                    cell.flux += self.quad['weight'][angle] * (cell.ang_flux[angle] + cell.ang_flux[angle + self.num_angles//4] + cell.ang_flux[angle + self.num_angles//2] + cell.ang_flux[angle + 3*self.num_angles//4])

        # compute average scalar flux
        RR_fuel = 0.0
        for y in range(self.width):
            for x in range(self.width):

                if self.cells[y*self.width+x].material.mat_type == 'fuel':
                    RR_fuel += self.cells[y*self.width+x].flux * self.cells[y*self.width+x].material.sigma_t * self.mesh_size**2

        return RR_fuel

    # compute the total rxn rate
    def computeRRTotal(self):

        # set old flux
        for y in range(self.width):
            for x in range(self.width):
                self.cells[y*self.width+x].flux = 0.0

        # compute scalar flux in each cell
        for y in range(self.width):
            for x in range(self.width):
                cell = self.cells[y*self.width+x]

                for angle in range(self.num_angles//4):
                    cell.flux += self.quad['weight'][angle] * (cell.ang_flux[angle] + cell.ang_flux[angle + self.num_angles//4] + cell.ang_flux[angle + self.num_angles//2] + cell.ang_flux[angle + 3*self.num_angles//4])

        # compute average scalar flux
        RR_total = 0.0
        for y in range(self.width):
            for x in range(self.width):
                RR_total += self.cells[y*self.width+x].flux * self.cells[y*self.width+x].material.sigma_t * self.mesh_size**2

        return RR_total

    # compute the L2 Norm of the scalar flux between iterations
    def computeL2Norm(self):

        eps = 0.0

        for y in range(self.width):
            for x in range(self.width):
                eps += (self.cells[y*self.width+x].flux - self.cells[y*self.width+x].old_flux)**2

        sqrt(eps)

        return eps


    # compute and normalize the scalar flux
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

                for angle in range(self.num_angles//4):
                    cell.flux += self.quad['weight'][angle] * (cell.ang_flux[angle] + cell.ang_flux[angle + self.num_angles//4] + cell.ang_flux[angle + self.num_angles//2] + cell.ang_flux[angle + 3*self.num_angles//4])


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
