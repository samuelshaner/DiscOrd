from discrete_ordinates import Mesh
from material import Material
import plotter
import quadrature
import time
import numpy as np

orders = [2,4,8,16,24]
spacings = [.005]
tol = .0001
num_iter = 50
geom = 'square'
plot_flux = False

for spacing in spacings:
    for order in orders:

        # create mesh
        mesh = Mesh(order=order, mesh_size=spacing, tolerance=tol, geometry=geom)

        print('CASE - mesh_size {} order {} geometry {}'.format(spacing, order, geom))

        # create fuel and moderator materials
        fuel = Material('fuel', 100.0, 1.0/(4.0*np.pi))
        moderator = Material('moderator', 0.25, 0.0)

        cell_width = int(1.26 / spacing)
        fuel_width = int(0.70 / spacing)
        # 1, 2, 3, 4, 5
        plot_cells = [int((cell_width/2.0)*cell_width+cell_width/2.0),
                      int((cell_width/2.0 + fuel_width/2-1)*cell_width+cell_width/2.0) + fuel_width/2-1,
                      int(cell_width*cell_width - 1),
                      int((cell_width/2.0)*cell_width+cell_width/2.0) + fuel_width/2 -1,
                      int((cell_width/2.0 + 1)*cell_width - 1)]

        # assign fuel and moderator materials to mesh
        mesh.setFuel(fuel)
        mesh.setModerator(moderator)
        mesh.makeMeshMaterials()

        # plot the materials in the mesh
        plotter.plotMaterial(mesh, spacing, plot_cells)

        # solve the vacuum boundary Sn problem
        mesh.solveSn(False, 1)

        start = time.time()

        # solve the reflective boundary Sn problem
        mesh.solveSn(True, num_iter)

        stop = time.time()

        print('Ran Sn solver with {} angles in {:.2f} seconds'.format(order*(order+2)//2, stop-start))

        if plot_flux:
            quad = quadrature.LevelSymmetricQuadrature().getQuadrature(sn_order)
            for i in plot_cells:
                print('plotting angular flux for cell ' + str(i))
                pttr.plotAngularFlux(mesh.cells[i], quad)

        print('----------------------------------------------------------------------')
