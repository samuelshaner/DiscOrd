import quadrature
import numpy

orders = numpy.arange(2, 24, 2)
quad = quadrature.LevelSymmetricQuadrature()

for order in orders:
    
    quadrature = quad.getQuadrature(order)
    
    print '-------------------------------------------------------------'
    print '  mu\t\t  eta\t\t  xi\t\t  weight'
    print '-------------------------------------------------------------'
    
    for angle in range(quadrature['num_angles_per_octant']):
        
        print '%f\t%f\t%f\t%f' % (quadrature['mu'][angle], \
                                  quadrature['eta'][angle], \
                                  quadrature['xi'][angle], \
                                  quadrature['weight'][angle])