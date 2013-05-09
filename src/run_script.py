
import os

orders = range(2,26,2)
spacing = .005

for i in orders:

	os.system('python discrete_ordinates.py -o ' + str(i) + ' -s ' + str(spacing))



