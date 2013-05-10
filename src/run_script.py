
import os

orders = [16,24]
spacing = [.005,.01]

for j in spacing:
	for i in orders:

		os.system('python discrete_ordinates.py -o ' + str(i) + ' -s ' + str(j))



