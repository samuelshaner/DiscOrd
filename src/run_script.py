
import os

orders = [2,4,8,16,24]
spacing = [.005]

for j in spacing:
	for i in orders:

		if i == 24:
			os.system('python discrete_ordinates.py -o ' + str(i) + ' -s ' + str(j) + ' -f True')
		else:
			os.system('python discrete_ordinates.py -o ' + str(i) + ' -s ' + str(j))

		


