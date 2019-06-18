import sys
import csv
from decimal import Decimal

timestep = list()
vacf = list()

fileName = sys.argv[1]

with open(fileName) as csv_file:
    csvReader = csv.reader(csv_file, delimiter=',')
    line = 0
    for row in csvReader:
        if(line != 0):
            timestep.append(int(row[0]))
            vacf.append(float(row[1]))
            line += 1
        else:
            line += 1

trapAreas = 0
for x in range(len(vacf)-1):
    trapAreas += (vacf[x] + vacf[x+1])

area = 10*trapAreas/2

diffCoeff = area/3

diffCoeff = diffCoeff * 1e-7

print('Coefficient of diffusion: ', '%.2E'%Decimal(str(diffCoeff)))

