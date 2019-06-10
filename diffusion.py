import pandas as pd
from decimal import Decimal

data = pd.read_csv('OUT')
timestep = data[u'timestep']
vacf = data[u' vacf']

trapAreas = 0
for x in range(len(vacf)-1):
    trapAreas += (vacf[x] + vacf[x+1])

area = 10*trapAreas/2

diffCoeff = area/3

diffCoeff = diffCoeff * 1e-7

print('Coefficient of diffusion: ', '%.2E'%Decimal(str(diffCoeff)))

