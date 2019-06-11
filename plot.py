import matplotlib.pyplot as plt
import csv
import sys

timestep = list()
vacf = list()

lower = int(sys.argv[1])
upper = int(sys.argv[2])

with open('OUT') as csv_file:
    csvReader = csv.reader(csv_file, delimiter=',')
    line = 0
    for row in csvReader:
        if(line != 0):
            timestep.append(int(row[0]))
            vacf.append(float(row[1]))
            line += 1
        else:
            line+=1

scalar = vacf[0]
for x in range(len(vacf)):
    vacf[x] /= scalar
    timestep[x] /= 100

if(lower==0 and upper==0):
    plt.plot(timestep, vacf, color="black")
else:
    plt.plot(timestep[lower:upper], vacf[lower:upper], color="black")

plt.axhline(y=0, color="blue")
plt.title("VACF for atoms of Argon (liq)")
plt.xlabel("Time in 1e-12s")
plt.ylabel("Velocity auto-correlation")
plt.savefig('vacf.png')