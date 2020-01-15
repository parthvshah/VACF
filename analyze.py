f = open('./HISTORY_atoms/HISTORY_CLEAN_108','r')

lines = f.readlines()

x = []
y = []
z = []

for line in lines:
    splitValues = line.split()
    x.append(float(splitValues[2]))
    y.append(float(splitValues[3]))
    z.append(float(splitValues[4]))

import statistics as stats
print("x - Max: ", max(x), "\tMin: ", min(x), "\tMean: ", stats.mean(x), "\tSD: ", stats.stdev(x))
print("y - Max: ", max(y), "\tMin: ", min(y), "\tMean: ", stats.mean(y), "\tSD: ", stats.stdev(y))
print("z - Max: ", max(z), "\tMin: ", min(z), "\tMean: ", stats.mean(z), "\tSD: ", stats.stdev(z))

