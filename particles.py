filenames = ["%04d"%(a,) for a in range(1,6913,1)]

filedata = {filename: open(filename, 'w') for filename in filenames}

with open("HISTORY_CLEAN_6912_l", "r") as f:
	for line in f:
		filedata["%04d"%(int(line.split()[1]),)].write(line)
f.close()
for file in filedata.values():
	file.close()
