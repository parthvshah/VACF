import sys
import time
import datetime
import subprocess as sp

jobs = list()

for i in range(1, 10, 1):
	try:
		jobs.append(sys.argv[i])
	except IndexError:
		pass

for job in jobs:
	while True:
		time.sleep(10)
		try:
			output = sp.check_output(['qsub', job], stderr=sp.STDOUT)
		except sp.CalledProcessError as e:
			output = e.output
	
		if('submit error' in  output):
			print '[Error] Job', job, 'failed submission at', datetime.datetime.now()
			continue
		else:
			print 'Job', job, 'submitted at', datetime.datetime.now(), 'with JOBID', output
			break
