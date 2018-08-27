#Author: Trent Cwiok
'''
To use this sample script:
1) Open the Command Prompt/Terminal
2) Have this file and 'slowDAQ_0' in the same directory
3) Navigate to this directory in the Command Prompt/Terminal
3) Make sure matplotlib is installed
	a) This is most easily done through pip
	b) Type "pip install pip"
	c) Then "pip install matplotlib"
4) Next, open a python command line by simply typing "python"
5) Finally, run the file by calling "execfile('Python Sample.py')"
'''

import numpy as np
import matplotlib.pyplot as plt

def ReadInFromFile(filename, headersize):
	instrument_dict = {}
	with open(filename, 'r') as read_in:
		header = [next(read_in).strip() for x in xrange(headersize)]
		instrument_list = header[1].split()
		for var in range(len(instrument_list)):
			instrument_dict[instrument_list[var]] = []

		while True:
			try:
				line = next(read_in)
				data = line.split()

				for i in range(len(instrument_list)):
					instrument_dict[instrument_list[i]].append(data[i])
					
			except StopIteration:
				break
			except Exception:
				print "Unknown Exception"
	for key in instrument_dict:
		instrument_dict[key] = np.asarray(instrument_dict[key])
	return instrument_dict

d = ReadInFromFile('slowDAQ_0.txt', 6)

plt.plot(d['elapsed_time'], d['PT4'])
plt.show()