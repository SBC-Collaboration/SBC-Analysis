#Author: Trent Cwiok
import matplotlib.pyplot as plt

#execfile('ReadText.py')
#execfile('ReadBinary.py')

#THIS FUNCTION CAN BE USED TO IDENTIFY THE TIME OF BUBBLE FORMATION
#############################################################################################

def PlotAcoustics(ev):
	time_step = ev['fastDAQ']['caldata']['dt'][0]
	origin = ev['fastDAQ']['caldata']['pretrigger_samples']*time_step
	PMTtrig_time_scale = [time_step*i-origin for i in range(len(ev['fastDAQ']['PMTtrig']))]

	plt.plot(PMTtrig_time_scale, ev['fastDAQ']['PMTtrig'], 'm-', PMTtrig_time_scale, ev['fastDAQ']['Piezo1'], 'b-')
	plt.show()
	return


#THIS FUNCTION CAN BE USED TO MANUALLY MATCH PMT TRACES TO PMT HITS
#############################################################################################

def MatchTraces():
	print "Loading fastDAQ..."
	fastDAQ = ReadBlock('fastDAQ_0.bin')
	print "Finished."

	print "Loading fastDAQ_cal..."
	fastDAQ_cal = ReadFile('fastDAQ_0_cal.txt')
	print "Finished."

	print "Loading PMT..."
	PMT = ReadBlock('PMTtraces.bin')
	print "Finished."

	time_step = fastDAQ_cal['dt'][0]
	origin = fastDAQ_cal['pretrigger_samples']*time_step
	PMTtrig_time_scale = [time_step*i-origin for i in range(len(fastDAQ['PMTtrig']))]

	PMTtrig_times = [PMT['t0_sec'][i][0]+PMT['t0_frac'][i][0] for i in range(len(PMT['t0_sec']))]
	for ii in range(len(PMTtrig_times)):
		PMTtrig_times[ii] -= PMTtrig_times[-1]
	amplitude = [100 for i in range(len(PMTtrig_times))]

	#Offset must be manually observed
	offset = -.13175

	counter = 0
	index = len(PMTtrig_times)-1
	while index > 1:
		for ii in range(len(PMTtrig_times)):
			PMTtrig_times[ii] += offset

		print "This is check number", counter+1, ".  Continue?"
		x = raw_input("Continue?\n")
		if x == ('no' or 'break' or 'n'):
			break
		plt.plot(PMTtrig_time_scale, fastDAQ['PMTtrig'], 'm-', PMTtrig_times, amplitude, 'bo')
		plt.axis([-.2,.2, 5, 150])
		plt.show()

		offset = PMTtrig_times[index] - PMTtrig_times[index-1]

		counter += 1
		index -= 1
	return
