#Author: Trent Cwiok
# 
# The purpose of this program is to, given a file name containing
# data collected from the instruments attached to the bubble chamber,
# create a data structure capable of storing the instruments and
# their readouts. The information for any given instrument should
# then be easily accessible and printable. Additional functionality
# for plotting any given instrument readout versus time is also
# included.

import pdb
import numpy as np
import matplotlib.pyplot as plt
#import SBCcode as sbc


np.set_printoptions(threshold = np.nan)


def ReadFile(filename):
	'''
	This function takes in the file's name to be read, as well as
	the size of the header containing metadata, and returns a
	dictionary where the individual instrument names are the keys,
	and the associated values are arrays of the value recorded by
	the instrument at each timestep. Returns None when errors.
	'''

	instrument_dict = {}

	#Get the number of lines in the file for appropriate array construction
	with open(filename, 'r') as f:
		for i, l in enumerate(f):
			pass
		#Subtract five to account for the header size
		num_lines = i - 5

	with open(filename, 'r') as read_in:
		#Read in the file header
		header = [next(read_in).strip() for x in xrange(6)]
		instrument_list = header[1].split()

		#Create keys in the dictionary for each variable and
		#initialize an empty array to store the data
		for var in instrument_list:
			shape = num_lines
			data_type = 'double'
			if '(' in var:
				dimensions = var[var.find("(")+1:var.find(")")]
				shape = str(num_lines)+','+dimensions
				shape = tuple(map(int, shape.split(',')))
			if var == 'run':
				data_type = 'string'
			instrument_dict[var] = np.zeros(shape, dtype = data_type, order = 'C')

		#This loop runs until the file end is reached, taking
		#in the data in each line, parsing it, and storing it
		#to the appropriate key in the dictionary
		index = 0
		while index < num_lines:
			line = next(read_in)
			data = line.split()

			if len(data) == len(instrument_list):
				for i in range(len(instrument_list)):
					instrument_dict[instrument_list[i]][index] = data[i]

			else:
				front = 0
				for var in instrument_list:
					data_type = 'double'
					if var == 'run':
						data_type = 'string'
					if '(' not in var:
						pdb.set_trace()
						instrument_dict[var][index] = data[front]
						front += 1

					#NOTE: This part of the loop is meant to handle
					#variables with multi-dimensional entries, but is
					#untested due to the absences of any test file
					else:
						shape = var[var.find('(')+1:var.find(')')]
						shape = tuple(map(int, shape.split(',')))
						iterations = 1
						for ele in iterations:
							iterations *= ele

						temp_buffer = np.zeros(counter, dtype = data_type, order = 'C')
						counter = 0
						while counter < iterations:
							temp_buffer[counter] = data[front]
							front += 1
							counter += 1
						temp_buffer = np.reshape(temp_buffer, shape, order = 'C')
						instrument_dict[var][index] = temp_buffer

			index += 1

	return instrument_dict


def DataTrim(dictionary, instrument):
	'''
	Given a dictionary constaining instrument data, it uses 
	TriggerLatch to trim away unwanted data points where the trigger is not latched
	'''
	tick=0
	
	
	instrument = str(instrument)
	for i in range (len(dictionary['TriggerLatch'])):
		if dictionary['TriggerLatch'][i] == 0:
			tick=tick+1
	trim_PT=np.zeros(tick)
	trim_time=np.zeros(tick)
	track =0
	for value in range(len(dictionary[instrument])):
		if dictionary['TriggerLatch'][value] == 0:
			trim_PT[track]=dictionary[instrument][value]
			trim_time[track]=dictionary['elapsed_time'][value]
			track=track+1
		
	return trim_PT
#
#DataTrim(d,'PT4')
 
def TrimAll(dictionary):
	
	d={}
	for i in range(1,10):
		index= str(i)
		d.update({'trimPT' +str(i): DataTrim(dictionary, 'PT'+index)})
	return d


def ShowIndex(dictionary):
	for key in dictionary:
		print key

step = []
for i in range(0,250):
	step.append(i)

def Pbin(dictionary, instrument, edge ):
	'''
	sort pressures into bins of 1 psi
	'''
	
	Bin=np.histogram(dictionary[instrument], bins=edge)
	BinTime=np.zeros(len(Bin[0]))
	for i in range(len(Bin[0])):
		x=Bin[0][i]*0.005
		BinTime[i]=x
		dictionary.update({'Bin'+instrument:Bin[1], 'Count'+instrument:Bin[0], 'BinTime'+instrument:BinTime})
	'''
	print Bin
	plt.hist(dictionary[instrument], bins=step)
	plt.show()
	'''
	return

def BinAll(dictionary, edge):
	for i in range(1,10):
		index=str(i)
		Pbin(dictionary,'trimPT'+index, edge) 
	return


def tGood(dictionary, instrument='PT4'):
	Pgood=0
	tg=0
	
	for i in range(1,len(dictionary['PT4'])):
		if abs(dictionary['PT4'][i]-dictionary['PressureSetpoint'][i]) <= 0.3:
			Pgood=Pgood +1
			tg = tg+((dictionary['elapsed_time'][i]-dictionary['elapsed_time'][i-1]))
	return tg

def tEvent(dictionary):
	Tevent=0
	dt = np.diff(dictionary['elapsed_time'][0:2])
	Tevent = np.sum(dictionary['TriggerLatch'] == 0) * dt

	latch_lowtohigh = np.nonzero(np.diff(dictionary['TriggerLatch'])==1)
	time_of_compression = dictionary['elapsed_time'][latch_lowtohigh[-1]]

	return Tevent

def PumpActiveTime(dictionary):
	tPumpPre=0
	tPumpPost=0
	for i in range(1,len(dictionary['PUMP'])/2):
		if dictionary['PUMP'][i]==1:
			tPumpPre = tPumpPre+((dictionary['elapsed_time'][i]-dictionary['elapsed_time'][i-1]))
	for i in range(len(dictionary['PUMP'])/2, len(dictionary['PUMP'])):
		if dictionary['PUMP'][i]==1:
			tPumpPost = tPumpPost+((dictionary['elapsed_time'][i]-dictionary['elapsed_time'][i-1]))		
	tPump=np.array([tPumpPre,tPumpPost])
	return tPump

def PumpActiveCycle(dictionary):
	CycleCountPre=0
	CycleCountPost=0
	for i in range(1,len(dictionary['PUMPcycles'])/2):
		dC=dictionary['PUMPcycles'][i]-dictionary['PUMPcycles'][i-1]
		if dC>0:
			CycleCountPre=CycleCountPre+dC
	for i in range(len(dictionary['PUMPcycles'])/2, len(dictionary['PUMPcycles'])):
		dC=dictionary['PUMPcycles'][i]-dictionary['PUMPcycles'][i-1]
		if dC>0:
			CycleCountPost=CycleCountPost+dC	
	CycleCount=np.array([CycleCountPre,CycleCountPost])
	return  CycleCount

def atTime(dictionary, time, instrument):
	index=0
	while time >= dictionary['elapsed_time'][index]:
		index=index+1
	if time == dictionary['elapsed_time'][index-1]:
		return dictionary[instrument][index-1]
	else:
		x=(dictionary[instrument][index]+dictionary[instrument][index-1])/2
		print index
		return x

def Tdata(dictionary, instrument):
	Tmin=min(dictionary[instrument])
	Tmax=max(dictionary[instrument])
	Tavg=np.mean(dictionary[instrument])
	return (Tmin,Tavg,Tmax)

def main(dictionary, edge=step):
	if not dictionary['slowDAQ']['loaded']:
		print "Failed to load slowDAQ dictionary, process terminated."
		emptyTempData=np.zeros((8,3))
		emptypump=np.array([0,0])
		emptytime=np.array([0])
		emptybin=np.zeros((9,max[edge]))
		emptydict={'PumpActiveCycle':emptypump, 'PumpActiveTime':emptypump, 'TempData':emptyTempData, 'tEvent':emptytime, 'tGood':emptytime, 'PressureBins':emptybin, 'PressureEdge':emptybin}
		return emptydict
	temp=TrimAll(dictionary['slowDAQ'])
	BinAll(temp,edge)
	TempData=np.ndarray(shape=(8,3), dtype=float, order='C')
	for i in range(1,9):
		TempData[i-1]=Tdata(dictionary['slowDAQ'], 'T'+str(i))
	PBins=np.ndarray(shape=(9,len(temp['BinTimetrimPT1'])), dtype=float, order='C')
	for i in range (1,10):
		PBins[i-1]=temp['BinTimetrimPT'+str(i)]
	PAC=PumpActiveCycle(dictionary['slowDAQ'])
	PAT=PumpActiveTime(dictionary['slowDAQ'])
	EventTime=tEvent(dictionary['slowDAQ'])
	GoodTime=tGood(dictionary['slowDAQ'])
	DataTrim={'PumpActiveCycle':PAC, 'PumpActiveTime':PAT, 'TempData':TempData, 'tEvent':EventTime, 'tGood':GoodTime, 'PressureBins':PBins, 'PressureEdge':temp['BintrimPT4']}
	print ShowIndex(DataTrim)
	return DataTrim
'''
d = ReadFile('13/slowDAQ_0.txt')

d['slowDAQ']=ReadFile('13/slowDAQ_0.txt')
d['slowDAQ']['loaded'] = True

r=main(d)
'''