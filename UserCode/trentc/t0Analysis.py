#Author: Trent Cwiok
from SBCcode.DataHandling.ReadBinary import ReadBlock as rb
from matplotlib import pyplot as plt

def main():
	dy = rb('event_output/DytranAnalysis.bin')
	aco = rb('event_output/AcousticAnalysis.bin')

	dy_t0 = dy['dytran_fitparams'][:,0]
	dy_t0 -= 400000
	dy_t0 *= 4e-7
	
	aco_t0 = aco['bubble_t0']

	t0_diff = dy_t0-aco_t0

	plt.plot(aco['bubble_loudness'][:,1], t0_diff, 'bo')
	#plt.axis([0,len(t0_diff), -.1, .1])
	plt.semilogx()
	plt.show()
	return