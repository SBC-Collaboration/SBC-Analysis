from SBCcode.DataHandling.GetSBCEvent import GetEvent as ge
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

def trapArea(b1, b2, h):
    return abs(h*((b1+b2)/2.))


def integrate(arr, dx):
    arr = arr.astype(np.float64, casting='safe')
    area = 0
    for ii in range(arr.size-1):
        y1 = arr[ii]
        y2 = arr[ii+1]
        area += trapArea(y1,y2,dx)
    return area

# Basedir is the outermost directory containing the folders with runs in them.
# A sample file directory should look like:
# 	Basedir (whatever it is called)
# 		20160623_0
# 			0
# 			1
# 			2
# 			3
# 		20160627_2
# 			0
# 			1
# 			2
# 	etc...

# You can either set basedir below directly to the absolute system path of where
# the PMT data is, or place this file in the same directory and os.path.join(folder_name)

basedir = os.path.join(os.getcwd(), 'PMT Test Runs', 'New PMT (New Base)')

for run in os.listdir(basedir):
    cwd = os.path.join(basedir, run)

    for event in [d for d in os.listdir(cwd) if os.path.isdir(os.path.join(cwd, d))]:
        ev = ge(cwd, int(event))
        if not ev['PMTtraces']['loaded']:
            continue
        
        shape = np.shape(ev['PMTtraces']['traces'])
        data_array = np.zeros((shape[0],2))
        #average_peak = np.zeros(20)
        num_peaks = 0
        for i in range(shape[0]):
            baseline = np.average(ev['PMTtraces']['traces'][i,1,:50])
            y = ev['PMTtraces']['traces'][i,1,:300] - baseline
            peak_height = y.min()
            peak_index = np.argmin(y)
            # if peak_height < -20:
            #     average_peak = average_peak+y[peak_index-10:peak_index+10]
            #     num_peaks += 1

            peak = y[peak_index-5:peak_index+5]
            area = integrate(peak, ev['PMTtraces']['dt'][0,0])

            data_array[i,0] = area
            data_array[i,1] = abs(peak_height)

        #average_peak /= float(num_peaks)
        # In the following line, data_array[:,0] will plot a histogram of areas
        # Change to data_array[:,1] for a histogram of peak heights
        n,bins,patches = plt.hist(data_array[:,0], 50)

        plt.yscale('log', nonposy='clip')
        plt.title(run+' '+event)
        #plt.plot(average_peak)
        plt.show()