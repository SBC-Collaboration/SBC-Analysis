
import pdb
import numpy as np


def tempdic(dictionary):
    latch_lowtohigh = np.nonzero(np.diff(dictionary['TriggerLatch']) == 1)
    dt = np.diff(dictionary['elapsed_time'][0:2])
    timestart = dt * latch_lowtohigh[-1] * (-1)
    timeend = (len(dictionary['elapsed_time']) - latch_lowtohigh[-1]) * dt
    zerotime = np.arange(timestart, timeend, dt, dtype=float)
    temp = {'timeaxis': zerotime}
    for i in range(1, 10):
        temp.update({'PT' + str(i): dictionary['PT' + str(i)]})
    return temp


def atTime(dictionary, time, instrument):
    index = 0
    while time >= dictionary['timeaxis'][index]:
        index = index + 1
    if time == dictionary['timeaxis'][index - 1]:
        return dictionary[instrument][index - 1]
    else:
        x = (dictionary[instrument][index] + dictionary[instrument][index - 1])/2
        return x


def ShowIndex(dictionary):
    for key in dictionary:
        print key + str(dictionary[key].shape) + str(dictionary[key].dtype)


def main(dictionary, timestamp):
    if not dictionary['slowDAQ']['loaded']:
        print "Failed to load slowDAQ dictionary, process terminated."
        empty = np.zeros(shape=(2, len(timestamp)), dtype=float, order='C')
        for i in range(len(timestamp)):
            empty[1][i] = timestamp[i]
        return empty
    temp = tempdic(dictionary['slowDAQ'])
    P = np.zeros(shape=(2, len(timestamp)), dtype=float, order='C')
    for i in range(len(timestamp)):
        x = atTime(temp, timestamp[i], 'PT4')
        P[0][i] = x
        P[1][i] = timestamp[i]
    return P
