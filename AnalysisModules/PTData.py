# Author: Trent Cwiok
#
# The purpose of this program is to, given a file name containing
# data collected from the instruments attached to the bubble chamber,
# create a data structure capable of storing the instruments and
# their readouts. The information for any given instrument should
# then be easily accessible and printable. Additional functionality
# for plotting any given instrument readout versus time is also
# included.

# import pdb
import numpy as np
# import matplotlib.pyplot as plt
# import SBCcode as sbc


np.set_printoptions(threshold=np.nan)


def DataTrim(dictionary, instrument):
    '''
    Given a dictionary constaining instrument data, it uses
    TriggerLatch to trim away unwanted data points where the
    trigger is not latched
    '''
    tick = 0

    instrument = str(instrument)
    for i in range(len(dictionary['TriggerLatch'])):
        if dictionary['TriggerLatch'][i] == 0:
            tick = tick + 1
    trim_PT = np.zeros(tick)
    trim_time = np.zeros(tick)
    track = 0
    if instrument in dictionary:
        for value in range(len(dictionary[instrument])):
            if dictionary['TriggerLatch'][value] == 0:
                trim_PT[track] = dictionary[instrument][value]
                trim_time[track] = dictionary['elapsed_time'][value]
                track = track + 1

    return trim_PT

# DataTrim(d,'PT4')


def TrimAll(dictionary):

    d = {}
    for i in range(1, 10):
        index = str(i)
        d.update({'trimPT' + str(i): DataTrim(dictionary, 'PT' + index)})
    return d


def ShowIndex(dictionary):
    for key in dictionary:
        print(key + str(dictionary[key].shape) + str(dictionary[key].dtype))


def Pbin(dictionary, instrument, edge):
    '''
    sort pressures into bins of 1 psi
    '''
    if instrument in dictionary:
        Bin = np.histogram(dictionary[instrument],
                           bins=edge, range=(edge[0], edge[-1]))
    else:
        Bin = np.histogram(np.float64([0]) + np.nan,
                           bins=edge, range=(edge[0], edge[-1]))
    BinTime = np.zeros(len(Bin[0]))
    for i in range(len(Bin[0])):
        x = Bin[0][i] * 0.005
        BinTime[i] = x
        dictionary.update({'Bin' + instrument: Bin[1],
                           'Count' + instrument: Bin[0],
                           'BinTime' + instrument: BinTime})
    '''
    print Bin
    plt.hist(dictionary[instrument], bins=step)
    plt.show()
    '''
    return


def BinAll(dictionary, edge):
    for i in range(1, 10):
        index = str(i)
        Pbin(dictionary, 'trimPT' + index, edge)
    return


def tGood(dictionary, instrument='PT4'):
    Pgood = np.float64(0)
    tg = np.float64(0)

    for i in range(1, len(dictionary[instrument])):
        if abs(dictionary[instrument][i] -
               dictionary['PressureSetpoint'][i]) <= 0.3:
            Pgood = Pgood + 1
            tg = tg + ((dictionary['elapsed_time'][i] -
                        dictionary['elapsed_time'][i - 1]))
    return tg


def tEvent(dictionary):
    Tevent = np.float64(0)
    dt = np.diff(dictionary['elapsed_time'][0:2])
    Tevent = np.sum(dictionary['TriggerLatch'] == 0) * dt

    latch_lowtohigh = np.nonzero(np.diff(dictionary['TriggerLatch']) == 1)
    time_of_compression = dictionary['elapsed_time'][latch_lowtohigh[-1]]

    return Tevent


def Pevent(dictionary, instrument):
    latch_lowtohigh = np.nonzero(np.diff(dictionary['TriggerLatch']) == 1)[0]
    pressure_of_compression = np.float64(0) + np.nan
    if instrument in dictionary and \
            latch_lowtohigh.shape[0] > 0 and \
            latch_lowtohigh[-1] >= 40:
        pressure_of_compression = \
            dictionary[instrument][latch_lowtohigh[-1] - 40]

    return pressure_of_compression


def PumpActiveTime(dictionary):
    tPumpPre = np.float64(0)
    tPumpPost = np.float64(0)
    if 'PUMP' in dictionary:
        for i in range(1, int(len(dictionary['PUMP']) / 2)):
            if dictionary['PUMP'][i] == 1:
                tPumpPre = tPumpPre + ((dictionary['elapsed_time'][i] -
                                        dictionary['elapsed_time'][i - 1]))
        for i in range(int(len(dictionary['PUMP']) / 2), len(dictionary['PUMP'])):
            if dictionary['PUMP'][i] == 1:
                tPumpPost = tPumpPost + ((dictionary['elapsed_time'][i] -
                                          dictionary['elapsed_time'][i - 1]))
    tPump = np.array([tPumpPre, tPumpPost], dtype=np.float64)
    return tPump


def PumpActiveCycle(dictionary):
    CycleCountPre = np.int32(0)
    CycleCountPost = np.int32(0)
    if 'PUMPcycles' in dictionary:
        for i in range(1, len(dictionary['PUMPcycles']) / 2):
            dC = dictionary['PUMPcycles'][i] - dictionary['PUMPcycles'][i - 1]
            if dC > 0:
                CycleCountPre = CycleCountPre + dC
        for i in range(len(dictionary['PUMPcycles']) / 2,
                       len(dictionary['PUMPcycles'])):
            dC = dictionary['PUMPcycles'][i] - dictionary['PUMPcycles'][i - 1]
            if dC > 0:
                CycleCountPost = CycleCountPost + dC
    CycleCount = np.array([CycleCountPre, CycleCountPost], dtype=np.int32)
    return CycleCount


def atTime(dictionary, time, instrument):
    index = 0
    while time >= dictionary['elapsed_time'][index]:
        index = index + 1
    if time == dictionary['elapsed_time'][index - 1]:
        return dictionary[instrument][index - 1]
    else:
        x = (dictionary[instrument][index] +
             dictionary[instrument][index - 1]) / 2
        print(index)
        return x


def Tdata(dictionary, instrument):
    Tmin = np.float64(0) + np.nan
    Tmax = np.float64(0) + np.nan
    Tavg = np.float64(0) + np.nan
    if instrument in dictionary:
        Tmin = min(dictionary[instrument])
        Tmax = max(dictionary[instrument])
        Tavg = np.mean(dictionary[instrument])
    return (Tmin, Tavg, Tmax)


def main(dictionary, edge=np.cumsum(np.ones(250)) - 1, targetPT='PT4'):
    if not dictionary['slowDAQ']['loaded']:
        print("Failed to load slowDAQ dictionary, process terminated.")
        emptyTempData = np.zeros((8, 3), dtype=np.float64)
        emptypump = np.zeros(2, dtype=np.float64)
        emptypump_int = np.zeros(2, dtype=np.int32)
        emptytime = np.zeros(1, dtype=np.float64)
        emptybin = np.zeros((9, len(edge) - 1), dtype=np.float64)
        emptyP = np.zeros(9, dtype=np.float64)
        emptydict = {'PumpActiveCycle': emptypump_int,
                     'PumpActiveTime': emptypump,
                     'TempData': emptyTempData,
                     'tEvent': emptytime,
                     'tGood': emptytime,
                     'PressureBins': emptybin,
                     'PressureEdge': edge,
                     'EventPressure': emptyP}
        return emptydict
    temp = TrimAll(dictionary['slowDAQ'])
    BinAll(temp, edge)
    TempData = np.ndarray(shape=(8, 3), dtype=float, order='C')
    for i in range(1, 9):
        TempData[i - 1] = Tdata(dictionary['slowDAQ'], 'T' + str(i))
    PBins = np.ndarray(shape=(9, len(temp['BinTimetrimPT1'])),
                       dtype=float, order='C')
    for i in range(1, 10):
        PBins[i - 1] = temp['BinTimetrimPT' + str(i)]
    PressData = np.zeros(9)
    for i in range(1, 10):
        PressData[i - 1] = Pevent(dictionary['slowDAQ'], 'PT' + str(i))
    PAC = PumpActiveCycle(dictionary['slowDAQ'])
    PAT = PumpActiveTime(dictionary['slowDAQ'])
    EventTime = tEvent(dictionary['slowDAQ'])
    GoodTime = tGood(dictionary['slowDAQ'], targetPT)
    DataTrim = {'PumpActiveCycle': PAC,
                'PumpActiveTime': PAT,
                'TempData': TempData,
                'tEvent': EventTime,
                'tGood': GoodTime,
                'PressureBins': PBins,
                'PressureEdge': temp['Bintrim'+targetPT],
                'EventPressure': PressData}
    # print(ShowIndex(DataTrim))
    # print(DataTrim['PressureBins'])
    return DataTrim


'''
d = ReadFile('13/slowDAQ_0.txt')

d['slowDAQ']=ReadFile('13/slowDAQ_0.txt')
d['slowDAQ']['loaded'] = True

r=main(d)
'''
