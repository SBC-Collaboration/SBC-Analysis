# Analyze Dytran

from scipy import optimize
# import matplotlib.pyplot as plt
import numpy as np


# from LowPass import lowPass
# from PiecewiseFit import LinearPiecewise
# from PiecewiseFit import Fit
# from SBCcode.DataHandling.ReadBinary import ReadBlock as rb
# from SBCcode.DataHandling.ReadText import ReadFile as rt
# import ipdb

def dytranAnalysis(event):
    # Inputs:
    #   event: Event data from GetEvent
    # Outputs: A dictionary of values we're interested in. For a full list, see 'default_output' below.
    default_output = dict(dytran_fitparams=np.zeros(4, dtype=np.float64)-1.)
    try:
        if not event['fastDAQ']['loaded'] or 'Dytran' not in event['fastDAQ'].keys():
            return default_output
        c_time = np.logical_and(event['fastDAQ']['time'][0] / 2 <= event['fastDAQ']['time'],
                                event['fastDAQ']['time'][-1] / 2 >= event['fastDAQ']['time'])
        g, _ = optimize.curve_fit(LinearPiecewise, event['fastDAQ']['time'][c_time], event['fastDAQ']['Dytran'][c_time],
                                  p0=[0, 0.15, 0, 2], bounds=([-0.12, -5, -100, -100], [0, 5, 100, 100]))
        # ipdb.set_trace()
        return dict(dytran_fitparams=g)
    except:
        return default_output


# dytranAnalysis1 is the original implementation, which is wrong.
def dytranAnalysis1(event):
    # d = rb(binary)
    # e = rt(text)

    # dy = d['Dytran']*e['Dytran_multiplier'] +e['Dytran_offset']

    default_output = dict(dytran_fitparams=np.zeros(4, dtype=np.float64))

    if not event['fastDAQ']['loaded'] or 'Dytran' not in event['fastDAQ'].keys():
        return default_output

    dy = event['fastDAQ']['Dytran']
    dytran = lowPass(dy)

    shortx = dytran[0][100:500]
    shorty = dytran[1][100:500]

    g = Fit(shortx, shorty)

    return dict(dytran_fitparams=g)


def plotTheDytran(event):
    if not event['fastDAQ']['loaded'] or 'Dytran' not in event['fastDAQ'].keys():
        return "Sorry, we can't seem to find the data you were looking for."
    dy = event['fastDAQ']['Dytran']
    dytran = lowPass(dy)
    shortx = dytran[0][100:500]
    shorty = dytran[1][100:500]

    g = Fit(shortx, shorty)

    xs = np.linspace(100000, 500000, 100000)


# plt.plot(dytran[0],dytran[1],"yo")
# plt.plot(xs,LinearPiecewise(xs,*g))

def lowPass(dytran):
    newXPoints = [x for x in range(1000, 700001,
                                   1000)]  # Array of new x axis points to be used, does not start at 0 to avoid 0,0 point

    averages = [None for x in range(700)]  # Array of average y values to be used

    counter = 0
    temp = []

    # reshaped_piezo = np.reshape(piezo, (700,1000), order='C')
    # Then take standard deviation of the reshaped rows?
    # Normalize by that?

    for h in range(0, 1000):
        temp.append(dytran[h])  # Edit this to take out the [0] once ReadBinary is properly updated

    standardDev = np.std(temp)

    for j in range(1, 700001):  # Starts at 1 to avoid a 0,0 point

        if j % 1000 != 0:  # Makes it so it adds a thousand values, then the else statement resets the counter and adds the average to the array
            counter += dytran[j]

        else:
            # averages.append(counter/1000/standardDev)
            averages[int((j / 1000) - 1)] = counter / 1000 / standardDev
            counter = 0

    newXPoints = np.array(newXPoints, dtype=float)
    averages = np.array(averages, dtype=float)
    return [newXPoints, averages]  # Returns an array of x values and the corresponding y values which are now averages


def LinearPiecewise(x, x0, y0, m1,
                    m2):  # x = x data, x0 and y0 are the coordinates of the 'turning point', m1 and m2 are the slopes of the 1st and second segments of the line, respectively
    return np.piecewise(x, [x < x0], [lambda x: m1 * x + y0 - m1 * x0, lambda x: m2 * x + y0 - m2 * x0])


# Returns an array of y values that have been computed according to the given expressions
# Since only one condition ([x <x0]) is given, the first lambda is executed when the condition is true, the second lambda is executed when the condition is false

def Fit(xdata, ydata):  # Input the x and y data that needs to be fitted, !!!Make sure x,y are np.arrays, not just lists
    params, covariance = optimize.curve_fit(LinearPiecewise, xdata,
                                            ydata)  # Returns the values of the parameters that best fit the data. Also returns covariance, but that's not important
    return params
