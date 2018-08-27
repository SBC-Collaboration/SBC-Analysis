
import pdb
import numpy as np


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
        header = [next(read_in).strip() for x in range(6)]
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
                data_type = 'S1'
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
                        data_type = 'S1'
                    if '(' not in var:
                        # pdb.set_trace()
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
    temp = tempdic(dictionary['slowDAQ_0'])
    P = np.zeros(shape=(2, len(timestamp)), dtype=float, order='C')
    for i in range(len(timestamp)):
        x = atTime(temp, timestamp[i], 'PT4')
        P[0][i] = x
        P[1][i] = timestamp[i]
    return P

'''    
l=[-10,-2,-1,-0.02, 0]

d=ReadFile('slowDAQ_0.txt')

print main(d,l)
'''