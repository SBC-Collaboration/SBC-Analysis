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

#np.set_printoptions(threshold=np.nan)


def ReadFile(filename):
    '''
    This function takes in the file's name to be read, as well as
    the size of the header containing metadata, and returns a
    dictionary where the individual instrument names are the keys,
    and the associated values are arrays of the value recorded by
    the instrument at each timestep. Returns None when errors.
    '''
    instrument_dict = {}

    # Get the number of lines in the file for appropriate array construction
    with open(filename, 'r') as f:
        for i, l in enumerate(f):
            pass
        # Subtract five to account for the header size
        num_lines = i - 5

    with open(filename, 'r') as read_in:
        # Read in the file header
        header = [next(read_in).strip() for x in range(6)]
        instrument_list = header[1].split()

        # Create keys in the dictionary for each variable and
        # initialize an empty array to store the data
        for var in instrument_list:
            shape = num_lines
            data_type = 'double'
            if '(' in var:
                dimensions = var[var.find("(") + 1:var.find(")")]
                shape = str(num_lines) + ',' + dimensions
                shape = tuple(map(int, shape.split(',')))
            if var == 'run':
                data_type = 'S1'
            instrument_dict[var] = np.zeros(shape, dtype=data_type, order='C')

        # This loop runs until the file end is reached, taking
        # in the data in each line, parsing it, and storing it
        # to the appropriate key in the dictionary
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

                    # NOTE: This part of the loop is meant to handle
                    # variables with multi-dimensional entries, but is
                    # untested due to the absences of any test file
                    else:
                        shape = var[var.find('(') + 1:var.find(')')]
                        shape = tuple(map(int, shape.split(',')))
                        iterations = 1
                        for ele in iterations:
                            iterations *= ele

                        temp_buffer = np.zeros(counter,
                                               dtype=data_type, order='C')
                        counter = 0
                        while counter < iterations:
                            temp_buffer[counter] = data[front]
                            front += 1
                            counter += 1
                        temp_buffer = np.reshape(temp_buffer, shape, order='C')
                        instrument_dict[var][index] = temp_buffer

            index += 1

    return instrument_dict


def InstrumentOutput(dictionary, instrument):
    '''
    Given a dictionary containing instrument:data as key:value
    pairs, and a specific instrument in that dictionary, prints
    the list of values recorded by that instrument.
    '''
    instrument = str(instrument)

    try:
        for value in range(len(dictionary[instrument])):
            print(dictionary[instrument][value])
        return

    except KeyError:
        print("Instrument" + str(instrument) +
              " does not exist in that dictionary, please try again.")
        return
    except TypeError:
        print("The first argument is not the proper data structure." +
              " Please pass a dictionary.")
        return


def PlotInstrument(dictionary, instrument):
    '''
    Given a dictionary containing instrument:data as key:value
    pairs, and a specific instrument in that dictionary, plots
    that data against the elapsed time and displays the results.
    '''
    instrument = str(instrument)

    try:
        # plt.plot(dictionary['elapsed_time'], dictionary[instrument])
        # plt.show()
        return

    except KeyError:
        print("Instrument" + str(instrument) +
              " does not exist in that dictionary, please try again.")
        return
    except TypeError:
        print("The first argument is not the proper data structure." +
              " Please pass a dictionary.")
        return
