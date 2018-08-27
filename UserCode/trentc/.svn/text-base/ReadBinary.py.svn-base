#Author: Trent Cwiok

import numpy as np
import sys
from collections import OrderedDict
import pdb
import matplotlib.pyplot as plt

np.set_printoptions(threshold = np.nan)

#Get the Endianness of the system and set a default file endianness
#If the two are equivalent, the bytes will have to be flipped to be read
system_endianness = sys.byteorder
file_endianness = 'little'

def ReadBlock(file_name):
	'''
	This function takes in a binary data file, and reads it in byte-wise, then recasts each variable
	to the proper data type, and stores it in a dictionary to be returned.
	'''
	variables_dict = OrderedDict()
	possible_data_types = {'char': 8, 'int8': 8, 'int16': 16, 'int32': 32, 'int64': 64, 'uint8': 8,
						   'uint16': 16, 'uint32': 32, 'uint64': 64, 'single': 32, 'double': 64, 'float128': 128}

	#Open file here
	with open(file_name, "rb") as read_in:
		#Check the Endianness flag of the block
		endianness = np.fromfile(read_in, dtype = np.uint32, count = 1)
		#0x01020304 = 16909060 in base 10
		if endianness != 16909060:
			file_endianness = 'big'
			print 'File Endianness Changed'

		#Check the length of the header string
		header_len = np.fromfile(read_in, dtype = np.uint16, count = 1)

		#Read in the header string and split it into different variable names
		header = np.fromfile(read_in, dtype = 'S'+str(header_len[0]), count = 1)
		header_str = header[0]
		header_components = header_str.split(';')

		#Create parallel arrays to store the variables types and numbers, then populate them
		meta_data = OrderedDict()
		for variable in xrange(0, len(header_components), 3):
			if header_components[variable]:
				meta_data[header_components[variable]] = [header_components[variable+1], header_components[variable+2]]

		#Get the number of data lines in the block
		num_lines = np.fromfile(read_in, dtype = np.int32, count = 1)[0]

		#Calculate the number of bytes per line
		bytes_per_line = 0
		for key in meta_data:
			if len(meta_data[key][1].split(',')) == 1:
				bytes_per_line += (possible_data_types[meta_data[key][0]] * int(meta_data[key][1]) / 8)
			else:
				temp_size = 1
				sizes = meta_data[key][1].split(',')  
				for ele in sizes:
					temp_size *= int(ele)
				bytes_per_line += (possible_data_types[meta_data[key][0]] * temp_size / 8)

		if num_lines <= 0:
			start_of_data = read_in.tell()
			read_in.seek(0,2)
			end_of_data = read_in.tell()
			read_in.seek(start_of_data, 0)
			num_lines = (end_of_data - start_of_data)/bytes_per_line
		
		for variable in xrange(0, len(header_components), 3):
			if header_components[variable]:
				variables_dict[header_components[variable]] = None

		uint8_buffer = np.fromfile(read_in, dtype = np.uint8, count = -1)
		uint8_buffer = np.reshape(uint8_buffer, (num_lines, bytes_per_line), order = 'C')
		
		start = 0
		for key in variables_dict:
			#If data shape is simple, the width is read directly
			if len(meta_data[key][1].split(',')) == 1:
				width = int(meta_data[key][1])
				if width==1:
					sizes = num_lines,
				else:
					sizes = (num_lines, width)
			#Else, boolean is flipped for later reshaping to occur, and the number of iterations is the product of each dimension
			else:
				width = 1
				sizes = meta_data[key][1].split(',')
				for ii in xrange(len(sizes)):
					sizes[ii] = int(sizes[ii])
					width *= sizes[ii]
				sizes.append(num_lines)
				sizes.reverse()
				sizes = tuple(sizes)
			width *= possible_data_types[meta_data[key][0]]/8

			temp = np.zeros((num_lines, width), dtype = np.uint8, order = 'C')
			temp[:,:] = uint8_buffer[::,start:start+width]
			variables_dict[key] = Cast(meta_data[key][0], temp)
			#Uncomment this line to save all data as a double type
			#variables_dict[key] = variables_dict[key].astype(np.float64, copy = False)
			
			variables_dict[key] = np.reshape(variables_dict[key], sizes, order = 'C')
			start += width

	return variables_dict


def Cast(variable_name, data):
	'''
	This function takes in the type to be cast to, as well as the data to be cast,
	and returns the data in its newly cast form.  Defaults to return None if the 
	data type is not a valid type to be cast to.
	'''
	if variable_name == 'char':
		return data.view(np.char)

	if variable_name == 'int8':
		return data.view(np.int8)

	if variable_name == 'int16':
		return data.view(np.int16)

	if variable_name == 'int32':
		return data.view(np.int32)

	if variable_name == 'int64':
		return data.view(np.int64)

	if variable_name == 'uint8':
		return data.view(np.uint8)

	if variable_name == 'uint16':
		return data.view(np.uint16)

	if variable_name == 'uint32':
		return data.view(np.uint32)

	if variable_name == 'uint64':
		return data.view(np.uint64)

	if variable_name == 'single':
		return data.view(np.float32)

	if variable_name == 'double':
		return data.view(np.float64)
 
	if variable_name == 'float128':
		return data.view(np.float128)

	else:
		return None