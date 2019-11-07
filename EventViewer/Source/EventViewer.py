import os
import re

## Directories

# raw_directory: where the data for an individual dataset is stored
# base_directory: where the datasets are stored, created from raw_directory
# scan_directory: where the handscans are stored
# reco_directory: where the .npy files (raw_events.npy, reco_events.npy) and merged_all.txt
# ped_directory: where the ped code is stored on the machine
# config_file_directory: where the configuration files are placed, should be a folder in the same directory as ped_directory
### also change self.ped_config_file_path_var to append desired initial dataset rather than default

# #for running on COUPP machine
#raw_directory = '/bluearc/storage/30l-16-data/'
#base_directory, end = re.compile('\\w*-\\w*-data').split(self.raw_directory)
#scan_directory = '/coupp/data/home/coupp/scan_output_30l-16/'
#reco_directory = '/coupp/data/home/coupp/recon/current/30l-16/output/'
#ped_directory = '/coupp/data/home/coupp/PEDsvn/'
#config_file_directory = os.path.join(self.ped_directory, 'configs')

###for running on local machine, change these based on local file location to set the correct data directories and initial dataset
raw_directory = '/bluearc/storage/SBC-17-data/'
base_directory, end = re.compile('\\w*-\\w*-data').split(raw_directory)
scan_directory = '/coupp/data/home/coupp/scan_output_SBC-17/'
# self.reco_directory = '/pnfs/coupp/persistent/grid_output/SBC-17/output/'
reco_directory = ""
ped_directory = '/nashome/j/jgresl/Desktop/EventViewer/PEDsvn'
config_file_directory = os.path.join(ped_directory, 'configs')
# #  for running on local machine, change these based on local file location to set the correct data directories and initial dataset
# self.raw_directory = './30l-16-data'
# self.base_directory, end = re.compile('\\w*-\\w*-data').split(self.raw_directory)
# self.scan_directory = './scan_output_30l-16'
# self.reco_directory = './30l-16-data/output'
# self.ped_directory = "./"
# self.config_file_directory = './configs'