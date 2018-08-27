import SBCcode as sbc
import ipdb

# data0 = sbc.read_bin('/bluearc/storage/recon/devel/SBC-17/output/AcousticAnalysis_all.bin')

# timing analysis array shape difference due to number of piezos
data1 = sbc.read_bin('/pnfs/coupp/persistent/grid_output/SBC-17/output/20170622_7/TimingAnalysis_20170622_7.bin')
data2 = sbc.read_bin('/pnfs/coupp/persistent/grid_output/SBC-17/output/20170708_3/TimingAnalysis_20170708_3.bin')
# data3 = sbc.read_bin('/pnfs/coupp/persistent/grid_output/SBC-17/output/20170802_2/TimingAnalysis_20170802_2.bin')
data3 = sbc.read_bin('/nashome/j/jzhang/SBCcode/UserCode/jzhang/20170802_2/TimingAnalysis_20170802_2.bin')
for key in data1.keys():
    print(key, data1[key].dtype, data1[key].shape, data2[key].dtype, data2[key].shape, data3[key].dtype, data3[key].shape)
# ipdb.set_trace()