import numpy as np
import SBCcode as sbc
import os
import re
from SBCcode.DataHandling.WriteBinary import WriteBinaryNtupleFile as wb
# import ipdb

modules = [
    'AcousticAnalysis_',
    'DytranAnalysis_',
    'EventAnalysis_',
    'HistoryAnalysis_',
    'ImageAnalysis_',
    'TimingAnalysis_',
    'PMTfastDAQalignment_']
# modules = ['PMTpulseAnalysis_']
# modules = ['ImageAnalysis_']
# modules = ['AcousticAnalysis_']
# modules = ['TimingAnalysis_']

# recondir = '/bluearc/storage/recon/devel/SBC-17/output'
recondir = '/pnfs/coupp/persistent/grid_output/SBC-17/output'
merge_dir = '/bluearc/storage/recon/devel/SBC-17/output'

runlist = os.listdir(recondir)
runlist = filter(lambda fn: (not re.search('^\d+_\d+$', fn) is None)
                            and os.path.isdir(os.path.join(recondir, fn))
                            and (len(os.listdir(os.path.join(recondir, fn))) > 0),
                 runlist)
# runlist = ['20170706_6']
# runlist = ['20170621_7','20170625_2']
print(runlist)
# one_piezo_list = [
#     '20170619_3',
#     '20170621_0',
#     '20170621_2',
#     '20170621_3',
#     '20170621_4',
#     '20170621_5',
#     '20170621_6',
#     '20170621_7',
#     '20170621_8',
#     '20170622_0',
#     '20170622_1',
#     '20170622_2',
#     '20170622_3',
#     '20170622_5',
#     '20170622_6',
#     '20170622_7',
#     '20170622_8',
#     '20170622_9',
#     '20170623_0',
#     '20170623_1',
#     '20170623_2']

# merge out by category to save memory
for module in modules:
    # bad_list = [
    #     '20170624_2',
    #     '20170624_4',
    #     '20170625_0',
    #     '20170625_1',
    #     '20170625_2',
    #     '20170704_3',
    #     '20170704_4',
    #     '20170705_0',
    #     '20170705_1',
    #     '20170705_2',
    #     '20170706_5',
    #     '20170713_3',
    #     '20170713_4',
    #     '20170713_5',
    #     '20170714_0',
    #     '20170714_1',
    #     '20170714_2',
    #     '20170715_0',
    #     '20170715_1',
    #     '20170715_2',
    #     '20170715_4',
    #     '20170716_0',
    #     '20170716_1',
    #     '20170716_2',
    #     '20170716_3',
    #     '20170716_5',
    #     '20170716_6',
    #     '20170716_7',
    #     '20170717_0']
    # if key == 'AcousticAnalysis_':
    #     bad_list += [
    #     '20170621_1',  '20170622_4',  '20170624_3',  '20170711_13', '20170706_6', '20170708_2', '20170719_11']
    #     bad_list = []

    # if key == 'ImageAnalysis_':
    #     bad_list = ['20170626_9', '20170703_3', '20170707_4']
    # elif key == 'DytranAnalysis_':
    #     bad_list = [
    #     '20170622_9',
    #     '20170624_4',
    #     '20170625_0',
    #     '20170625_1',
    #     '20170704_3',
    #     '20170704_4',
    #     '20170705_0',
    #     '20170705_1',
    #     '20170705_2',
    #     '20170706_5']
    # elif key == 'EventAnalysis_':
    #     bad_list = ['20170621_1'  '20170622_4'  '20170624_3']
    # elif key == 'PMTfastDAQalignment_':
    #     bad_list = ['20170621_1'  '20170622_4'  '20170624_3']
    bad_list = []
    print("Loading " + module)

    merge_out = []
    shapes0 = []

    for runname in runlist:

        if runname in set(bad_list):
            print(runname + ' is in bad_list')
            continue

        runid_str = runname.split('_')
        runid = np.int32(runid_str)
        runsn = runid[0] * 1000 + runid[1]

        if (runsn >= 20170619003) and (runsn < 20170901000):
            fpath = os.path.join(recondir, runname, module + runname + '.bin')
            if os.path.exists(fpath):
                if os.stat(fpath).st_size > 0:
                    data = sbc.read_bin(fpath)

                    # # check array sizes
                    # shapes = [data[x].shape for x in data.keys()]
                    # if len(shapes0) < 1:
                    #     shapes0 = shapes
                    # print(runname + "\t" + str(shapes))

                    # Pad 0's to fields without Piezo2
                    if module == 'AcousticAnalysis_' and len(data['piezo_list'].shape) == 1:
                        size = [data['piezo_list'].shape[0], 2]

                        tmp = data['piezo_list']
                        data['piezo_list'] = np.zeros(size, dtype=np.int32)
                        data['piezo_list'][:, 0] = tmp

                        tmp = data['bubble_t0']
                        data['bubble_t0'] = np.zeros(size, dtype=np.float64)
                        data['bubble_t0'][:, 0] = tmp

                        tmp = data['peak_t0']
                        data['peak_t0'] = np.zeros(size, dtype=np.float64)
                        data['peak_t0'][:, 0] = tmp

                        size = list(data['piezoE'].shape)
                        size[1] += 1

                        tmp = data['piezoE']
                        data['piezoE'] = np.zeros(size, dtype=np.float64)
                        # ipdb.set_trace()
                        data['piezoE'][:, 0, :, :] = tmp[:, 0, :, :]

                    if module == 'TimingAnalysis_' and len(data['PMTmatch_t0'].shape) == 1:
                        var_names = ['CAMstate', 'PMTmatch_area', 'PMTmatch_area_nobs', 'PMTmatch_baseline', 'PMTmatch_baserms', 'PMTmatch_coinc', 'PMTmatch_ix', 'PMTmatch_lag', 'PMTmatch_max', 'PMTmatch_min', 'PMTmatch_pulse_area', 'PMTmatch_pulse_height', 'PMTmatch_pulse_t10', 'PMTmatch_pulse_t90', 'PMTmatch_pulse_tend', 'PMTmatch_pulse_tpeak', 'PMTmatch_pulse_tstart', 'PMTmatch_t0', 'nPMThits_fastdaq', 'nVetohits_fastdaq', 't_nearestPMThit', 't_nearestVetohit']
                        for var_name in var_names:
                            if len(data[var_name].shape) == 1:
                                data[var_name] = np.stack((data[var_name],
                                                       np.zeros(data[var_name].shape, data[var_name].dtype)), axis=1)
                            elif len(data[var_name].shape) > 1:
                                data[var_name] = np.concatenate((data[var_name],
                                                           np.zeros(data[var_name].shape, data[var_name].dtype)),
                                                          axis=1)
                    if module == 'TimingAnalysis_': # fix int32/int64 problem
                        var_name = 'PMTmatch_ix'
                        data[var_name] = np.int64(data[var_name])

                    shapes = [(x, data[x].dtype, data[x].shape) for x in data.keys()]
                    if len(shapes0) < 1:
                        shapes0 = shapes
                    print(runname + "\t" + str(shapes))
                    # ipdb.set_trace()
                    merge_out.append(data)
                else:
                    print("zero size file: " + fpath)
            else:
                print("nonexis file: " + fpath)

    merge_name = 'all'
    rowdef = 1

    if module in set(['PMTpulseAnalysis_', 'PMTpheAnalysis_']):
        rowdef = 7
    if module in set(['HumanGetBub_']):
        rowdef = 8
    print("Writing " + module)
    wb(os.path.join(merge_dir, module + merge_name + '.bin'), merge_out,
       rowdef=rowdef, initialkeys=['runid', 'ev'], drop_first_dim=True)
