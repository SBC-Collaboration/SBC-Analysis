##

import os
import time
import copy
import gc

import numpy as np
import matplotlib.pyplot

import SBCcode

def add_mask(d, new_mask, do_copy=False):
    # Inputs:
    #   d: A dictionary where the values are boolean masks
    #   new_mask: The new_mask. Will apply logical and between this mask and all values of d
    #   do_copy: Boolean telling us whether or not to copy d before modifying
    # Outputs
    if do_copy:
        d = copy.deepcopy(d)
    for source in d.keys():
        d[source] = d[source] & new_mask
    return d


def display_masks(d):
    # Inputs:
    #   d: A dictionary where the values are boolean masks
    # Outputs: None. Just prints some useful information
    for source, mask in d.items():
        print("{}: {} elements after masking".format(source, sum(mask)))
    return

if __name__ == "__main__":
    print("If you're running this yourself and get a warning about an invalid comparison,\n"
          "it's because sometimes comparing nans with values will raise that warning. You\n"
          "can disable warning with 'python -W ignore scriptname.py' when you run it.\n\n\n")
    t0 = time.time()

    recon_path = "/pnfs/coupp/persistent/grid_output/SBC-17/output"
    alt_recon_path = "/pnfs/coupp/persistent/grid_output/SBC-17-T0Test3/output"
    raw_path = "/bluearc/storage/SBC-17-data/"
    sources = {"Bi207": [1207],
               "BiBe": [19207],
               "Cf252": [10252, 20252],
               "bkg": [0]}

    ## 0: Load _all binary files
    timing_path      = os.path.join(recon_path, "TimingAnalysis_all.bin")
    acoustic_path    = os.path.join(alt_recon_path, "AcousticAnalysis_all.bin")
    history_path     = os.path.join(recon_path, "HistoryAnalysis_all.bin")
    ev_analysis_path = os.path.join(recon_path, "EventAnalysis_all.bin")
    xyz_path         = os.path.join(recon_path, "SimpleXYZ_all.bin")

    timing_data      = SBCcode.read_bin(timing_path)
    acoustic_data    = SBCcode.read_bin(acoustic_path)
    history_data     = SBCcode.read_bin(history_path)
    ev_analysis_data = SBCcode.read_bin(ev_analysis_path)
    xyz_data         = SBCcode.read_bin(xyz_path)

    t_read = time.time()
    print("\tTook {t:.4f}s to read files.".format(t=t_read-t0))
    print()

    ## 1: Generate dictionary of masks
    print("\tGenerating Event_run_type Masks...")
    masks = {}
    for source, runtype in sources.items():
        mask = np.zeros(ev_analysis_data["Event_run_type"].shape, dtype=bool)
        for rt in runtype:
            mask = mask | (ev_analysis_data["Event_run_type"] == rt)
        masks[source] = mask

    ## 2. Pressure Cutoff
    print("\tGenerating pressure masks using pt6...", end="")
    max_pressure = 25 # psi
    add_mask(masks, history_data['EventPressure'][0:, 5] < max_pressure)
    print("\t\tSum = {}".format(sum(history_data["EventPressure"][0:, 5] < max_pressure)))

    ## 3. Bubble existence mask
    print("\tPerforming 'did-a-bubble-happen' masks...", end="")
    z = xyz_data["bubZ"]
    z_mask = np.invert(np.isnan(z))
    add_mask(masks, z_mask) # So np.nan -> False value in our mask. Everything else is True
    print("\tSum = {}".format(sum(z_mask)))

    ## 4. Timing Lag Cutoff
    # We should only look at events where their timing lag is within
    # some criteria, based on the bubble's z-position. For 30-psi cutoff,
    # y = (20.65087x - 41.03) * 1e-6
    # and our 98% band was +- 49.166419 *1e-6
    print("\tPerforming time lag masks...", end="")
    lag = timing_data["PMTmatch_lag"][:,1]
    slope = 20.65087
    intercept = -41.03
    max_lag = 49.166419e-6 # ~50 micro seconds
    time_cutoff_mask = (slope*z+intercept-max_lag <= lag*1e6) | (slope*z+intercept+max_lag >= lag*1e6)
    add_mask(masks, time_cutoff_mask)
    print("\t\t\t\tSum = {}".format(sum(time_cutoff_mask)))
    print()
    display_masks(masks)
    t_masks = time.time()
    print("\tTook {t:.4f}s to create masks for sources.".format(t=t_masks - t_read))

    ## 5. Get nPHE for each event and save to a numpy file since it'll probably take a while to load.
    gc.disable()
    nPHE = {}
    cached_events = {}
    gc_total_time = 0
    n_completed = 0
    #cached_files = np.array(["nPHE_with_masks.npy", ])
    save_to = "nPHE_with_masks.npy"
    #cached_files_to_check = cached_files[np.where([os.path.isfile(f) for f in cached_files])]
    #cache_storage = [np.load(f) for f in cached_files_to_check]
    n_total = sum([sum(v) for v in masks.values()])
    for source, mask in masks.items():
        runids = timing_data["runid"][mask]
        events = timing_data["ev"][mask]
        n_shape = runids.shape
        nPHE_source = np.zeros(n_shape) + np.nan
        for n in range(n_shape[0]):
            tgc = time.time()
            gc.collect()
            gc_total_time += (time.time()-tgc)
            runidstr = str(runids[n][0])+"_"+str(runids[n][1])
            print("{evstr} -- {source} -- {pSource}% of Source -- {pTotal} of Total". \
                  format(evstr=runidstr+"_"+str(events[n]),
                         source=source,
                         pSource=100*n/n_shape[0],
                         pTotal=100*n_completed/n_total))


            # Because this takes so long, create a cache so we can lookup later in case things overlap,
            # however it's unlikely that it ever will since you can't have two sources. Generality doesn't
            # hurt, though! Alternatively, in case we change our masks we can check to see if we've already
            # loaded nPHE somewhere.
            # if runidstr+"_"+str(events[n]) in cached_events.keys():
            #     print("\tLooking up {} in local cache!".format(runidstr+str(events[n])))
            #     nPHE[n] = cached_events[runidstr+"_"+str(events[n])]
            #     continue
            # for cs in cache_storage:
            #     print("\tLooking up {} in file cache at {}".format(runidstr + str(events[n]), save_to))

            pmtpa_path = os.path.join(recon_path, runidstr, "PMTpulseAnalysis_{runid}.bin".format(runid=runidstr))
            try:
                nPHE_source[n] = SBCcode.read_bin(pmtpa_path, max_file_size=750)["pmt_nphe"][events[n]]
                #cached_events[runidstr+"_"+str(events[n])] = nPHE[n]
            except OSError:
                pass
            n_completed += 1
        nPHE[source] = nPHE_source
    t_nphe = time.time()

    print("\tTook {:.4f}s to create lists of nPHE. {:.4f}s was spent garbage collecting.".\
          format(t_nphe-t_masks, gc_total_time))
    np.save(arr=nPHE, file=save_to)
