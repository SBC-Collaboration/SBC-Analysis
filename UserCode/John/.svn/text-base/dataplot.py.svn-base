import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
import scipy

from SBCcode.DataHandling.GetSBCEvent import GetEvent as GE
from SBCcode.DataHandling.ReadBinary import ReadBlock as RB


if __name__ == "__main__":
    runid = "20170709_8"
    event = 84

    raw = "/bluearc/storage/SBC-17-data/"
    rec = "/pnfs/coupp/persistent/grid_output/SBC-17-T0Test3/output"


    a = GE(os.path.join(raw, runid), event, "fastDAQ", "PMTtraces")
    acoustic = RB(os.path.join(rec, runid, "AcousticAnalysis_{runid}.bin".format(runid=runid)))
    print("fastDAQ")
    print("\t", a["fastDAQ"].keys())
    print("\nPMTtraces")
    print("\t", a["PMTtraces"].keys())
    print("\nAcoustic")
    print("\t", acoustic.keys())


    # piezo and led
    #plt.plot(a["fastDAQ"]["time"], a["fastDAQ"]["Piezo1"])
    print("Getting t0 for: {} -- {}".format(acoustic["runid"][event], acoustic["ev"][event]))
    t0 = 1e3*acoustic["bubble_t0"][event][1]
    fig, ax = plt.subplots(nrows=3, ncols=1)
    ax[0].plot((1e3)*a["fastDAQ"]["time"], a["fastDAQ"]["Piezo2"])
    ax[0].plot((1e3)*a["fastDAQ"]["time"], -(a["fastDAQ"]["CAMgate"]*0.21-0.25)-.73)
    ax[0].axvline(t0, color="m", linewidth=3.5)
    ax[0].set_xlim([-0.15*1e3, 0.045*1e3])
    ax[0].set_ylim([-0.7, 0.7])
    ax[0].set_ylabel("Acoustic Amplitude") # Will move and center in powerpoint


    zoom_mask = np.abs(a["fastDAQ"]["time"]-t0*1e-3) <= .6e-3
    ax[1].plot((1e3)*a["fastDAQ"]["time"][zoom_mask], a["fastDAQ"]["Piezo2"][zoom_mask])
    ax[1].axvline(t0, color="m", linewidth=3.5)
    ax[1].set_xlim([-25, -24.3])


    # pmt
    #ax[2].plot(a["PMTtraces"][])
    def return_pmt_time(trace_t0, align_t0):
        # Returns the calibrated time of the PMT trace. Both variables are 2-tuples, with the first element being the
        # "seconds" component, and the 2nd element being the "fractional" component.
        return (trace_t0[0] - align_t0[0]) + (trace_t0[1] - align_t0[1])

    align_data = RB("/pnfs/coupp/persistent/grid_output/SBC-17/output/{runid}/PMTfastDAQalignment_{runid}.bin".\
                    format(runid=runid))
    print("\nAlign Data")
    print("\t", align_data.keys())
    PMT_data = a["PMTtraces"]

    # pmt_data = RB("/pnfs/coupp/persistent/grid_output/SBC-17/output/{runid}/PMTpulseAnalysis_{runid}.bin".\
    #                   format(runid=runid))
    # print("\nPMT Data")
    # print("\t", pmt_data.keys())

    n_max = PMT_data["t0_sec"][:, 0].shape[0]
    times = []
    align_t0_sec = align_data["PMT_trigt0_sec"][event]
    align_t0_frac = align_data["PMT_trigt0_frac"][event]
    for n in range(n_max):
        # Improve this by using a binary search method in the future.
        # Do this now as a proof of concept
        trace_t0_sec = PMT_data["t0_sec"][n, 0]
        trace_t0_frac = PMT_data["t0_frac"][n, 0]
        times.append(np.abs(return_pmt_time((trace_t0_sec, trace_t0_frac), (align_t0_sec, align_t0_frac))))

    times = np.array(times)
    min_indicies = times.argsort()[:10]-14
    min_index = min_indicies[0]
    min_rest = min_indicies[1:]

    for n in min_rest:
        xd = 1e9*np.arange(a["PMTtraces"]["traces"].shape[2]) * a["PMTtraces"]["dt"][n, 0]
        yd_fine = a["PMTtraces"]["traces"][n, 0, :] * a["PMTtraces"]["v_scale"][n, 0] + \
                  a["PMTtraces"]["v_offset"][n, 0]
        ax[2].plot(xd, 1e2*yd_fine, color="k", linewidth=0.3)

    xd = 1e9*np.arange(a["PMTtraces"]["traces"].shape[2]) * a["PMTtraces"]["dt"][min_index, 0]
    yd_fine = a["PMTtraces"]["traces"][min_index, 0, :] * a["PMTtraces"]["v_scale"][min_index, 0] + \
              a["PMTtraces"]["v_offset"][min_index, 0]
    ax[2].plot(xd, 1e2*yd_fine, color="r", linewidth=1)
    ax[2].set_ylabel("PMT Amplitude (mV)")
    ax[2].set_xlim([0,700])
    ax[2].set_ylim([-6, -1])
    #ix = np.argmax(pmt_data["ev"] == event)
    #nphe = pmt_data["pmt_nphe"][ix+min_index]

    # pulse_mask = (100 <= xd) & (xd <= 200)
    # nphe=np.log10(-*np.trapz(yd_fine[pulse_mask], dx=a["PMTtraces"]["dt"][min_index,0]))
    #print("nphe = {}".format(nphe))

    ax[0].set_xlabel("Time (ms)")
    ax[1].set_xlabel("Time (ms)")
    ax[2].set_xlabel("Time (ns)")


    plt.tight_layout()
    plt.savefig("Plots/raw_data.png")
    plt.savefig("Plots/raw_data.svg")
    #plt.show()






