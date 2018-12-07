# John Gresl
# j.gresl12@gmail.com

import os
from collections import OrderedDict

import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

import SBCcode
from SBCcode.Tools import SBCtools







if __name__ == "__main__":
    raw_directory = "/bluearc/storage/SBC-18-data/"



    collection_dict = OrderedDict()
    collection_dict[0]=("20181207_0", "Everything Off (Just Noise)", "green", 1)



    for value in collection_dict.values():
        if not value[-1]:
            continue
        current_path = os.path.join(raw_directory, value[0])
        events = SBCtools.BuildEventList(current_path)
        for current_event in events:
            pmt_data = SBCcode.get_event(current_path, current_event,
                                         "PMTtraces", max_file_size=1300)["PMTtraces"]
            n_triggers = pmt_data["traces"].shape[0]
            for current_trigger in range(n_triggers):
                cutoff = 600
                timebase = (10e6*np.arange(pmt_data["traces"].shape[2])*pmt_data["dt"][current_trigger, 0])[:cutoff]
                sig1 = (pmt_data["traces"][current_trigger, 0, :] * pmt_data["v_scale"][current_trigger, 0] +\
                       pmt_data["v_offset"][current_trigger, 0])[:cutoff]
                sig2 = (pmt_data["traces"][current_trigger, 1, :] * pmt_data["v_scale"][current_trigger, 1] +\
                       pmt_data["v_offset"][current_trigger, 1])[:cutoff]
                ### To fit the signals to a sin wave
                ## Attempt to fit y=A*sin(B x + C) + D
                # D can be estimated by averaging the entire signal
                D1 = np.mean(sig1)
                D2 = np.mean(sig2)
                # A can be estimated by taking the maximum displacement from the mean
                A1 = np.max(np.abs(sig1)-np.abs(D1))*0.8
                A2 = np.max(np.abs(sig2)-np.abs(D2))*0.8
                # By looking at the data, we know the time is about 200 ns. T=2pi/B -> B=2pi/T
                B1 = np.pi*1.6
                B2 = np.pi*1.6
                # No guess for the phase.
                C1 = 0
                C2 = 0.1
                # We can see our initial guess
                guess1 = A1*np.sin(B1*timebase+C1)+D1
                guess2 = A2*np.sin(B2*timebase+C2)+D2
                ## Now the fit
                fit1 = lambda x: x[0]*np.sin(x[1]*timebase+x[2]) + x[3] - sig1
                fit2 = lambda x: x[0]*np.sin(B2*timebase+x[2]) + x[3] - sig2
                A1fit, B1fit, C1fit, D1fit = leastsq(fit1, [A1, B1, C1, D1])[0]
                A2fit, B2fit, C2fit, D2fit = leastsq(fit2, [A2, B2, C2, D2])[0]

                done_fit1 = A1fit*np.sin(B1fit*timebase+C1fit)+D1fit
                done_fit2 = A2fit*np.sin(B2*timebase+C2fit)+D2fit

                plt.plot(timebase, sig2, "bo")
                plt.plot(timebase, guess2, "r-")
                plt.plot(timebase, done_fit2, "g-")
                plt.show()








    pass