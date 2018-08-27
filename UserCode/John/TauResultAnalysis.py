## John Gresl


import os
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt



RESULTS_DIR = "/nashome/j/jgresl/TauResults/"

def open_tau_results(path):
    out = OrderedDict()
    out["runid"] = []
    out["ev"] = []
    out["tau"] = []
    if path == None:
        return out
    with open(path) as tau_results:
        for line in tau_results.readlines():
            if line[0] != "#":
                data = line.split(",")
                out["runid"].append(data[0].strip())
                out["ev"].append(data[1].strip())
                out["tau"].append(data[2].strip())
    return out

def dictionary_append(d1, d2):
    ## Appends the dictionary d1 to d2. They must have the same keys!
    ## Returns the modified dictionary.
    keys = list(d1.keys())
    for k in keys:
        d2[k].extend(d1[k])
    return d2


def dictionary_print(d):
    for k,v in d.items():
        print("\t{}:   {}".format(k, v))


def load_all_results(root_directory, filter=True):
    # Loads all results and saves them into a directory
    results = open_tau_results(path=None)  # <-- To get the initial dictionary to append to
    for root, dirs, files in os.walk(RESULTS_DIR):
        for sub_dir in dirs:
            for f in os.listdir(os.path.join(RESULTS_DIR, sub_dir)):
                print(f)
                path = os.path.join(RESULTS_DIR, sub_dir, f)
                results = dictionary_append(results, open_tau_results(path))
    if filter:
        results = filter_bad_fits(results)
    return results


def filter_bad_fits(d):
    ## Given a dictionary, d, will go through and eliminate all entries with tau = -1
    keys = list(d.keys())
    number_filtered = 0
    for n, tau in enumerate(d["tau"][:]):
        if tau == "-1":
            for k in keys:
                d[k].pop(n-number_filtered)
            number_filtered += 1
    print("*****Removed {} bad fits*****".format(number_filtered))
    return d


def _to_float(arr):
    ## Tries to turn every element in arr into a float. If it cannot, turns it into a -1 instead
    out = []
    for el in arr:
        try:
            out.append(float(el))
        except ValueError:
            print("Error converting", el)
            out.append(-1)
    return out

def what_tau_to_use(hist_data):
    return hist_data[1][np.argmax(hist_data[0])]+0.5*(hist_data[1][1]-hist_data[1][0])




if __name__ == "__main__":
    from matplotlib import rc


    results = load_all_results(RESULTS_DIR)
    float_taus = _to_float(results["tau"])
    bins = np.linspace(min(float_taus), max(float_taus), 13)
    hist = plt.hist(float_taus, bins=bins)
    plt.xlabel(r"$\tau$ (ms)", fontsize=16)
    plt.ylabel("N", fontsize=16)
    tau = what_tau_to_use(hist)
    plt.xticks(plt.xticks()[0], plt.xticks()[0]*1000)
    print("Number of Taus = {}".format(len(float_taus)))
    print("Best fit Tau = {} ms".format(tau))
    print("Average Tau = {} ms".format(np.average(float_taus)))
    txt = r"# $\tau$: {}"+"\n"+r"Avg $\tau$: {:.3f} ms"+"\n"+r"Peak $\tau$: {:.3f} ms"
    txt = txt.format(len(float_taus), 1e3*np.average(float_taus), 1e3*tau)
    plt.text(0.007, 180, txt, bbox={"facecolor": "red", "alpha": 0.5, "pad": 10}, fontsize=25)

    plt.savefig("Plots/TauHist.png")
    plt.savefig("Plots/TauHist.svg")
    plt.show()



