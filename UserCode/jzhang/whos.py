# Mimic matlab whos function
import numpy as np
import collections
import pprint


def whos_long(*args):
    pp = pprint.PrettyPrinter(indent=4)
    for var in args:
        pp.pprint(var)


def who(*args):
    sequentialTypes = [list, tuple]
    for var in args:
        t = type(var)
        print(t)
        print type(var),
        if t == np.ndarray:
            print var.dtype, var.shape,
        elif t in sequentialTypes:
            print len(var)
            whos(var[0])
        elif t == dict or t ==  collections.OrderedDict:
            print
            for key in var.keys():
                print repr(key).rjust(40), " : ", whos(var[key])
        elif t == str:
            print len(var),
    return ''


def whos(var, n=0):
    sequentialTypes = [list, tuple]
    t = type(var)
    # print(t)
    if t == np.ndarray:
        n = n + var.nbytes
    elif t in sequentialTypes:
        for sub_var in var:
            n = n + whos(sub_var, 0)
    elif t == dict or t ==  collections.OrderedDict:
        for key in var.keys():
            # print repr(key).rjust(40), " : ", whos(var[key])
            n = n + whos(var[key], 0)

    return n
