import pickle
import scipy.io
import sys
import os

""" Convert pickle file to matlab file. One pickle per file and matlab variable name is the
    basename of the pickle file name.
    Usage: python pickle2mat.py <pickle file name> <mat file name>
"""

def main(args):
    pk_file = args[1]
    mat_file = args[2]
    data = pickle.load(open(pk_file, "rb"))
    varname = os.path.splitext(os.path.basename(pk_file))[0]
    scipy.io.savemat(mat_file, mdict={varname: data})


if __name__ == '__main__':
    main(sys.argv)
