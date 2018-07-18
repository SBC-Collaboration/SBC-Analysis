import numpy as np
from functools import reduce
import collections

endianflag = np.uint32(0x01020304)


def WriteBinaryNtupleFile(filename, data,
                          rowdef=0, initialkeys=[], drop_first_dim=False):

    with open(filename, 'wb') as fid:

        # check whether data is a dict or a list of dicts
        if type(data) is dict or type(data) is collections.OrderedDict:
            keylist = data.keys()
            if all([data[k].shape[0] == data[keylist[0]].shape[0]
                   for k in keylist]):
                drop_first_dim = True
            data = [data]

        # return if list is empty
        if (type(data) is list) and (len(data) < 1):
            print("WriteBinaryNtupleFile::Error: Input list is empty, no data to write.")
            return

        # Now get some info on keys and stuff
        whole_keylist = data[0].keys()
        keylist = initialkeys[:]
        keylist.extend(filter(lambda key:
                              key not in initialkeys, whole_keylist))
        keylist = tuple(keylist)
        if drop_first_dim:
            keyshapes = {key: data[0][key].shape[1:] for key in keylist}
            for key in keylist:
                if len(keyshapes[key]) == 0:
                    keyshapes[key] = (1,)
        else:
            keyshapes = {key: data[0][key].shape for key in keylist}

        keyshape_strings = {key:
                            reduce(lambda a, s: a + str(s) + ',',
                                   keyshapes[key][-1::-1], '')
                            for key in keylist}
        keytypes = {key: data[0][key].dtype for key in keylist}

        headerstring = reduce(lambda a, k:
                              a + k + ';' +
                              str(keytypes[k]) + ';' +
                              keyshape_strings[k][:-1] + ';',
                              keylist, '')
        # OK, now time to write some files
        endianflag.tofile(fid)
        np.uint16(len(headerstring)).tofile(fid)
        fid.write(headerstring.encode('ascii'))

        np.int32(rowdef).tofile(fid)

        # And now the header is done, time to write the actual file
        for d in data:
            if drop_first_dim:
                for ix in range(d[keylist[0]].shape[0]):
                    for key in keylist:
                        if d[key][ix].size > 1:
                            d[key][ix].ravel(order='C').tofile(fid)
                        else:
                            d[key][ix].tofile(fid)
            else:
                for key in keylist:
                    if d[key].size > 1:
                        d[key].ravel(order='C').tofile(fid)
                    else:
                        d[key].tofile(fid)
    return
