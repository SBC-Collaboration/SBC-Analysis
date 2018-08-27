import numpy as np

# execfile('readbinary.py')

# d = ReadBlock('File location/PMTtraces.bin')

# s = ReadBlock('File location/fastDAQ_0.bin')


def PMTandFastDAQalignment(d):
    default_output = dict(PMTt0_sec=np.float64(-1),
                          PMTt0_frac=np.float64(-1))

    if not (d['fastDAQ']['loaded'] and d['PMTtraces']['loaded']):
        return default_output

    output_dict = dict()

    pmttrig_coarse = np.mean(np.reshape(d['fastDAQ']['PMTtrig'][10:],
                                        (-1, 10)), axis=1)

    coarse_dt = dt * 10

    pmthittimes_coarse = np.mean(np.reshape(pmthittimes, (-1,10)),axis=1)

    q_coarse = np.convolve(pmthittimes_coarse, np.fliplr(pmttrig_coarse), mode='full')


    a = np.argmax(q_coarse)
    b = np.argmax(q_coarse) + len(d['fastDAQ']['PMTtrig'])

    q_fine = []

    pmthittimes_length=np.intp(np.floor(len(d['fastDAQ']['PMTtrig'])))
    pmthittimes=np.zeros((pmthittimes_length))
    pmthittimes[-1]=1

    pmt_ix = -2
    while pmt_ix>-len(d['PMTtraces']['t0_sec']):
        pmthittimes[pmthittimes.size+np.intp(np.round((d['PMTtraces']['t0_sec'][pmt_ix,0]-d['PMTtraces']['t0_sec'][-1,0]+d['PMTtraces']['t0_frac'][pmt_ix,0]-d['PMTtraces']['t0_frac'][-1,0])/dt))]=1
        pmt_ix-=1
        if pmthittimes.size+np.intp(np.round((d['PMTtraces']['t0_sec'][pmt_ix,0]-d['PMTtraces']['t0_sec'][-1,0]+d['PMTtraces']['t0_frac'][pmt_ix,0]-d['PMTtraces']['t0_frac'][-1,0])/dt))<0:
                break

    if b > (len(d['fastDAQ']['PMTtrig']) + len(pmthittimes)):
        q_fine = np.convolve(pmthittimes[b - 10:],np.fliplr(d['fastDAQ']['PMTtrig'][len(d['fastDAQ']['PMTtrig']-a-10)]),mode='full')

    if a < 0:
        q_fine = np.convolve(pmthittimes[:b + 10],np.fliplr(d['fastDAQ']['PMTtrig'][len(d['fastDAQ']['PMTtrig']) - a - 10:]),mode='full')

    else:
        q_fine = np.convolve(pmthittimes[a - 10: b + 10],d['fastDAQ']['PMTtrig'],mode='full')


    t0 = d['PMTtraces']['t0_sec'][:, 0][-1] + d['PMTtraces']['t0_frac'][:, 0][-1] - np.argmax(q_fine) * dt - d['fastDAQ']['PMTtrig'][-1]

    output_dict['PMTt0_sec']=np.floor(t0)

    output_dict['PMTt0_frac']=t0-np.floor(t0)


    return output_dict








































