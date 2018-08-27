import numpy as np


def ClosestApproach(hgb, lt):
    ib_cut = hgb['ibubimage'] < 2
    goodcut = hgb['nbubimage'][ib_cut] == 2
    b1_cut = (hgb['ibubimage'] == 1) * (hgb['nbubimage'] == 2)
    b2_cut = (hgb['ibubimage'] == 2) * (hgb['nbubimage'] == 2)

    n_ev = np.sum(ib_cut)

    output = dict(runid=hgb['runid'][ib_cut],
                  ev=hgb['ev'][ib_cut],
                  bubpos=np.nan+np.zeros((n_ev, 3), dtype=np.float64),
                  bubpos_adj=np.nan+np.zeros((n_ev, 3), dtype=np.float64),
                  bubX=np.nan+np.zeros((n_ev), dtype=np.float64),
                  bubY=np.nan+np.zeros((n_ev), dtype=np.float64),
                  bubZ=np.nan+np.zeros((n_ev), dtype=np.float64),
                  bubR2=np.nan+np.zeros((n_ev), dtype=np.float64),
                  bubD=np.nan+np.zeros((n_ev), dtype=np.float64))

    i1 = hgb['ipix'][b1_cut]
    i2 = hgb['ipix'][b2_cut]
    j1 = hgb['jpix'][b1_cut]
    j2 = hgb['jpix'][b2_cut]
    c1 = np.intp(hgb['cam'][b1_cut]) - 1
    c2 = np.intp(hgb['cam'][b2_cut]) - 1

    i1_low = np.intp(np.floor(i1))
    i2_low = np.intp(np.floor(i2))
    j1_low = np.intp(np.floor(j1))
    j2_low = np.intp(np.floor(j2))

    i1_rem = i1-i1_low
    i2_rem = i2-i2_low
    j1_rem = j1-j1_low
    j2_rem = j2-j2_low

    r1 = (i1_rem*j1_rem*np.transpose(lt[c1, i1_low, j1_low]) +
          (1-i1_rem)*j1_rem*np.transpose(lt[c1, i1_low-1, j1_low]) +
          i1_rem*(1-j1_rem)*np.transpose(lt[c1, i1_low, j1_low-1]) +
          (1-i1_rem)*(1-j1_rem)*np.transpose(lt[c1, i1_low-1, j1_low-1]))
    r2 = (i2_rem*j2_rem*np.transpose(lt[c2, i2_low, j2_low]) +
          (1-i2_rem)*j2_rem*np.transpose(lt[c2, i2_low-1, j2_low]) +
          i2_rem*(1-j2_rem)*np.transpose(lt[c2, i2_low, j2_low-1]) +
          (1-i2_rem)*(1-j2_rem)*np.transpose(lt[c2, i2_low-1, j2_low-1]))

    d1 = np.diff(r1, axis=0)[0]
    d2 = np.diff(r2, axis=0)[0]

    d1 = d1 / np.sqrt(np.sum(d1*d1, axis=0))
    d2 = d2 / np.sqrt(np.sum(d2*d2, axis=0))

    costheta = np.sum(d1*d2, axis=0)

    l1 = (np.sum((r1[0]-r2[0])*(costheta*d2 - d1), axis=0) /
          (1 - costheta*costheta))
    l2 = (np.sum((r2[0]-r1[0])*(costheta*d1 - d2), axis=0) /
          (1 - costheta*costheta))

    p1 = np.transpose(r1[0] + l1*d1)
    p2 = np.transpose(r2[0] + l2*d2)
    output['bubpos'][goodcut] = .5*(p1+p2)
    output['bubD'][goodcut] = np.sqrt(np.sum((p1-p2)*(p1-p2), axis=1))
    output['bubpos_adj'][:, 0] = output['bubpos'][:, 0] * .95 / .9
    output['bubpos_adj'][:, 1] = output['bubpos'][:, 1] + 0.05
    output['bubpos_adj'][:, 2] = ((output['bubpos'][:, 2] - 1.725) * 2.85 /
                                  (2.85 + 0.1*output['bubpos'][:, 1])) - 1.425
    output['bubX'] = output['bubpos_adj'][:, 0]
    output['bubY'] = output['bubpos_adj'][:, 1]
    output['bubZ'] = output['bubpos_adj'][:, 2]
    output['bubR2'] = (output['bubX']*output['bubX'] +
                       output['bubY']*output['bubY'])

    return output
