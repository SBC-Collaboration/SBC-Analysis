import numpy as np
import SBCcode as sbc
import os
import warnings


''' This module calculates goodness of fit for the global
    PICO calibration set

    The key function to be called (e.g. from emcee) is
    PICOcalChisq(theta).  Parameters defining the fit are
    all present in the introductory block that is executed
    when the module is imported.
'''


''' The following code is run when the library is imported -- it sets
    fenceposts in probability and thermodynamic threshold,
    it loads the MonteCarlo outputs, it defines nuisance parameters,
    and it defines the data to be fit.
'''
warnings.simplefilter("ignore", RuntimeWarning, 459, True)

single_ET = True

# First define fenceposts for Epts
p_fenceposts = np.array([0, .2, .5, .8, 1],
                        dtype=np.float64)
if not single_ET:
    # threshold_fenceposts = np.array([2, 3.2, 3.3, 6.1, 8.1, 13],
    # threshold_fenceposts = np.array([2, 3.2, 4.3, 6.1, 8.1, 13],
    threshold_fenceposts = np.array([1.8, 3.2, 4.3, 5.3],
                                    dtype=np.float64)
else:
    threshold_fenceposts = np.array([3.2], dtype=np.float64)
#    threshold_fenceposts = np.array([1.8], dtype=np.float64)

species_list = np.array([12, 19],
                        dtype=np.int32)
n_Epts = 0 #  p_fenceposts.size * threshold_fenceposts.size * species_list.size

pico2l_bestfit = [6.3, 4.3, 6.9, 4.9, 7.0, 5.9, 7.1, 6.2, 9.1, 6.9]

Epts_initialguess = np.array([[[4, 3],
                               [5.86, 4.72],
                               [7.6, 6.13],
                               [8.47, 6.93],
                               [10.66, 8.63],
                               [16, 14]],
                              [[5, 4],
                               [6.28, 5.08],
                               [8.16, 6.60],
                               [11.72, 8.74],
                               [15.2, 11.3],
                               [17, 14.1]],
                              [[6, 5],
                               [7.26, 6.05],
                               [9.43, 7.86],
                               [13.14, 10.00],
                               [16.8, 12.5],
                               [18, 15]],
                              [[6.1, 5.5],
                               [7.28, 6.52],
                               [9.54, 8.5],
                               [16.33, 10.89],
                               [18.9, 14],
                               [21, 16]],
                              [[8, 6.5],
                               [8.92, 7.59],
                               [10.73, 10.43],
                               [16.73, 13],
                               [19.8, 15],
                               [22, 17]]],
                             dtype=np.float64)

if single_ET:
    #Epts_initialguess = Epts_initialguess[:, 0, :]
    Epts_initialguess = pico2l_bestfit
    
# now define our experiment list
experiment_list = ['2013_97',
                   '2013_61',
                   # '2013_40',
                   '2014_97',
                   '2014_61',
                   '2014_50',
                   '2014_34',
                   'pico2l_2013_lt',
                   'pico2l_2013_ht',
                   # 'pico2l_2015',
                   # 'SbBe1',
                   # 'SbBe4',
                   # 'SbBe4_1inPb',
                   'SbBe4_2inPb',
                   #'pico60_Cf_run15',
                   ]

# Now find where stuff lives
topdir_searchlocations = ['/Users/cdahl/data/MCMCdata',
                          '/nashome/c/cdahl/MCMCdata',
                          '/home/mjin1989/Documents/LL/MCMCdata',
                          '/home/mjn693/Documents/LL/MCMCdata',
                          ]

for topdir in topdir_searchlocations:
    if os.path.isdir(topdir):
        break

# now load the simulation outputs
simfile_list = [os.path.join(topdir, exp, 'simout.bin')
                for exp in experiment_list]
neutron_sims = [sbc.read_bin(simfile) for simfile in simfile_list]
# neutron_sims is now a list of dictionaries, each dictionary
# with fields id(n), pos(n,3), Er(n), and species(n)


# now roll our hidden variable dice
np.random.seed()
n_hiddenvars = 11  # we can try this many times (-1) to make each bubble
for nsim in neutron_sims:
    nsim['hiddenvars'] = np.random.rand(nsim['Er'].size, n_hiddenvars)
# hiddenvars is shape (n, n_hiddenvars)

# now load the data tables (thresholds, counts, exposures, livetimes)
datafile_list = [os.path.join(topdir, exp, 'data.bin')
                 for exp in experiment_list]
neutron_data = [sbc.read_bin(datafile) for datafile in datafile_list]
# neutron_data is now a list of dictionaries, each dictionary
# with fields E_T(m), max_mult(m), exp(m), lt(m),
#             bkg_rate(m,k), counts(m,k), nuisance(m,3,b,k)

## creating masks to deal with problematic datasets
#for i in range(len(neutron_data)):
#    neutron_data[i]['mask'] = np.ones(neutron_data[i]['counts'].shape) 
    
## adding individual mask
## 2014_34
#neutron_data[5]['mask'][:,:]  = 0

## SbBe
#neutron_data[8]['mask'][:,0] = 0
    

if single_ET:
    for nd in neutron_data:
        # ETcut = (nd['E_T'] > 1.1) * (nd['E_T'] < 2.5)
        ETcut = (nd['E_T'] > 2.5) * (nd['E_T'] < 3.9)
        for k in nd.keys():
            nd[k] = nd[k][ETcut]
else:
    for nd in neutron_data:
        ETcut = (nd['E_T'] < 6)
        for k in nd.keys():
            nd[k] = nd[k][ETcut]
            
            
n_nuisance_all = [nd['nuisance'].shape[2] for nd in neutron_data]
if not np.all(np.diff(np.array(n_nuisance_all)) == 0):
    print('Inconsistent numbers of nusicance parameters between experiments')

n_nuisance = n_nuisance_all[0]
n_params = n_Epts + n_nuisance

# also need to set cuts -- start out with true fiducial cuts and
# false veto cuts, can overwrite cuts for individual experiments
# later as needed.  Cuts can be made on both pos and hidden variables
fidfun_list = [lambda pos, hv: np.ones(hv.shape, dtype=np.bool)
               for exp in experiment_list]

vetofun_list = [lambda pos, hv: np.zeros(hv.shape, dtype=np.bool)
                for exp in experiment_list]

# now overwrite specific fidfun's and vetofun's, using the ordering
# in experiment_list above
# PICO0.1 2013 having a foam between C3F8 & buffer,
# vetofun to kill any event with a bubble in foam.
vetofun_list[0] = lambda pos, hv: pos[:, 2] > 0.995  # 2013_97
vetofun_list[1] = lambda pos, hv: pos[:, 2] > 1.1438  # 2013_61

fidfun_list[6] = lambda pos, hv: hv < .7757  # pico2l_2013_lt
fidfun_list[7] = lambda pos, hv: hv < .6241  # pico2l_2013_ht


def PICOcalLL(theta, whichnuisance=np.ones(n_nuisance, dtype=np.bool)):
    ''' Top level log-likelihod function in this module

        This function inputs parameters and outputs a total
        log-likelihood.  Optional argument whichnuisance is a
        boolean that identifies which nuisance parameters are
        included in the input parameter list theta (the rest
        are held at zero)
        '''

#    Epts = np.reshape(theta[:n_Epts], (p_fenceposts.size,
#                                       threshold_fenceposts.size,
#                                       species_list.size))

    Epts = np.reshape(pico2l_bestfit, (p_fenceposts.size,
                                       threshold_fenceposts.size,
                                       species_list.size))
    
    if not CheckEpts(Epts):
        return -np.inf

    xpts = Epts / threshold_fenceposts[:, np.newaxis]

    eb = np.zeros(n_nuisance, dtype=np.float64)
    eb[whichnuisance] = theta[n_Epts:]

    total_chisq = np.sum(eb * eb)

    for i_exp in range(len(experiment_list)):
        
        if neutron_data[i_exp]['E_T'].size == 0:
            continue
        
        nu = SimulatedCounts(xpts, eb, i_exp)
        # nu.shape = neutron_Data[i_exp]['counts'].shape = (m, k)
        nu[nu < .1] = .1
        # that's to prevent bad -inf's from showing up...  really ought
        # to be handled in bkg_rate though.
        this_chisq = PoissonChisq(neutron_data[i_exp]['counts'], nu)
        #this_chisq = PoissonChisq(neutron_data[i_exp]['counts'], nu) * neutron_data[i_exp]['mask']
        bub_mult = np.cumsum(np.ones(nu.shape, dtype=np.int32), axis=1)
        bub_keep = bub_mult <= neutron_data[i_exp]['max_mult'][:, np.newaxis]
        total_chisq = total_chisq + np.sum(this_chisq[bub_keep].ravel())

    return -0.5 * total_chisq


def PICOcalLL_prior(theta, whichnuisance=np.ones(n_nuisance, dtype=np.bool)):
    ''' Prior log-likelihod function in this module

        This function inputs parameters and outputs a total
        log-likelihood.  Optional argument whichnuisance is a
        boolean that identifies which nuisance parameters are
        included in the input parameter list theta (the rest
        are held at zero)
        '''

#    Epts = np.reshape(theta[:n_Epts], (p_fenceposts.size,
#                                       threshold_fenceposts.size,
#                                       species_list.size))

#    if not CheckEpts(Epts):
#        return -np.inf
    
    eb = np.zeros(n_nuisance, dtype=np.float64)
    eb[whichnuisance] = theta[n_Epts:]

    total_chisq = np.sum(eb * eb)

    return -0.5 * total_chisq


def PICOcalLL_post(theta, whichnuisance=np.ones(n_nuisance, dtype=np.bool)):
    ''' Posterior log-likelihod function in this module

        This function inputs parameters and outputs a total
        log-likelihood.  Optional argument whichnuisance is a
        boolean that identifies which nuisance parameters are
        included in the input parameter list theta (the rest
        are held at zero)
        '''

#    Epts = np.reshape(theta[:n_Epts], (p_fenceposts.size,
#                                       threshold_fenceposts.size,
#                                       species_list.size))

    Epts = np.reshape(pico2l_bestfit, (p_fenceposts.size,
                                       threshold_fenceposts.size,
                                       species_list.size))
    
    if not CheckEpts(Epts):
        return -np.inf

    xpts = Epts / threshold_fenceposts[:, np.newaxis]

    eb = np.zeros(n_nuisance, dtype=np.float64)
    eb[whichnuisance] = theta[n_Epts:]

    total_chisq = 0

    for i_exp in range(len(experiment_list)):

        if neutron_data[i_exp]['E_T'].size == 0:
            continue
            
        nu = SimulatedCounts(xpts, eb, i_exp)
        # nu.shape = neutron_Data[i_exp]['counts'].shape = (m, k)
        nu[nu < .1] = .1
        # that's to prevent bad -inf's from showing up...  really ought
        # to be handled in bkg_rate though.
        this_chisq = PoissonChisq(neutron_data[i_exp]['counts'], nu)
        #this_chisq = PoissonChisq(neutron_data[i_exp]['counts'], nu) * neutron_data[i_exp]['mask']
        bub_mult = np.cumsum(np.ones(nu.shape, dtype=np.int32), axis=1)
        bub_keep = bub_mult <= neutron_data[i_exp]['max_mult'][:, np.newaxis]
        total_chisq = total_chisq + np.sum(this_chisq[bub_keep].ravel())

    return -0.5 * total_chisq


def PICOcalResiduals(theta, whichnuisance=np.ones(n_nuisance, dtype=np.bool)):
    ''' Residuals (signed log-likelihood) function in this module

        This function inputs parameters and outputs a total
        log-likelihood.  Optional argument whichnuisance is a
        boolean that identifies which nuisance parameters are
        included in the input parameter list theta (the rest
        are held at zero)
        '''

    Epts = np.reshape(theta[:n_Epts], (p_fenceposts.size,
                                       threshold_fenceposts.size,
                                       species_list.size))

    if not CheckEpts(Epts):
        return -np.inf

    xpts = Epts / threshold_fenceposts[:, np.newaxis]

    eb = np.zeros(n_nuisance, dtype=np.float64)
    eb[whichnuisance] = theta[n_Epts:]

    chilist = eb[whichnuisance]

    for i_exp in range(len(experiment_list)):
        nu = SimulatedCounts(xpts, eb, i_exp)
        # nu.shape = neutron_Data[i_exp]['counts'].shape = (m, k)
        bub_mult = np.cumsum(np.ones(nu.shape, dtype=np.int32), axis=1)
        bub_keep = bub_mult <= neutron_data[i_exp]['max_mult'][:, np.newaxis]
        new_chilist = PoissonChi(neutron_data[i_exp]['counts'], nu)
        chilist = np.append(chilist,
                            new_chilist[bub_keep])
    return chilist


def CheckEpts(Epts):
    ''' Checks if this set of Epts is allowed or not

        Checks that Epts is everywhere above thermodynaimc
        threshold, monotoninically increasing in p and t axes,
        and monotonically decreasing in s axis
        
        also adding a 3MeV hard cap (maximum carbon recoil for AmBe) 
        to prevent MCMC walker from getting lost
        '''

    if np.any((Epts < threshold_fenceposts[:, np.newaxis]).ravel()):
        return False

    if np.any(np.diff(Epts, axis=0).ravel() < 0):
        return False

    if np.any(np.diff(Epts, axis=1).ravel() < 0):
        return False

    if np.any(np.diff(Epts, axis=2).ravel() > 0):
        return False

    if np.any(Epts.ravel()) > 3000:
        return False
    
    return True


def SimulatedCounts(xpts, eb, i_exp):
    ''' calculate expected counts given simulated data and nuisances

        xpts:  (p,t,s) nd array of recoil energies
        eb:   (b,) nd array of low-level nuisance parameters
        i_exp:  scalar identifying which experiment to lookup

        output:  nu, size given by neutron_data[i_exp]['counts']
        '''

    exp_rescale = np.sum(neutron_data[i_exp]['nuisance'][:, 0, :, :] *
                         eb[:, np.newaxis], axis=1)
    lt_rescale = np.sum(neutron_data[i_exp]['nuisance'][:, 1, :, :] *
                        eb[:, np.newaxis], axis=1)
    thr_rescale = np.sum(neutron_data[i_exp]['nuisance'][:, 2, :, :] *
                         eb[:, np.newaxis], axis=1)
    # all the rescales have shape (m, k) where m is # of thresholds

    bkg_counts = neutron_data[i_exp]['bkg_rate'] *\
        neutron_data[i_exp]['lt'][:, np.newaxis] *\
        np.exp(lt_rescale)
    # bkg_counts.shape = (m, k)

    ET_eff = neutron_data[i_exp]['E_T'] * np.exp(thr_rescale[:, 0])
    # ET_eff.shape = (m,)

    p = EfficiencyInterpolation(neutron_sims[i_exp]['Er'],
                                neutron_sims[i_exp]['species'],
                                ET_eff, xpts)
    # p.shape = (m, n), where n is number of recoils in sim

    bubbles = neutron_sims[i_exp]['hiddenvars'][:, 1:] < p[:, :, np.newaxis]
    # bubbles.shape = (m, n, t), were t is number of trials per bubble

    fidcut = fidfun_list[i_exp](neutron_sims[i_exp]['pos'],
                                neutron_sims[i_exp]['hiddenvars'][:, 0])
    vetocut = vetofun_list[i_exp](neutron_sims[i_exp]['pos'],
                                  neutron_sims[i_exp]['hiddenvars'][:, 0])
    # cuts have shape (n,)

    # now, for each trial and threshold, we need to figure out
    # how many bubbles each simulated neutron creates, and
    # adjust for fiducial and veto cuts
    ev_midposts = np.nonzero(np.diff(neutron_sims[i_exp]['id']))[0]
    ev_posts = np.zeros(ev_midposts.size + 2, dtype=np.intp)
    ev_posts[1:-1] = ev_midposts + 1
    ev_posts[-1] = neutron_sims[i_exp]['id'].size

    bubcount = np.zeros((bubbles.shape[0],
                         bubbles.shape[1] + 1,
                         bubbles.shape[2]))
    np.cumsum(bubbles, axis=1, out=bubcount[:, 1:, :])
    nbub = np.diff(bubcount[:, ev_posts, :], axis=1)

    vetoed_bubbles = vetocut[:, np.newaxis] * bubbles
    vetocount = np.zeros((bubbles.shape[0],
                          bubbles.shape[1] + 1,
                          bubbles.shape[2]))
    np.cumsum(vetoed_bubbles, axis=1, out=vetocount[:, 1:, :])
    vetoed = np.diff(vetocount[:, ev_posts, :], axis=1) > 0
    nbub[vetoed] = 0

    fid_bubbles = fidcut[:, np.newaxis] * bubbles
    fidcount = np.zeros((bubbles.shape[0],
                         bubbles.shape[1] + 1,
                         bubbles.shape[2]))
    np.cumsum(fid_bubbles, axis=1, out=fidcount[:, 1:, :])
    infid = np.diff(fidcount[:, ev_posts, :], axis=1) > 0
    singles_out_of_fid = (nbub == 1) * (~infid)
    nbub[singles_out_of_fid] = 0
    # now we have nbub all set
    # with nbub.shape = (m, i, t) where i is number of neutrons

    #debugging when i_exp == 7, code fails as ET_eff = []
    
    #from IPython.core.debugger import Tracer; Tracer()() 
    
    # now we'll count events by multiplicity!
    nbub_all = nbub.reshape(ET_eff.size, -1)
    # nbub_all.shape = (m, t*i)

    nbub_test = np.array(range(bkg_counts.shape[1])) + 1
    # nbub_test.shape = (k,)

    sim_counts = np.sum(nbub_all[:, np.newaxis, :] ==
                        nbub_test[np.newaxis, :, np.newaxis],
                        axis=2)
    # nbub_hist.shape = (m, k)

    sim_counts[:, -1] += np.sum(nbub_all > bkg_counts.shape[1],
                                axis=1)

    trials_rescale = 1.0 / bubbles.shape[2]

    nu = bkg_counts + (sim_counts *
                       np.exp(exp_rescale) *
                       neutron_data[i_exp]['exp'][:, np.newaxis] *
                       trials_rescale)
    return nu


def EfficiencyInterpolation(E_r, s, E_T, xpts):
    ''' This calculates bubble nucleation efficiencies from xpts

        Inputs: E_r, 1-D ndarray of recoil energies, length n
                s,   1-D ndarray of recoil species, length n
                E_T, 1-D ndarray of detector thresholds, length m
                xpts, 3-D ndarray of E_r/E_T, shape is
                len(p_fenceposts), len(threshold_fenceposts), len(species_list)

        Outputs: 2-D ndarray of nucleation probabilities, shape m, n
        '''

    if threshold_fenceposts.size == 1:
        xps = xpts
    else:
        ET_ix = np.searchsorted(threshold_fenceposts, E_T)
        # ET_ix.shape = (m,)

        # Set extrapolation mode for E_T
        ET_ix[ET_ix < 1] = 1
        ET_ix[ET_ix >= threshold_fenceposts.size] =\
            threshold_fenceposts.size - 1

        # do this interp by hand so we get the index broadcasting right
        xps = xpts[:, ET_ix - 1, :] +\
            (xpts[:, ET_ix, :] - xpts[:, ET_ix - 1, :]) *\
            ((E_T[:, np.newaxis] -
              threshold_fenceposts[ET_ix - 1, np.newaxis]) /
             (threshold_fenceposts[ET_ix, np.newaxis] -
              threshold_fenceposts[ET_ix - 1, np.newaxis]))

    E_ps = xps * E_T[:, np.newaxis]
    # E_ps.shape = (len(p_fenceposts), m, len(species_list))

    # Now calculate probabilities for each species, threshold pair
    p = np.zeros((E_T.size, E_r.size))
    for i_ET in range(p.shape[0]):
        for i_s, this_s in enumerate(species_list):
            s_cut = (s == this_s)
            p[i_ET, s_cut] = np.interp(E_r[s_cut],
                                       E_ps[:, i_ET, i_s],
                                       p_fenceposts,
                                       left=0, right=1)
    return p


def PoissonChisq(n, nu):
    ''' Calculates Poisson chi-square '''
    chisq_0 = 2 * (nu - n)
    chisq_1 = -2 * n * np.log(nu / n)
    chisq_1[np.isnan(chisq_1)] = 0
    return chisq_0 + chisq_1


def PoissonChi(n, nu):
    ''' Calculates signed sqrt of Poisson chi-square '''
    return np.sqrt(PoissonChisq(n, nu)) * np.sign(nu - n)
