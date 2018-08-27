import SBCcode.MonteCarlo.PICOcalGlobalLikelihood_v2 as pcgl
import numpy as np
import pickle
import os


def v2_theta(theta):
    var = theta.copy()
    Epts = np.reshape(var[:pcgl.n_Epts], (pcgl.p_fenceposts.size,
                                          pcgl.threshold_fenceposts.size,
                                          pcgl.species_list.size))

    dEpts = np.zeros(Epts.shape, dtype=Epts.dtype)
    # ipdb.set_trace()
    dEpts[1:, :, :] = np.diff(Epts, 1, axis=0)
    dEpts[0, :, :] = Epts[0, :, :]

    var[:pcgl.n_Epts] = np.reshape(dEpts, (pcgl.p_fenceposts.size *
                                           pcgl.threshold_fenceposts.size *
                                           pcgl.species_list.size))

    return var


mcmc_dir = '/home/mjn693/Documents/LL/Chain/3_2_single_Russell/'

mcmc_files = [
    'Single_truncated_lnlikelihood.p',
    'Single_truncated_chain.p',
    'Single_truncated_lnprobability.p']

mcmc_loglikelihood = pickle.load(open(os.path.join(mcmc_dir, mcmc_files[0]), 'rb'))
mcmc_chain = pickle.load(open(os.path.join(mcmc_dir, mcmc_files[1]), 'rb'))
mcmc_logpost = pickle.load(open(os.path.join(mcmc_dir, mcmc_files[2]), 'rb'))

n_walker, n_step, n_var = mcmc_chain.shape

chi2 = np.zeros((n_walker, n_step))

n_walker = 5
n_step = 5

for i_walker in range(n_walker):
    for i_step in range(n_step):
        theta = mcmc_chain[i_walker, i_step, :].ravel()
        chi2[i_walker, i_step] = pcgl.PICOcalLL(v2_theta(theta))

ratio = chi2 / mcmc_logpost
print(ratio[:n_walker, :n_step])
# ipdb.set_trace()
