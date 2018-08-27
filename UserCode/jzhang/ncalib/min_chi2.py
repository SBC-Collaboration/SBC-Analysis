import SBCcode.MonteCarlo.PICOcalGlobalLikelihood_v2_jz as pcgl
import numpy as np
from iminuit import Minuit, describe, Struct
from iminuit.util import make_func_code
# import sys
import ipdb


# Minuit target functions can only take scalar arguments, so fake the function signature
class Chi2Functor:
    def __init__(self, f, n_params):
        self.f = f

        varnames = []
        for i in range(n_params):
            varnames.append('p{:d}'.format(i))

        self.func_code = make_func_code(varnames)  # fake signature
        self.func_defaults = None  # this keeps np.vectorize happy

    def __call__(self, *args):
        # pass args to the original function
        chi2 = -2 * self.f(np.array([*args], dtype=np.float64))
        return chi2


# # Fake function signature
# def chi2_0(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
#            p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23,
#            p24):
#     testparams = np.array([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
#                            p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24])
#     return -2 * pcgl.PICOcalLL(testparams)


# v1 theta to v2 input for Minuit
def v2_input_v1(v1_theta):
    var = v1_theta.copy()
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

    pa = dict(('p{:d}'.format(i), x) for i, x in enumerate(var))
    pa.update(dict(('error_p{:d}'.format(i), 1.0) for i, x in enumerate(var)))
    pa.update({'errordef': 1.0})
    return pa


# # v2 theta to v1 input for Minuit
# def v1_input_v2(v2_theta):
#     var = v2_theta.copy()
#     dEpts = np.reshape(var[:pcgl.n_Epts], (pcgl.p_fenceposts.size,
#                                           pcgl.threshold_fenceposts.size,
#                                           pcgl.species_list.size))
#
#     Epts = np.cumsum(dEpts, axis=0)
#
#     var[:pcgl.n_Epts] = np.reshape(Epts, (pcgl.p_fenceposts.size *
#                                            pcgl.threshold_fenceposts.size *
#                                            pcgl.species_list.size))
#
#     pa = dict(('p{:d}'.format(i + 1), x) for i, x in enumerate(var))
#     pa.update(dict(('error_p{:d}'.format(i + 1), 1.0) for i, x in enumerate(var)))
#     pa.update({'errordef': 1.0})
#     return pa

# v2 theta to v2 input for Minuit
def v2_input_v2(v2_theta):
    var = v2_theta.copy()

    pa = dict(('p{:d}'.format(i), x) for i, x in enumerate(var))
    pa.update(dict(('error_p{:d}'.format(i), 1.0) for i, x in enumerate(var)))
    pa.update({'errordef': 1.0})
    return pa


def v2_theta_v1(v1_theta):
    var = v1_theta.copy()
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


def v1_theta_v2(v2_theta):
    var = v2_theta.copy()
    dEpts = np.reshape(var[:pcgl.n_Epts], (pcgl.p_fenceposts.size,
                                           pcgl.threshold_fenceposts.size,
                                           pcgl.species_list.size))
    Epts = np.cumsum(dEpts, axis=0)
    var[:pcgl.n_Epts] = np.reshape(Epts, (pcgl.p_fenceposts.size *
                                          pcgl.threshold_fenceposts.size *
                                          pcgl.species_list.size))
    return var


chi2 = Chi2Functor(pcgl.PICOcalLL, pcgl.n_params)

v1_theta = np.array(
    [8.1508, 3.8139, 9.3993, 3.8495, 12.5022, 5.1765, 14.4248, 6.4186, 22.0866, 12.5607, -0.2518, -2.3705, 2.1569,
     -0.6839, -1.9952, -0.4400, 1.7738, -3.5711, 1.2013, -1.1032, 1.3775, -1.8020, -0.1159, 2.3455])

v1_theta = np.array(
    [5.5238, 3.5621, 8.1342, 3.7994, 9.7790, 4.5712, 13.6781, 5.1958, 17.7044, 12.5956, 0.7342, -0.2849, -0.4499,
     -0.5762, 0.0693, -1.2320, 0.4503, 0.3148, 1.6189, -0.0172, -0.3546, 0.8604, 0.0126, 0.1009])

lnp = pcgl.PICOcalLL(v2_theta_v1(v1_theta))
print(lnp)

a0 = chi2(8.130482516771455, 3.8138960973263827, 1.2485265107853758, 0.03559997115378979, 3.1030182471540675,
          1.3270049596236873, 1.922785910687116, 1.0765548611236906, 7.663286150748873, 6.142091665490981,
          -0.2518002325608806, -2.37049679567336, -0.011037064034669112, -0.6838756106641036, -1.9951985398411285,
          -0.18733848097600925, 0.33817512620994966, -1.033389116224571, 1.6272382073109755, -0.0033786023815845453,
          0.004218659268035917, -0.00551871079563987, -0.022818011602259484, 0.6425853568603498)
print(a0)

a1 = chi2(8.150800, 3.813900, 1.248500, 0.035600, 3.102900, 1.327000, 1.922600, 1.242100, 7.661800, 6.142100, -0.251800,
          -2.370500, 2.156900, -0.683900, -1.995200, -0.440000, 1.773800, -3.571100, 1.201300, -1.103200, 1.377500,
          -1.802000, -0.115900, 2.345500)
print(a1)

b2 = (5.522868138136657, 3.558935274360199, 2.610396143923082, 0.23730023172241915, 1.6444299668577755,
      0.7976446271170996, 3.8997258524425864, 0.6944576278195789, 4.027183598141905, 7.374498997417081,
      -0.32207681121177373, -0.12409145077800868, -0.38747651006546263, -0.026483786144845867, 0.06948497661324038,
      0.19157039640893028, 0.3347527553076381, -0.5505222424384045, 1.0935073414668803, -0.0003084409295731206,
      -0.006358904375787696, 0.015429219904940047, 0.00033782192358825723, 0.6266110967879278)
a2 = chi2(*b2)

print(a2)

pp = {'p0': 3.8794122437980905, 'p1': 3.457310665131314, 'p2': 2.9806499862726366, 'p3': 0.05668130791705585, 'p4': 2.06264457851852, 'p5': 1.0324908638121335, 'p6': 1.7210021058702543, 'p7': 0.7213258859204852, 'p8': 4.215707642903244, 'p9': 7.39447946369445, 'p10': 0.0510140468908305, 'p11': 0.07057098056657085, 'p12': -0.31428499999693665, 'p13': -0.11584364946523429, 'p14': -0.22103318512268863, 'p15': -0.11075865636106251, 'p16': 0.3387391651328604, 'p17': -0.7342477056197758, 'p18': 1.0221583648291954, 'p19': -1.3874830921855485e-09, 'p20': -2.7995334575865994e-08, 'p21': 6.458292411954612e-08, 'p22': -0.031681578270365905, 'p23': 0.6360259968736162};
vv = [pp[key] for key in pp.keys()]

print(chi2(*vv))

#
# b0 = v2_input_v2(np.array(b2))
# m = Minuit(chi2, **b0)
# m.migrad()
# print(m.values)
# print(m.errors)
# print(m)
# # ipdb.set_trace()
