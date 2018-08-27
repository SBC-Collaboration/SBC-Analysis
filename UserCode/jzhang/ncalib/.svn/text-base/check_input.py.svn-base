import SBCcode.MonteCarlo.PICOcalGlobalLikelihood_v2 as pcgl
import numpy as np

# var_dict = vars(pcgl)
# for name in var_dict.keys():
#     # print(name, var_dict[name].__class__.__name__)
#     print(name, type(var_dict[name]))


print('-' * 80)
print('pcgl.p_fenceposts')
print(pcgl.p_fenceposts)

print('-' * 80)
print('pcgl.threshold_fenceposts')
print(pcgl.threshold_fenceposts)

print('-' * 80)
print('pcgl.single_ET')
print(pcgl.single_ET)


print('-' * 80)
print('pcgl.ETcut')
print(pcgl.ETcut)

print('-' * 80)
print('pcgl.Epts_initialguess', pcgl.Epts_initialguess.shape)
print(pcgl.Epts_initialguess)

# print(pcgl.counts)
# print(pcgl.i)
# print(pcgl.k)
print('-' * 80)
print('pcgl.max_recoils')
print(pcgl.max_recoils)
print('-' * 80)
print('pcgl.n_Epts')
print(pcgl.n_Epts)
print('-' * 80)
print('pcgl.n_hiddenvars')
print(pcgl.n_hiddenvars)
print('-' * 80)
print('pcgl.n_nuisance')
print(pcgl.n_nuisance)
print('-' * 80)
print('pcgl.n_params')
print(pcgl.n_params)

# print(nd <class 'collections.OrderedDict'>
# print(nsim <class 'collections.OrderedDict'>


# print(pcgl.unique)

print('-' * 80)
print('pcgl.species_list')
print(pcgl.species_list)
print('-' * 80)
print('pcgl.n_nuisance_all')
print(pcgl.n_nuisance_all)

print('-' * 80)
print('pcgl.topdir')
print(pcgl.topdir)

print('-' * 80)
print('pcgl.experiment_list')
[print(repr(x).rjust(20),i) for i, x in enumerate(pcgl.experiment_list)]
print('-' * 80)
print('pcgl.datafile_list')
[print(repr(x).rjust(20),i) for i, x in enumerate(pcgl.datafile_list)]
print('-' * 80)
print('pcgl.simfile_list')
[print(repr(x).rjust(20),i) for i, x in enumerate(pcgl.simfile_list)]


# print(pcgl.datafile_list)
# print(pcgl.simfile_list)
# print(pcgl.topdir_searchlocations)

# neutron_data
# data filed shapes
for i in range(len(pcgl.neutron_data)):
    print('-' * 80)
    print('neutron_data[{0:d}]'.format(i))
    print('-' * 80)
    for key in pcgl.neutron_data[i].keys():
        print(repr(key).rjust(15), repr(pcgl.neutron_data[i][key].dtype).ljust(17), pcgl.neutron_data[i][key].shape)

# data values
for i in range(len(pcgl.neutron_data)):
    print('-' * 80)
    print('neutron_data[{0:d}]'.format(i))
    print('-' * 80)
    for key in sorted(pcgl.neutron_data[i].keys()):
        if key not in ['nuisance']:
            if len(pcgl.neutron_data[i][key].shape) == 2:
                for j in range(pcgl.neutron_data[i][key].shape[1]):
                    print(repr(key).rjust(15), pcgl.neutron_data[i][key][:,j])
            else:
                print(repr(key).rjust(15), pcgl.neutron_data[i][key])
      #          'E_T' (1,)
      # 'sim_n_counts' (1,)
      #     'max_mult' (1,)
      #  'UdeM_counts' (1,)
      #          'exp' (1,)
      #           'lt' (1,)
      #     'bkg_rate' (1, 3)
      #       'counts' (1, 3)

# neutron_sims
# data shapes
for i in range(len(pcgl.neutron_sims)):
    print('-' * 80)
    print('neutron_sims[{0:d}]'.format(i))
    print('-' * 80)
    for key in sorted(pcgl.neutron_sims[i].keys()):
        print(repr(key).rjust(13), repr(pcgl.neutron_sims[i][key].dtype).ljust(17), pcgl.neutron_sims[i][key].shape)


# exp = (He3 counts in data)/(He3 counts in sim)
# np.set_printoptions(suppress=False)
np.set_printoptions(formatter={'all':lambda x: ' {0:.3e}'.format(x)})
print('-' * 80)
print('exp = (He3 counts in data)/(He3 counts in sim), summary table')
for i in range(len(pcgl.neutron_data)):
    print(repr(pcgl.experiment_list[i]).rjust(20), pcgl.neutron_data[i]['exp'])
np.set_printoptions()


# nuisance parameters
# first pass, collect unique matrix elements
print('-' * 80)
print('Independent errors')
nu_table = [];
for i_nu in range(pcgl.neutron_data[0]['nuisance'].shape[2]):
    tmp = []
    for i_exp in range(len(pcgl.neutron_data)):
        var = np.unique(pcgl.neutron_data[i_exp]['nuisance'][:,:,i_nu,:])
        non_zero = np.nonzero(var)
        if non_zero[0].size > 0:
            tmp.append(var[non_zero])
    # print(tmp)
    if len(tmp) > 0:
        nu_table.append(np.unique(np.concatenate(tmp)).tolist())
    else:
        nu_table.append([])

print(nu_table)


# second pass, find out the indices of the unique matrix elements
index_table = []
for i_nu in range(pcgl.neutron_data[0]['nuisance'].shape[2]):
    each_nu = []
    if len(nu_table[i_nu]) > 0:
        for i_var in range(len(nu_table[i_nu])):
            type_list = []
            data_list = []
            mult_list = []
            for i_exp in range(len(pcgl.neutron_data)):
                non_zero = np.nonzero( nu_table[i_nu][i_var] == pcgl.neutron_data[i_exp]['nuisance'][:, :, i_nu, :])
                # print(non_zero)
                if len(non_zero[0]) > 0:
                    type_list.append(non_zero[1])
                    if np.unique(non_zero[1]).size > 1:
                        Warning('Multiple type share the same error')
                    data_list.append(i_exp)
                    mult_list.append(non_zero[2])
                    if non_zero[2].size < pcgl.neutron_data[i_exp]['nuisance'].shape[3]:
                        Warning('Multiplicity mismatch')
            # print('-' * 80)
            # print(type_list)
            # print(data_list)
            # print(mult_list)
            # print('-' * 80)
            if len(type_list) > 0 and len(data_list) > 0 and len(mult_list) > 0:
                each_nu.append([np.unique(np.concatenate(type_list)).tolist(),
                               data_list,
                               np.unique(np.concatenate(mult_list)).tolist()])

    index_table.append(each_nu)

# [print(i) for i in index_table]

print('-' * 80)
print('Nuisance parameter table')
[print(x, '-->', y) for x, y in zip(nu_table, index_table)]
