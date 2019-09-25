import numpy as np
import SBCcode as sbc

# defining some constants
SI_denominator = 3*12.011 + 8*18.998403163
SI_C_numerator = 3*12.011
SI_F_numerator = 8*18.998403163

default_num_mass = 47 

## taking desired mass_list, MCMC chain & lnprobability, WIMPspectra as input
## Output a table of SI_wimp and SD_wimp sensitivity versus lnprobability
def get_table(mass_array,test_chain, test_lnprobability, WIMPspectra,bin_length,bin_width, threshold = 3.2, which_mass=np.ones(default_num_mass, dtype=np.bool)):

    table = np.zeros((np.sum(which_mass),test_chain.shape[1]*test_chain.shape[2],3))
    
    for i_band in range(np.sum(which_mass)):
        for i_walker in range(test_chain.shape[1]):
            for i_step in range(test_chain.shape[2]):

                SI_wimp = 0
                SD_wimp = 0

                F = (np.cumsum([np.exp(test_chain[0,i_walker,i_step,1:10:2])])+3.2 ) * threshold/3.2
                C = (np.cumsum([np.exp(test_chain[0,i_walker,i_step,0:10:2])])+3.2 ) * threshold/3.2
                F_interp = np.zeros(bin_length)
                C_interp = np.zeros(bin_length)

                F_interp = np.interp(1.0 + np.arange(0,bin_length*bin_width,bin_width) ,F,[0, .2, .5, .8, 1.0])
                C_interp = np.interp(1.0 + np.arange(0,bin_length*bin_width,bin_width) ,C,[0, .2, .5, .8, 1.0])

                SD_wimp += np.dot(F_interp, WIMPspectra['SD_F_table'][which_mass,:][i_band])* bin_width 
                SI_wimp += (SI_F_numerator/SI_denominator * np.dot(F_interp, WIMPspectra['SI_F_table'][which_mass,:][i_band]) + SI_C_numerator/SI_denominator * np.dot(C_interp, WIMPspectra['SI_C_table'][which_mass,:][i_band])) * bin_width

                table[i_band, i_walker*test_chain.shape[2] + i_step,0] = SI_wimp
                table[i_band, i_walker*test_chain.shape[2] + i_step,1] = SD_wimp
                table[i_band, i_walker*test_chain.shape[2] + i_step,2] = test_lnprobability[0, i_walker , i_step ]


    #pickle.dump(table, open( "/home/mjn693/Documents/LL/Python_objects/reparam/production/WIMP/table%d.p"%i_epoch, "wb" ))



                
                
    return table


def get_table_multi(mass_array,test_chain, test_lnprobability, WIMPspectra,bin_length,bin_width, reparam_threshold = np.array([0,0]),threshold = np.array([2.45,3.2]), chisq_hard_cap=65, target_threshold = np.array([2.45,3.29]), which_mass=np.ones(default_num_mass, dtype=np.bool)):

    table = np.zeros((np.sum(which_mass),test_chain.shape[1]*test_chain.shape[2],3, threshold.shape[0]))
    
    for i_band in range(np.sum(which_mass)):
        for i_walker in range(test_chain.shape[1]):
            for i_step in range(test_chain.shape[2]):
                for i_threshold in range(threshold.shape[0]):


                    SI_wimp = 0
                    SD_wimp = 0

                    F = (np.cumsum([np.exp(test_chain[0,i_walker,i_step,1+2*i_threshold:10*threshold.shape[0]:2*threshold.shape[0]])])+reparam_threshold[i_threshold] ) * target_threshold[i_threshold]/threshold[i_threshold]
                    C = (np.cumsum([np.exp(test_chain[0,i_walker,i_step,0+2*i_threshold:10*threshold.shape[0]:2*threshold.shape[0]])])+reparam_threshold[i_threshold] ) * target_threshold[i_threshold]/threshold[i_threshold]
                    F_interp = np.zeros(bin_length)
                    C_interp = np.zeros(bin_length)

                    #from IPython.core.debugger import Tracer; Tracer()() 
                    
                    F_interp = np.interp(1.0 + np.arange(0,bin_length*bin_width,bin_width) ,F,[0, .2, .5, .8, 1.0])
                    C_interp = np.interp(1.0 + np.arange(0,bin_length*bin_width,bin_width) ,C,[0, .2, .5, .8, 1.0])

                    SD_wimp += np.dot(F_interp, WIMPspectra['SD_F_table'][which_mass,:][i_band])* bin_width 
                    SI_wimp += (SI_F_numerator/SI_denominator * np.dot(F_interp, WIMPspectra['SI_F_table'][which_mass,:][i_band]) + SI_C_numerator/SI_denominator * np.dot(C_interp, WIMPspectra['SI_C_table'][which_mass,:][i_band])) * bin_width

                    table[i_band, i_walker*test_chain.shape[2] + i_step,0,i_threshold] = SI_wimp
                    table[i_band, i_walker*test_chain.shape[2] + i_step,1,i_threshold] = SD_wimp
                    table[i_band, i_walker*test_chain.shape[2] + i_step,2,i_threshold] = test_lnprobability[0, i_walker , i_step ]             
                
    return table

## checking knobs monotonicity
##
def get_table_multi_mono(mass_array,test_chain, test_lnprobability, WIMPspectra,bin_length,bin_width, threshold = np.array([2.45,3.2]), chisq_hard_cap=65, target_threshold = np.array([2.45,3.29]), which_mass=np.ones(default_num_mass, dtype=np.bool)):

    table = np.zeros((np.sum(which_mass),test_chain.shape[1]*test_chain.shape[2],3, threshold.shape[0]))
    
    for i_band in range(np.sum(which_mass)):
        for i_walker in range(test_chain.shape[1]):
            for i_step in range(test_chain.shape[2]):
                
                F0  =(np.cumsum([np.exp(test_chain[0,i_walker,i_step,1+2*0:10*threshold.shape[0]:2*threshold.shape[0]])])+threshold[0] ) * target_threshold[0]/threshold[0]
                
                F1 = (np.cumsum([np.exp(test_chain[0,i_walker,i_step,1+2*1:10*threshold.shape[0]:2*threshold.shape[0]])])+threshold[1] ) * target_threshold[1]/threshold[1]
                
                if any(F0>F1):
                    continue
                
                for i_threshold in range(threshold.shape[0]):


                    SI_wimp = 0
                    SD_wimp = 0

                    F = (np.cumsum([np.exp(test_chain[0,i_walker,i_step,1+2*i_threshold:10*threshold.shape[0]:2*threshold.shape[0]])])+threshold[i_threshold] ) * target_threshold[i_threshold]/threshold[i_threshold]
                    C = (np.cumsum([np.exp(test_chain[0,i_walker,i_step,0+2*i_threshold:10*threshold.shape[0]:2*threshold.shape[0]])])+threshold[i_threshold] ) * target_threshold[i_threshold]/threshold[i_threshold]
                    F_interp = np.zeros(bin_length)
                    C_interp = np.zeros(bin_length)

                    #from IPython.core.debugger import Tracer; Tracer()() 
                    
                    F_interp = np.interp(1.0 + np.arange(0,bin_length*bin_width,bin_width) ,F,[0, .2, .5, .8, 1.0])
                    C_interp = np.interp(1.0 + np.arange(0,bin_length*bin_width,bin_width) ,C,[0, .2, .5, .8, 1.0])

                    SD_wimp += np.dot(F_interp, WIMPspectra['SD_F_table'][which_mass,:][i_band])* bin_width 
                    SI_wimp += (SI_F_numerator/SI_denominator * np.dot(F_interp, WIMPspectra['SI_F_table'][which_mass,:][i_band]) + SI_C_numerator/SI_denominator * np.dot(C_interp, WIMPspectra['SI_C_table'][which_mass,:][i_band])) * bin_width

                    table[i_band, i_walker*test_chain.shape[2] + i_step,0,i_threshold] = SI_wimp
                    table[i_band, i_walker*test_chain.shape[2] + i_step,1,i_threshold] = SD_wimp
                    table[i_band, i_walker*test_chain.shape[2] + i_step,2,i_threshold] = test_lnprobability[0, i_walker , i_step ]             
                
    return table



def get_epoch_starting_points(mass_array, bin_number, table, test_chain,which_mass=np.ones(default_num_mass, dtype=np.bool)):

    epoch_starting_points = np.zeros(24) 

    for i_band in range(np.sum(which_mass)):
        for i_species in range(2):
            for i_bin in range(bin_number):

                b = np.max(table[i_band,:,i_species])
                a = np.min(table[i_band,:,i_species])

                bin_size = (b-a)/bin_number
                index = np.asarray(np.where((table[i_band,:,0] >= a + i_bin*bin_size) & (table[i_band,:,0] <a + (i_bin+1)*bin_size)))
                unique, unique_indices = np.unique(table[i_band,index,2],return_index=True)

                if unique.size != 0: 
                    epoch_starting_points = np.vstack((epoch_starting_points,test_chain[0,index[0,unique_indices[-1]]//test_chain.shape[2],index[0,unique_indices[-1]]%test_chain.shape[2]])) 

    epoch_starting_points = np.delete(epoch_starting_points,0,axis=0)

    if epoch_starting_points.shape[0]%2 == 1:
        epoch_starting_points = np.insert(epoch_starting_points,0, epoch_starting_points[0,:],axis = 0)

    epoch_starting_points = np.expand_dims(epoch_starting_points,axis=0)
    
    return epoch_starting_points 

def get_WIMP_LF(mass_array, bin_number, table,which_mass=np.ones(default_num_mass, dtype=np.bool)):

    WIMP_LF = np.zeros((np.sum(which_mass), bin_number,2,2)) 

    for i_band in range(np.sum(which_mass)):
        for i_bin in range(bin_number):
            for i_species in range(2):

                b = np.max(table[i_band,:,i_species])
                a = np.min(table[i_band,:,i_species])

                bin_size = (b-a)/bin_number
                index = np.asarray(np.where((table[i_band,:,i_species] >= a + i_bin*bin_size) & (table[i_band,:,i_species] <a + (i_bin+1)*bin_size)))
                unique, unique_indices = np.unique(table[i_band,index,2],return_index=True)

                WIMP_LF[i_band,i_bin,i_species,0] = table[i_band, index[0,unique_indices[-1]],i_species]
                WIMP_LF[i_band,i_bin,i_species,1] = table[i_band, index[0,unique_indices[-1]],2]
    
    return WIMP_LF 

def get_WIMP_LF_multi(mass_array, bin_number, table,threshold = [2.45, 3.29],chisq_hard_cap=65, which_mass=np.ones(default_num_mass, dtype=np.bool)):

    WIMP_LF = np.zeros((np.sum(which_mass), bin_number,2,2, len(threshold))) 

    for i_band in range(np.sum(which_mass)):
        for i_bin in range(bin_number):
            for i_species in range(2):
                for i_threshold in range(len(threshold)):
                    
                    clean_index = np.asarray(np.where((table[i_band,:,2,i_threshold]*-2 <= chisq_hard_cap)))
                    
                    b = np.max(table[i_band,clean_index,i_species,i_threshold])
                    a = np.min(table[i_band,clean_index,i_species,i_threshold])

                    bin_size = (b-a)/bin_number
                    index = np.asarray(np.where((table[i_band,:,i_species,i_threshold] >= a + i_bin*bin_size) & (table[i_band,:,i_species,i_threshold] <a + (i_bin+1)*bin_size)))
                    unique, unique_indices = np.unique(table[i_band,index,2,i_threshold],return_index=True)
                    
                    if unique.size==0:
                        continue
                    
                    WIMP_LF[i_band,i_bin,i_species,0,i_threshold] = table[i_band, index[0,unique_indices[-1]],i_species,i_threshold]
                    WIMP_LF[i_band,i_bin,i_species,1,i_threshold] = table[i_band, index[0,unique_indices[-1]],2,i_threshold]
    
    return WIMP_LF 

def get_WIMP_LF_multi_special(mass_array, bin_number, table, use_mass_index = 2, threshold = [2.45, 3.29],chisq_hard_cap=65, which_mass=np.ones(default_num_mass, dtype=np.bool)):

    WIMP_LF = np.zeros((np.sum(which_mass), bin_number,2,2, len(threshold))) 

    for i_band in range(np.sum(which_mass)):
        for i_bin in range(bin_number):
            for i_species in range(2):
                for i_threshold in range(len(threshold)):
                    
                    
                    
                    clean_index = np.asarray(np.where((table[use_mass_index,:,2,i_threshold]*-2 <= chisq_hard_cap)))
                    
                    b = np.max(table[use_mass_index,clean_index,i_species,i_threshold])
                    a = np.min(table[use_mass_index,clean_index,i_species,i_threshold])

                    bin_size = (b-a)/bin_number
                    index = np.asarray(np.where((table[use_mass_index,:,i_species,i_threshold] >= a + i_bin*bin_size) & (table[use_mass_index,:,i_species,i_threshold] <a + (i_bin+1)*bin_size)))
                    unique, unique_indices = np.unique(table[use_mass_index,index,2,i_threshold],return_index=True)
                    
                    if unique.size==0:
                        continue
                    
                    WIMP_LF[i_band,i_bin,i_species,0,i_threshold] = table[i_band, index[0,unique_indices[-1]],i_species,i_threshold]
                    WIMP_LF[i_band,i_bin,i_species,1,i_threshold] = table[i_band, index[0,unique_indices[-1]],2,i_threshold]
    
    return WIMP_LF 

def get_WIMP_LF_multi_vertical(mass_array, bin_number, table,threshold = [2.45, 3.29],chisq_hard_cap=65, which_mass=np.ones(default_num_mass, dtype=np.bool)):

    WIMP_LF = np.zeros((np.sum(which_mass), bin_number*2,2,2, len(threshold))) 

    for i_band in range(np.sum(which_mass)):
        for i_species in range(2):
            for i_threshold in range(len(threshold)):
                b = chisq_hard_cap
                a = np.min(table[i_band,:,2,i_threshold]*-2)

                for i_bin in range(bin_number):    

                    bin_size = (b-a)/bin_number
                    index = np.asarray(np.where((table[i_band,:,2,i_threshold]*-2 >= a + i_bin*bin_size) & (table[i_band,:,2,i_threshold]*-2 <a + (i_bin+1)*bin_size)))
                    unique, unique_indices = np.unique(table[i_band,index,i_species,i_threshold],return_index=True)
                    WIMP_LF[i_band,i_bin,i_species,0,i_threshold] = table[i_band, index[0,unique_indices[-1]],i_species,i_threshold]
                    WIMP_LF[i_band,i_bin,i_species,1,i_threshold] = table[i_band, index[0,unique_indices[-1]],2,i_threshold]

                    WIMP_LF[i_band,bin_number*2 -1 - i_bin,i_species,0,i_threshold] = table[i_band, index[0,unique_indices[0]],i_species,i_threshold]
                    WIMP_LF[i_band,bin_number*2 -1 - i_bin,i_species,1,i_threshold] = table[i_band, index[0,unique_indices[0]],2,i_threshold]
    return WIMP_LF 


def get_log_bin(x_min,x_max, n_bins):
    
    if x_min==0:
        x_min = x_max/1000000
    
    bin_size = np.exp(np.log(x_max/x_min) / n_bins)

    return bin_size 


def get_nodes(test_chain, thresholds = [2.45, 3.2]):
    
    nodes = np.zeros((test_chain.shape[0],test_chain.shape[1],test_chain.shape[2],test_chain.shape[3]))
    
    for i_walker in range(test_chain.shape[1]):
        for i_step in range(test_chain.shape[2]):
            
            nodes[0,i_walker,i_step,20:] = test_chain[0,i_walker,i_step,20:]
            
            nodes[0,i_walker,i_step,0:20:4] = np.cumsum(np.exp(test_chain[0,i_walker,i_step,0:20:4])) + thresholds[0]
            nodes[0,i_walker,i_step,1:20:4] = np.cumsum(np.exp(test_chain[0,i_walker,i_step,1:20:4])) + thresholds[0]
            nodes[0,i_walker,i_step,2:20:4] = np.cumsum(np.exp(test_chain[0,i_walker,i_step,2:20:4])) + thresholds[1]
            nodes[0,i_walker,i_step,3:20:4] = np.cumsum(np.exp(test_chain[0,i_walker,i_step,3:20:4])) + thresholds[1]
            
            

    return nodes 

