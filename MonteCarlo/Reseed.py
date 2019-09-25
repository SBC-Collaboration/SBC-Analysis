#author m.jin, 2018-03
#package for reseeding from iterative MCMC runs


import numpy as np


#compute new starting seeds table selected by WIMP sensitivity
def WIMP(chain,lnprobability,E_band, bin_number, F_table,C_table,starting_energy,interval): 
    

    epoch_starting_points = np.zeros(24)    
    
    
    table = np.zeros((E_band.shape[0],chain.shape[1]*chain.shape[2],3))

    for i_band in range(E_band.shape[0]):
        for i_walker in range(chain.shape[1]):
            for i_step in range(chain.shape[2]):
        
                SI_wimp = 0
                SD_wimp = 0
            
                F = np.cumsum([np.exp(chain[0,i_walker,i_step,1:10:2])])+3.2
                C = np.cumsum([np.exp(chain[0,i_walker,i_step,0:10:2])])+3.2
                F_interp = np.zeros(171)
                C_interp = np.zeros(171)
        
                for i_interp in range(171):
                    F_interp[i_interp] = np.interp(starting_energy + i_interp*interval ,F,[0, .2, .5, .8, 1.0])
                    C_interp[i_interp] = np.interp(starting_energy + i_interp*interval ,C,[0, .2, .5, .8, 1.0])
                
                    SD_wimp += .778/19/19 * F_interp[i_interp]*F_table['table'][i_band,i_interp] 
                    SI_wimp += SI_F_numerator/SI_denominator * F_interp[i_interp]*F_table['table'][i_band,i_interp] + SI_C_numerator/SI_denominator * C_interp[i_interp]*C_table['table'][i_band,i_interp]
            
                table[i_band, i_walker*chain.shape[2] + i_step,0] = SI_wimp
                table[i_band, i_walker*chain.shape[2] + i_step,1] = SD_wimp
                table[i_band, i_walker*chain.shape[2] + i_step,2] = lnprobability[0, i_walker , i_step ]
        

        
#for i_band in range(E_band.shape[0]):
        for i_species in range(2):
            for i_bin in range(bin_number):
                
                b = np.max(table[i_band,:,i_species])
                a = np.min(table[i_band,:,i_species])
            
                bin_size = (b-a)/bin_number
            
                index = np.asarray(np.where((table[i_band,:,0] >= a + i_bin*bin_size) & (table[i_band,:,0] <a + (i_bin+1)*bin_size)))
                unique, unique_indices = np.unique(table[i_band,index,2],return_index=True)
            
                if unique.size != 0: 
                    
                    #index[unique_indices[-1]]//chain.shape[2]
                    #index[unique_indices[-1]]%chain.shape[2]
                    
                    epoch_starting_points = np.vstack((epoch_starting_points,chain[0,index[0,unique_indices[-1]]//chain.shape[2],index[0,unique_indices[-1]]%chain.shape[2]])) 
            
    ## seeding for next epoch
    #for i_walker in range(nwalkers):                  
    #    epoch_starting_points[0,i_walker,:] = chain[0,unique_indices[-i_walker]//epoch_nstep,unique_indices[-i_walker]%epoch_nstep,:] 
    
    epoch_starting_points = np.delete(epoch_starting_points,0,axis=0)
    
    if epoch_starting_points.shape[0]%2 == 1:
        epoch_starting_points = np.insert(epoch_starting_points,0, epoch_starting_points[0,:],axis = 0)
        
    epoch_starting_points = np.expand_dims(epoch_starting_points,axis=0)
        
    
    return epoch_starting_points;


## sdfs

def get_table(E_band,test_chain, test_probability):

    SI_denominator = 3*12.011 + 8*18.998403163
    SI_C_numerator = 3*12.011
    SI_F_numerator = 8*18.998403163

    #E_band = np.array([5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,400,1600,5600])

    #E_band = np.array([20])

    table = np.zeros((E_band.shape[0],test_chain.shape[1]*test_chain.shape[2],3))

    for i_band in range(E_band.shape[0]):
        for i_walker in range(test_chain.shape[1]):
            for i_step in range(test_chain.shape[2]):
        
                SI_wimp = 0
                SD_wimp = 0
            
                F = np.cumsum([np.exp(test_chain[0,i_walker,i_step,1:10:2])])+3.2
                C = np.cumsum([np.exp(test_chain[0,i_walker,i_step,0:10:2])])+3.2
                F_interp = np.zeros(171)
                C_interp = np.zeros(171)

                for i_interp in range(171):
                    F_interp[i_interp] = np.interp(3.0 + i_interp*.1 ,F,[0, .2, .5, .8, 1.0])
                    C_interp[i_interp] = np.interp(3.0 + i_interp*.1 ,C,[0, .2, .5, .8, 1.0])
                
                    SD_wimp += .778/19/19 * F_interp[i_interp]*F_table['table'][i_band,i_interp] 
                    SI_wimp += SI_F_numerator/SI_denominator * F_interp[i_interp]*F_table['table'][i_band,i_interp] + SI_C_numerator/SI_denominator * C_interp[i_interp]*C_table['table'][i_band,i_interp]
            
                SI_wimp += SI_F_numerator/SI_denominator * integral_table['integral_table'][i_band,1] + SI_C_numerator/SI_denominator * integral_table['integral_table'][i_band,0] 
                SD_wimp += .778/19/19 * integral_table['integral_table'][i_band,1]
            
                table[i_band, i_walker*test_chain.shape[2] + i_step,0] = SI_wimp
                table[i_band, i_walker*test_chain.shape[2] + i_step,1] = SD_wimp
                table[i_band, i_walker*test_chain.shape[2] + i_step,2] = test_lnprobability[0, i_walker , i_step ]
                
                
                
    return table
        



