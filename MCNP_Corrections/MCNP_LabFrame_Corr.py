# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 16:07:15 2018

@author: theodorebaker
"""

import numpy as np
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Lorentz Boost with Momentum and Theta inputs -- all in keV
def boost(p_photon, theta, mass_target = 9.012183078): # reaction-dependent
    amu = 931494; # keV/c^2 per amu, per physics.bu.edu
    mass_neutron = 1008664.91585; # listed in "keV/c^2" 1008664.91585
    Q = 1661; # Be9_mass - Be8_mass; # listed in keV # reaction-dependent
    
    # Masses
    M = mass_target*amu;
    m_1 = mass_neutron;
    m_2 = M - m_1 + Q;

    # Rotational Lorentz "phi"
    phi = .5*np.log((2*p_photon/M) + 1);
    
    # Lorentz Boost Matrix  
    dim = 4
    Lambda = [[0 for x in range(dim)] for y in range(dim)];
    Lambda[0][0] = np.cosh(phi);
    Lambda[0][1] = -np.sinh(phi);
    Lambda[1][0] = -np.sinh(phi);
    Lambda[1][1] = np.cosh(phi);
    Lambda[2][2] = 1;
    Lambda[3][3] = 1;
    
    p_4lab = np.array([[p_photon],[p_photon],[0],[0]]);  
    
    p_4cm_photon = np.matmul(Lambda, p_4lab);
    
    # Center-of-Momentum Collision
    p_prime = p_4cm_photon[0];
    p_dprime = (((2*m_1*m_2)/(m_1+m_2))*(p_prime - Q +(p_prime**2/(2*M))))**.5;
    
    # inverse-Lorentz Matrix
    Lambda_prime = [[0 for x in range(dim)] for y in range(dim)];
    Lambda_prime[0][0] = np.cosh(phi);
    Lambda_prime[0][1] = np.sinh(phi);
    Lambda_prime[1][0] = np.sinh(phi);
    Lambda_prime[1][1] = np.cosh(phi);
    Lambda_prime[2][2] = 1;
    Lambda_prime[3][3] = 1;
    
    p_4cm_f = [[0 for x in range(1)] for y in range(dim)]; # length(theta)
    p_4cm_f[0][0] = (m_1**2+p_dprime**2)**.5;
    p_4cm_f[1][0] = p_dprime*np.cos(theta);
    p_4cm_f[2][0] = p_dprime*np.sin(theta);
    p_4cm_f[3][0] = 0;
    
    p_lab = np.matmul(Lambda_prime, p_4cm_f);
    spat_sum = p_lab[1]**2 + p_lab[2]**2 + p_lab[3]**2;
    
    recoil_energy = np.true_divide(spat_sum, 2*m_1);
    
    # alpha (deflection angle in lab frame) calculation
    spatial_p = np.array([[p_lab[1]],[p_lab[2]],[p_lab[3]]])
    alpha = np.arctan2(spatial_p[1],spatial_p[0]);
    
    return float(recoil_energy), float(alpha)

# End of Momentum Correction Function

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# lab_re = 0
# lab_alpha = 0

def corrections_MCNP(input_file, output_file):
    to_b_read = input_file
    to_b_written = output_file
    
    re_lab_final = []
    alpha_final = []
    lab_re = 0
    lab_alpha = 0

    with open(to_b_read, 'r') as f:
        with open(to_b_written, 'w') as g:
            index = 0
            input_energy = 0 # energy of photon; reaction-dependent
            
            for line in f:                
                if index % 3 == 0:
                    # no corrections necessary
                    g.write(line)
                
                if index % 3 == 1:
                    # no corrections necessary   
                    input_energy = 1000*float(line[2:9])
                    g.write(line)
                
                if index % 3 == 2:
                    # corrections definitely necessary
                    g.write(line[:5])

                    cos_theta = float(line[37:45])
                    CoM_theta = np.arccos(cos_theta)
                    lab_re, lab_alpha = boost(input_energy,CoM_theta)
 
                    g.write(str(round(lab_re/1000, 3)))
                
                    g.write(line[10:37])
                
                    cos_alpha = np.cos(lab_alpha)
                
                    if cos_alpha >= 0:
                        g.write(" " + str(round(cos_alpha, 5)))
                    
                    else:
                        g.write(str(round(cos_alpha, 5)))
                    
                    g.write(line[45:])
                
                
                re_lab_final.append(lab_re) 
                alpha_final.append(lab_alpha)

            
                index = index + 1