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
def boost(p_photon, cos_theta, Qval, mass_target = np.float64(9.01218307)): # reaction-dependent
    amu = np.float64(931494.0954) # keV/c^2 per amu, per physics.bu.edu
    mass_neutron = np.float64(1.0086649158) # amu
    known_Qvals = np.float64([ # keV
                              1664.5371, # 9Be
                              ])
    known_targets = np.float64([ # amu
                                9.01218307, # 9Be
                                ])
    known_gammas = np.float64([ # keV
                               1690.971, # 124Sb
                               2090.930, # 124Sb
                               1770.228, # 207Bi
                               1836.063, # 88Y
                               2734.0, # 88Y
                               1674.725, # 58Co
                               1677.85, # 148Eu
                               2173.28, # 148Eu
                               ])

    # Identify target by Q value
    Qdiff = np.abs(known_Qvals - Qval)
    target_ix = np.argmin(Qdiff)
    if Qdiff[target_ix] < 1:
        Qval = known_Qvals[target_ix]
        mass_target = known_targets[target_ix]

    # Masses
    M = mass_target*amu;
    m_1 = mass_neutron*amu;
    m_2 = M - m_1 + Qval;

    # Identify photon line
    photon_Ediff = np.abs(known_gammas - p_photon)
    photon_ix = np.argmin(photon_Ediff)
    if photon_Ediff[photon_ix] < 1:
        p_photon = known_gammas[photon_ix]

    # Rotational Lorentz "phi"
    phi = .5*np.log((2*p_photon/M) + 1);
    
    """
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
    """

    Ep_com = p_photon * np.sqrt(M/(2*p_photon+M)) # photon energy in com frame
    EM_com = np.sqrt(M*M + Ep_com*Ep_com) # target nucleus energy in com frame
    p2_com = ( # |p|^2 for either particle in com frame
              Ep_com**4 + (Ep_com**2 - 0.5*Qval**2)**2 +
              2*Qval*(M-m_1)*(.5*Qval**2 - Ep_com**2 - m_1*M) +
              2*Ep_com*EM_com*(Ep_com**2 + (M-m_1)*(m_1-Qval) - .5*Qval**2) +
              Ep_com**2 * (M**2 - 2*m_1**2 + 2*M*m_1) +
              Qval**2 * (M**2 + m_1**2 - 3*m_1*M)) / ((Ep_com + EM_com)**2)
    p_com = np.zeros(3, dtype=np.float64) # 4-momentum of neutron in com frame, 4th element=0
    p_com[0] = np.sqrt(m_1**2 + p2_com)
    p_com[1] = np.sqrt(p2_com)*cos_theta
    sin_theta = np.sqrt(1 - cos_theta**2)
    p_com[2] = np.sqrt(p2_com)*sin_theta

    p = np.zeros(3, dtype=np.float64) # 4-momentum of neutron in lab frame, 4th element=0
    p[0] = p_com[0]*np.cosh(phi) + p_com[1]*np.sinh(phi)
    p[1] = p_com[1]*np.cosh(phi) + p_com[0]*np.sinh(phi)
    p[2] = p_com[2]

    p2 = p[1]**2 + p[2]**2 # |p|^2 of neutron in lab frame

    recoil_energy = (.5 * p2 / m_1) * (1 - .25 * p2/(m_1**2)) # 1st order rel correction

    recoil_costheta = p[1] / np.sqrt(p2)
    return recoil_energy, recoil_costheta

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
                    input_energy = 1000*np.float64(line[2:9])
                    input_Qval = -1000*np.float64(line[20:29])
                    g.write(line)
                
                if index % 3 == 2:
                    # corrections definitely necessary
                    g.write(line[:5])

                    cos_theta = np.float64(line[37:45])
                    # CoM_theta = np.arccos(cos_theta)
                    lab_re, cos_alpha = boost(input_energy, cos_theta, input_Qval)
 
                    g.write(str(round(lab_re/1000, 6)))
                
                    g.write(line[10:37])
                
                    # cos_alpha = np.cos(lab_alpha)
                
                    if cos_alpha >= 0:
                        g.write(" " + str(round(cos_alpha, 5)))
                    
                    else:
                        g.write(str(round(cos_alpha, 5)))
                    
                    g.write(line[45:])
                
                
                re_lab_final.append(lab_re) 
                alpha_final.append(lab_alpha)

            
                index = index + 1
