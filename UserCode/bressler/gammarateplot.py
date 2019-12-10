#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 10:01:02 2019

@author: bressler
"""

import matplotlib.pyplot as plt
#-45C data
Q_bg45 = [1.09,1.26,1.48,1.77]
r_bg45 = [3.2,5.75,5.63,3.15]
r_err_bg45 = [1.85,2.6,2.3,1.4]

Q_YBeO45 = [1.09,1.26,1.37,1.48]
r_YBeO45 = [26.5,19.0,22.1,25.0]
r_err_YBeO45 = [5.3,4.2,4.4,5.0]

Q_YAl45 = [1.09,1.26,1.37,1.48]
r_YAl45 = [8.93,9.76,10.6,9.57]
r_err_YAl45 = [2.1,2.2,2.3,2.3]
 

plt.figure()
plt.errorbar(Q_bg45,r_bg45,r_err_bg45,linestyle='none',capsize=0,elinewidth=3,
             marker='o',label="background")
plt.errorbar(Q_YBeO45,r_YBeO45,r_err_YBeO45,linestyle='none',capsize=0,elinewidth=3,
             marker='o',label="Y-BeO")
plt.errorbar(Q_YAl45,r_YAl45,r_err_YAl45,linestyle='none',capsize=0,elinewidth=3,
             marker='o',label="Y-Al")
plt.xlabel('NR threshold [keV]',fontsize=20)
plt.ylabel('Bubble Rate [mHz]',fontsize=20)
plt.legend(fontsize=18)
plt.grid()
plt.show

#-43C data:
Q_bg43 = [0.85,0.9148,0.982,1.0567]
r_bg43 = [4.24,4.66,3.61,4.40]
r_err_bg43 = [.58,1.07,.85,1.01]

r_YBeO43 = [21.1,16.5,19.1,19.4]
r_err_YBeO43 = [2.7,3.7,3.8,3.9]

r_YAl43 = [4.9,3.9,4.7,3.2]
r_err_YAl43 = [.6,.9,.9,.7]

r_Co43 = [5.38,3.20,0,2.62]
r_err_Co43 = [2.0,2.3,0,2.6]

r_BiBe43 = [5.43,2.28,5.91,5.99]
r_err_BiBe43 = [0.93,1.02,1.64,1.66]

plt.figure()
plt.errorbar(Q_bg43,r_bg43,r_err_bg43,linestyle='none',capsize=0,elinewidth=3,
             marker='o',label="background -43C")
plt.errorbar(Q_bg43,r_YBeO43,r_err_YBeO43,linestyle='none',capsize=0,elinewidth=3,
             marker='o',label="Y-BeO -43C")
plt.errorbar(Q_bg43,r_YAl43,r_err_YAl43,linestyle='none',capsize=0,elinewidth=3,
             marker='o',label="Y-Al -43C")
plt.errorbar(Q_bg43,r_Co43,r_err_Co43,linestyle='none',capsize=0,elinewidth=3,
             marker='o',label = "Cobalt -43C")
plt.errorbar(Q_bg43,r_BiBe43,r_err_BiBe43,linestyle='none',capsize=0,elinewidth=3,
             marker='o',label = "Bi-Be -43C")
plt.xlabel('NR threshold [keV]',fontsize=20)
plt.ylabel('Bubble Rate [mHz]',fontsize=20)
plt.legend(fontsize=18)
plt.grid()
plt.show

