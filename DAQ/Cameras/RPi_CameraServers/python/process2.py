#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 12:07:33 2020

@author: pi
"""
import numpy as np

background = ls[0,:,0:399]
current = ls[1,:,0:399]
results1 = np.subtract(background,current)