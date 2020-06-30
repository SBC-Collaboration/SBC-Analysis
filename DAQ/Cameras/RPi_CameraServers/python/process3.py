#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 12:09:32 2020

@author: pi
"""

import numpy as np

background = ls[0,:,400:799]
current = ls[1,:,400:799]
results2 = np.subtract(background,current)