#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 15:33:04 2020

@author: pi
"""

from setuptools import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize("count.pyx"))