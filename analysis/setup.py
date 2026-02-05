#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 08:38:35 2023

@author: louisbrezin
"""

from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize('readFile.pyx'))
