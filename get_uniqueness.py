#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 12:56:28 2024

@author: sudhaus
"""

import numpy as num
import sys


run_folder = '/media/hog/docCrucial1T/msc_grond_noemi/volume_opti/runs/roenne_real00044_hog.grun'
models = num.fromfile(run_folder + '/harvest/' + 'models', dtype='<f8') 


nmods = num.shape(models)[0]
uniqueness = num.shape(num.unique(models))[0]/nmods
print('Uniqueness of ', run_folder, 'is: ',uniqueness)