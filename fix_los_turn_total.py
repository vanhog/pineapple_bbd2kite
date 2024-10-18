#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 13:21:29 2024

@author: hog
"""



import kite
import numpy as np
filen = '/media/hog/docCrucial1T/msc_grond_noemi/volume_opti/data/events/roenne/insar_float32/tl5_l2b_d_139_01_mscaoi_f32_foccov'
sc= kite.Scene.load(filen)

sc.phi = np.pi * (sc.phiDeg - 180.)/180.
sc.phiDeg = sc.phiDeg - 180.
sc.theta *= -1.
sc.thetaDeg *= -1.
sc.los.unitE *= -1.
sc.los.unitN *= -1
sc.los.unitU *= -1.

#sc.spool()
filen_out = filen + '_losturn'
sc.save(filen_out)

#save scene ... and exit