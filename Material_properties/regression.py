#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 14:16:06 2023

@author: admin-shuang
"""
import numpy as np

T = np.array([90,100,150,200,250,
              273,293,300,350,400,
              500,600,700,800,900,
              1000,1100,1200,1300,1400,
              1500,1600])
rho = np.array([58.5,59.5,64.5,69.1,73.6,
                75.4,77.7,77.7,81.5,85.2,
                91.7,97.7,103.1,108.0,112.1,
                115.7,118.9,121.7,124.3,126.8,
                129.2,131.3])

coefs=np.polyfit(T[:],rho[:],2)

Y = coefs[2]+coefs[1]*T+coefs[0]*T**2
print (coefs)
print (Y)

