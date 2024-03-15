#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 14:16:06 2023

@author: admin-shuang
"""
import matplotlib.pyplot as plt
import numpy as np
import os
plt.rcParams["font.family"] = "Times New Roman"

A = np.array([25,35,50,70,95,120,150,185,240])
C = np.array([9.875,13,14.125,19.75,28,35.25,41.25,51,67])

coefs=np.polyfit(A[:],C[:],1)

Y = coefs[1]+coefs[0]*A
print (coefs)
print (Y)

fig = plt.figure(figsize=(6,4))
plt.scatter(A,C,marker="x",label = 'Data')
x = np.linspace(20,250,100)
y = coefs[1]+coefs[0]*x
plt.plot(x,y,c = 'red', label = 'Regression')
plt.xlabel('Cross-sectional area (mm$^2$)')
plt.ylabel('Unit cost (USD/m)')
plt.savefig(os.getcwd()+'/Regression.png',dpi=100)
plt.close(fig)


C_material = np.array([6.240475732,11.56564116,16.64496313,21.46260803,26.96848792,33.85083778])
C_pipe = np.array([34.32835821,52.23880597,77.6119403,89.55223881,105.9701493,141.7910448])

coefs=np.polyfit(C_material[:],C_pipe[:],1)

Y = coefs[1]+coefs[0]*A
print (coefs)
print (Y)

fig = plt.figure(figsize=(6,4))
plt.scatter(C_material,C_pipe,marker="x",label = 'Data')
x = np.linspace(6,35,100)
y = coefs[1]+coefs[0]*x
plt.plot(x,y,c = 'red', label = 'Regression')
plt.xlabel('Material cost (USD/m)')
plt.ylabel('Pipe cost (USD/m)')
plt.savefig(os.getcwd()+'/Regression2.png',dpi=100)
plt.close(fig)