#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 14:54:52 2023

@author: admin-shuang
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from pvlib import pvsystem
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"

if __name__=='__main__':
    
    Output_03 = pd.read_csv(os.getcwd()+'/pv_system_output_0.3.csv').fillna(0).iloc[0:24*7]
    Output_05 = pd.read_csv(os.getcwd()+'/pv_system_output_0.5.csv').fillna(0).iloc[0:24*7]
    Output_07 = pd.read_csv(os.getcwd()+'/pv_system_output_0.7.csv').fillna(0).iloc[0:24*7]
    Output_09 = pd.read_csv(os.getcwd()+'/pv_system_output_0.9.csv').fillna(0).iloc[0:24*7]
    Output = pd.read_csv(os.getcwd()+'/pv_system_output_fixed.csv').fillna(0).iloc[0:24*7]
    
    
    # plot the output
    X = np.linspace(1,len(Output_09),len(Output_09))
    fig = plt.figure(figsize=(6,4))

    plt.plot(X, Output_03['POA Global Irradiance'], label = 'GCR=0.3',linewidth=0.5)
    #plt.plot(X, Output_05['POA Global Irradiance'], label = 'GCR=0.5',linewidth=0.5)
    plt.plot(X, Output_07['POA Global Irradiance'], label = 'GCR=0.7',linewidth=0.5)
    #plt.plot(X, Output_09['POA Global Irradiance'], label = 'GCR=0.9',linewidth=0.5)
    plt.plot(X, Output['POA Global Irradiance'], label = 'Fixed tilt',linewidth=0.5,color = 'black')
    
    plt.xlabel('Time (h)')
    plt.ylabel('Total irradiance (W/m2)')
    plt.xlim(0,len(Output_03))
    plt.ylim(0,1300)
    plt.legend(ncol=4,loc='upper center')
    plt.tight_layout()
    plt.savefig(os.getcwd()+'/output_G.png',dpi=500)
    plt.close(fig)
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    '''
    Input_voltage = np.linspace(0,45,20)
    Input_current = np.zeros(10)
    
    from scipy.interpolate import interp1d
    interp_func = interp1d(curve_info['v'][0], curve_info['i'][0], kind='cubic')
        
    Input_current = interp_func(Input_voltage)
    '''
    
    '''
    fig = plt.figure(figsize=(6,4))
    for i, case in conditions.iterrows():
        #print (i,case)
        label = (
            "$G_{eff}$ " + f"{case['Geff']} $W/m^2$\n"
            "$T_{cell}$ " + f"{case['Tcell']} $C$"
        )
        plt.plot(curve_info['v'][i], curve_info['i'][i], label=label)
        #v_mp = curve_info['v_mp'][i]
        #i_mp = curve_info['i_mp'][i]
        #plt.plot([v_mp], [i_mp], ls='', marker='o', c='k')
    
    plt.legend(loc=(1.0, 0))
    plt.xlabel('Module voltage [V]')
    plt.ylabel('Module current [A]')
    plt.scatter(Input_voltage,Input_current)
    plt.title(parameters['Name'])
    #plt.show()
    #plt.gcf().set_tight_layout(True)
    plt.tight_layout()
    plt.savefig(os.getcwd()+'/Figure-VI2.png',dpi=500)
    plt.close(fig)

    print(pd.DataFrame({
        'i_sc': curve_info['i_sc'],
        'v_oc': curve_info['v_oc'],
        'i_mp': curve_info['i_mp'],
        'v_mp': curve_info['v_mp'],
        'p_mp': curve_info['p_mp'],
    }))
    '''
    
    