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
    
    def calculate_intersection(voltage_values,current_values,slope):
        linear_current_values = slope * voltage_values
        
        # Find the absolute difference between the two sets of current values
        difference = np.abs(current_values - linear_current_values)
        
        # The index of the minimum difference points to the intersection
        min_index = np.argmin(difference)
        
        # Intersection point
        intersecting_voltage = voltage_values[min_index]
        intersecting_current = current_values[min_index]
        return (intersecting_voltage,intersecting_current)
    
    # technical parameters for the PV module
    parameters = {
        'Name': 'Trina Solar TSM-670DE21',
        'A_c': 3.08,
        'N_s': 44,
        'I_sc_ref': 18.62,
        'V_oc_ref': 46.1,
        'I_mp_ref': 17.55,
        'V_mp_ref': 38.2,
        'alpha_sc': 0.0059584, 
        'beta_oc': -0.098654,
        'a_ref': 1.5816,
        'I_L_ref': 18.6439,
        'I_o_ref': 4.01246E-12,
        'R_s': 0.165157,
        'R_sh_ref': 128.445,
        'Adjust': 12.2293,
        'gamma_r': -0.301,
        'Version': '2023.10.31',
    }
    
    # input the irradiance and the temperature
    start = 0
    delta_t = 24*7
    conditions = pd.read_csv(os.getcwd()+'/pv_system_output.csv').fillna(0).iloc[start:start+delta_t]
    conditions.columns = ['Time', 'Geff', 'Tcell']
    
    output_MPPT = np.array([])
    output_OP = np.array([])
    for t in range(start,start+delta_t):
        if conditions['Geff'][t] <10:
            output_MPPT = np.append(output_MPPT,0)
            output_OP = np.append(output_OP,0)
            continue
        
        # calculate the five parameters at the given input condition
        IL, I0, Rs, Rsh, nNsVth = pvsystem.calcparams_desoto(
            conditions['Geff'][t], # The irradiance (W/m2) that is converted to photocurrent.
            conditions['Tcell'][t], # The average cell temperature of cells within a module in C.
            alpha_sc=parameters['alpha_sc'], # The short-circuit current temperature coefficient of the module in units of A/C.
            a_ref=parameters['a_ref'], # The product of the usual diode ideality factor (n, unitless), number of cells in series (Ns), and cell thermal voltage at reference conditions, in units of V.
            I_L_ref=parameters['I_L_ref'], # The light-generated current (or photocurrent) at reference conditions, in amperes.
            I_o_ref=parameters['I_o_ref'], # The dark or diode reverse saturation current at reference conditions, in amperes.
            R_sh_ref=parameters['R_sh_ref'], # The shunt resistance at reference conditions, in ohms.
            R_s=parameters['R_s'], # The series resistance at reference conditions, in ohms.
            EgRef=1.121, # The energy bandgap at reference temperature in units of eV. 1.121 eV for crystalline silicon. 
            dEgdT=-0.0002677
        )
        
        # plug the parameters into the single diode equation and solve for IV curves:
        curve_info = pvsystem.singlediode(
            photocurrent=IL,
            saturation_current=I0,
            resistance_series=Rs,
            resistance_shunt=Rsh,
            nNsVth=nNsVth,
            #ivcurve_pnts=100,
            method='lambertw'
        )
        curve_v = np.linspace(0,curve_info['v_oc'],100)
        curve_i = pvsystem.i_from_v(resistance_shunt=Rsh, resistance_series=Rs, nNsVth=nNsVth, voltage=curve_v, 
                          saturation_current=I0, photocurrent=IL, method='lambertw')
        
        # maximum power point
        v_mp = curve_info['v_mp']
        i_mp = curve_info['i_mp']
        
        # calculation for the V-I limitation due to the fixed load
        R_pipe = 0.0119188 # pipe resistance
        N_panel = 10 # number of panels in a string
        slope1 = 1/(R_pipe*N_panel*1**2)
        slope2 = 1/(R_pipe*N_panel*5**2)
        
        # calculate the intersection points with V-I curve
        v2,i2 = calculate_intersection(curve_v,curve_i,slope2)
        v1,i1 = calculate_intersection(curve_v,curve_i,slope1)
        
        if v2<=v_mp:
            v_ope = v2
            i_ope = i2
        else:
            v_ope = v_mp
            i_ope = i_mp
        
        output_MPPT = np.append(output_MPPT,v_mp*i_mp)
        output_OP = np.append(output_OP,v_ope*i_ope)
        
        '''
        # plot the calculated curves:
        fig = plt.figure(figsize=(6,4))
        plt.plot(curve_v, curve_i)
        plt.plot(curve_v,curve_v*slope1,label = 'k=1')
        plt.plot(curve_v,curve_v*slope2,label = 'k=5')
        plt.plot([v_mp], [i_mp], ls='', marker='o', c='k',label = 'MPP')
        plt.plot([v_ope], [i_ope], ls='', marker='x', c='red',label = 'OP')
        plt.xlabel('Module voltage [V]')
        plt.ylabel('Module current [A]')
        plt.xlim(0,50)
        plt.ylim(0,22)
        plt.legend(loc=(0.8, 0.7))
        plt.title('G = %s W/m2, T = %s C'%(int(conditions['Geff'][t]),int(conditions['Tcell'][t])))
        plt.tight_layout()
        plt.savefig(os.getcwd()+'/Figure-VI_%s.png'%t,dpi=500)
        plt.close(fig)
        '''
        
    print (sum(output_MPPT),sum(output_OP))
    
    # plot the output
    X = np.linspace(1,len(output_MPPT),len(output_MPPT))
    fig = plt.figure(figsize=(6,4))

    plt.plot(X, output_MPPT, label = 'MPPT')
    plt.plot(X, output_OP, label = 'OP')
    
    plt.xlabel('Time (h)')
    plt.ylabel('Power output (W)')
    plt.xlim(0,len(output_MPPT))
    plt.ylim(0,700)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.getcwd()+'/output.png',dpi=500)
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
    
    