#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 14:54:52 2023

@author: admin-shuang
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import pvlib
from pvlib import pvsystem
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"

class DEH_MPPT():
    def __init__(self,folder):
        self.folder = folder
    
    def annual_simulation(self,GCR):
        # load weather data
        filename = os.getcwd()+'/AUS_ACT.Canberra.949260_IWEC.epw'
        weather, metadata = pvlib.iotools.read_epw(filename)
        
        # input location info
        coordinates = (-35.2809, 149.1300, 'Canberra', 577, 'Australia/Sydney')
        times  = times = weather.index
        latitude, longitude, name, altitude, timezone = coordinates
        
        # calculate solar positions
        solpos = pvlib.solarposition.get_solarposition(
            time=times,
            latitude=latitude,
            longitude=longitude,
            altitude=altitude,
            temperature=weather['temp_air'],
            pressure=pvlib.atmosphere.alt2pres(altitude),
        )
        
        # get tracking positions for a single axis tracking system
        tracker_data = pvlib.tracking.singleaxis(solpos['apparent_zenith'], solpos['azimuth'],
                                             axis_tilt=0, axis_azimuth=180, 
                                             max_angle=45, backtrack=True, 
                                             gcr=GCR)
        
        dni_extra = pvlib.irradiance.get_extra_radiation(weather.index)
        
        # calculate total irradiance
        total_irradiance = pvlib.irradiance.get_total_irradiance(
            tracker_data['surface_tilt'],
            tracker_data['surface_azimuth'],
            solpos['apparent_zenith'],
            solpos['azimuth'],
            weather['dni'],
            weather['ghi'],
            weather['dhi'],
            dni_extra=dni_extra,
            model='haydavies',
        )
        
        # calculate cell temperatures
        temperature_model_parameters = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS['sapm']['open_rack_glass_glass']
        cell_temperature = pvlib.temperature.sapm_cell(
            total_irradiance['poa_global'],
            weather["temp_air"],
            weather["wind_speed"],
            **temperature_model_parameters,
        )
        
        # ASSUME that the effective irradiance is equal to the total irradiance
        effective_irradiance_proxy = total_irradiance['poa_global']
        
        # output
        data_to_export = pd.DataFrame({
            'POA Global Irradiance': total_irradiance['poa_global'],
            'Cell Temperature': cell_temperature
        }, index=times)
        csv_filename = 'pv_system_output.csv'
        data_to_export.to_csv(csv_filename)
    
    def DEH(self,N_panel,N_panel_unit,GCR,T_in,T_out):
        
        
    
    def MPPT(self):
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
        delta_t = 8760
        conditions = pd.read_csv(os.getcwd()+'/pv_system_output.csv').fillna(0).iloc[start:start+delta_t]
        conditions.columns = ['Time', 'Geff', 'Tcell']
        
        output_MPPT = np.array([])
        output_OP = np.array([])
        for t in range(start,start+delta_t):
            print (t)
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
            N_panel = 20 # number of panels in a string
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
        


if __name__=='__main__':
    GCR = 0.3
    N_panel = 10 # number of panels per string
    N_panel_unit = 120
    T_in = 200
    T_out = 600
    
    DEH_system = DEH_MPPT(folder = os.getcwd())
    DEH_system.DEH(N_panel,N_panel_unit,GCR,T_in,T_out)
    #DEH_system.annual_simulation(GCR)
    
        