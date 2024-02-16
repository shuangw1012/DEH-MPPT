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

if __name__=='__main__':
    
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
                                         gcr=0.3)
    
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
    