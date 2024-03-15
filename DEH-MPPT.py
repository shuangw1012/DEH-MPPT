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
from Material_properties.HC import Na

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
    
    def rho_Na(self,T):
        rho_1 = 7.756 + 0.02054 * T + 0.00003481 * T**2 # µΩ·cm
        rho = rho_1*10**-8 # resistivity, Ω·m
        return rho
    
    def select_pipe(self,d_min):
        D_inner = np.array([7.8, 10.418, 13.848, 17.124, 22.454, 27.862, 36.622, 42.722, 54.792, 
                            66.934, 82.804, 95.504, 108.204, 120.192, 134.492, 161.472, 
                            186.162, 211.562, 236.962, 264.668, 314.706, 342.9, 393.7, 444.5, 
                            495.3, 546.1, 596.9,644.55,695.35,746.15,796.95,847.75,898.55])/1000
        Th = np.array([1.245, 1.651, 1.651, 2.108, 2.108, 2.769, 2.769, 2.769, 2.769, 3.048, 3.048, 
                       3.048, 3.048, 3.404, 3.404, 3.404, 3.759, 3.759, 3.759, 4.191, 4.572, 6.35, 
                       6.35, 6.35, 6.35, 6.35, 6.35,7.925,7.925,7.925,7.925,7.925,7.925])/1000
        d_inner = D_inner[d_min-D_inner<0][0]
        th = Th[d_min-D_inner<0][0]
        return (d_inner,th)
    
    def DEH(self,DEH_parameters,panel_parameters):
        GCR = DEH_parameters['GCR']
        N_panel = DEH_parameters['N_panel']
        N_string = DEH_parameters['N_string']
        T_in = DEH_parameters['T_in']
        T_out = DEH_parameters['T_out']
        eta_c = DEH_parameters['eta_c']
        
        material_f = Na() 
        T_ref = (T_in+T_out)/2
        T_ref_K = (T_in+T_out)/2+273.15
        rho = self.rho_Na(T=T_ref) # resistivity
        material_f = Na() 
        
        P0 = panel_parameters['P0']
        w_panel = panel_parameters['w_panel']
        
        P_heater = N_panel*P0
        P_heater_tot = N_panel*N_string*P0*eta_c
        m_dot = P_heater_tot/material_f.Cp(T_ref_K)/(T_out-T_in)  # mass flow rate, kg/s
        V_dot = m_dot/material_f.rho(T_ref_K)
        V_limit = 3 #m/s
        d_min = np.sqrt(4*V_dot/V_limit/np.pi)
        d_inner,th = self.select_pipe(d_min)
        V = V_dot/(np.pi/4*d_inner**2) # flow velocity
        
        L_pipe = w_panel/GCR
        R_pipe = rho*L_pipe/(np.pi/4*d_inner**2)
        P_max = np.pi/4*d_inner**2*V_limit*material_f.rho(T_ref_K)*material_f.Cp(T_ref_K)*(T_out-T_in)/N_panel/N_string

        DEH_parameters.update( {'d_inner': d_inner, 'th': th, 'R_pipe': R_pipe, 
                                'P_max': P_max, 'L_pipe': L_pipe,
                                'm_dot': m_dot, 'T_ref_K': T_ref_K})
        print (12*670/material_f.Cp(T_ref_K)/(T_out-T_in))
        print ('P_heater:' + str(P_heater))
        print ('rho:' + str(rho))
        print ('L:' + str(L_pipe))
        print ('d:' + str(d_inner))
        print ('R:' + str(R_pipe))
        print (material_f.Cp(T_ref_K))
        
        return DEH_parameters
    
    def calculate_intersection(self,voltage_values,current_values,slope):
        linear_current_values = slope * voltage_values
        
        # Find the absolute difference between the two sets of current values
        difference = np.abs(current_values - linear_current_values)
        
        # The index of the minimum difference points to the intersection
        min_index = np.argmin(difference)
        
        # Intersection point
        intersecting_voltage = voltage_values[min_index]
        intersecting_current = current_values[min_index]
        return (intersecting_voltage,intersecting_current)
    
    def MPPT(self,DEH_parameters,parameters):
        
        d_inner = DEH_parameters['d_inner']
        R_pipe = DEH_parameters['R_pipe']
        P_max = DEH_parameters['P_max']
        k_max = DEH_parameters['k_max']
        
        plot_total = True
        plot_hourly = False
        
        # input the irradiance and the temperature
        start = 0#+48+8
        delta_t = 8760
        conditions = pd.read_csv(os.getcwd()+'/pv_system_output.csv').fillna(0).iloc[start:start+delta_t]
        conditions.columns = ['Time', 'Geff', 'Tcell']
        
        output_MPPT = np.array([])
        output_OP = np.array([])
        for t in range(start,start+delta_t):
            #print (t)
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
                method='lambertw'
            )
            curve_v = np.linspace(0,curve_info['v_oc_ref'],100)
            curve_i = pvsystem.i_from_v(resistance_shunt=Rsh, resistance_series=Rs, nNsVth=nNsVth, voltage=curve_v, 
                              saturation_current=I0, photocurrent=IL, method='lambertw')
            
            # maximum power point
            v_mp = curve_info['v_mp']
            i_mp = curve_info['i_mp']
            
            # calculation for the V-I limitation due to the fixed load
            #R_pipe = 0.0119188 # pipe resistance
            #N_panel = 20 # number of panels in a string
            slope1 = 1/(R_pipe*N_panel*1**2)
            slope2 = 1/(R_pipe*N_panel*k_max**2)
            
            # calculate the intersection points with V-I curve
            v2,i2 = self.calculate_intersection(curve_v,curve_i,slope2)
            v1,i1 = self.calculate_intersection(curve_v,curve_i,slope1)
            
            if v2<=v_mp:
                v_ope = v2
                i_ope = i2
            else:
                v_ope = v_mp
                i_ope = i_mp
            
            P_ope = min(v_ope*i_ope*eta_c,P_max)
            #print (P_ope,P_max)
            #print (v_ope,i_ope,P_ope)
            #print (v_mp,i_mp,v_mp*i_mp*eta_c)
            
            if plot_hourly == True:
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
                plt.title('t = %s h, G = %s W/m2, T = %s C'%(t%24,int(conditions['Geff'][t]),int(conditions['Tcell'][t])))
                plt.tight_layout()
                plt.savefig(os.getcwd()+'/Figure-VI_%s.png'%t,dpi=500)
                plt.close(fig)
            
            output_MPPT = np.append(output_MPPT,v_mp*i_mp*eta_c)
            output_OP = np.append(output_OP,P_ope)
        
        print (sum(output_MPPT),sum(output_OP))
        
        if plot_total == True:
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
        
    def DEH_cost(self,DEH_parameters,panel_parameters):
        d_inner,th,L_pipe,m_dot,T_ref_K,P_tot = DEH_parameters['d_inner'],DEH_parameters['th'],DEH_parameters['L_pipe'],DEH_parameters['m_dot'],DEH_parameters['T_ref_K'],DEH_parameters['P_tot']
        N_panel,N_string,GCR,k_max = DEH_parameters['N_panel'],DEH_parameters['N_string'],DEH_parameters['GCR'],DEH_parameters['k_max']
        I_sc,l_panel,w_panel,P0 = panel_parameters['I_sc_ref'],panel_parameters['l_panel'],panel_parameters['w_panel'],panel_parameters['P0']
        material = 'Al' # cable material
        
        # cable cost
        C_cable = 0
        for i in range(min(20,N_panel)):
            d_cable_1 = self.current_carrying_cables(I_sc*i,material=material)
            l_cable_1 = l_panel*2
            C_cable+=self.cal_cable_cost(d_cable_1,material=material)*l_cable_1
            
        for i in range(min(20,N_panel),N_panel): # chang to busbar
            A_cable,shape_cable = self.current_carrying_busbars(I_sc*i,material=material)
            l_cable_1 = l_panel*2
            C_cable += self.cal_busbar_cost(A_cable,material=material)*l_cable_1
        
        C_cable = C_cable*N_string
        
        # pipe cost
        C_pipe = self.cal_pipe_heating_cost(d_inner,th)*L_pipe*N_string
        
        # transportation pipe cost
        N_unit = np.floor(P_tot/(N_panel*N_string*P0))+1
        X_unit = l_panel*(N_panel)+1.5
        Y_unit = w_panel/GCR*(N_string)
        
        N_vert = int(np.floor(np.sqrt(X_unit/Y_unit*N_unit)))
        N_horiz = int(np.floor(N_unit/N_vert)+1)
        #print (N_vert,N_horiz,X_unit,Y_unit)
        C_pipe_trans,L_pipe_trans = self.cal_pipe_trans_cost(m_dot,N_panel,N_string,l_panel,w_panel,GCR,T_ref_K,N_vert,N_horiz,X_unit,Y_unit)
        
        C_converter = 0.035*P0*N_panel*N_string
        
        # busbar cost
        I_pipe = I_sc*k_max
        L_busbar = w_panel/GCR*N_string
        A_busbar,shape = self.current_carrying_busbars(I_pipe,material=material)
        C_busbar = self.cal_busbar_cost(A_busbar,material=material)*L_busbar
        
        DEH_parameters.update( {'N_vert': N_vert, 'N_horiz': N_horiz})
        P_heater_tot = P0*N_panel*N_string
        print ('Output power: ' + str(int(P_heater_tot)))# + ' W')
        print ('Inner diameter: ' + str(round(d_inner*1000,2)))# + ' A')
        print ('Number of panels: ' + str(N_panel))# + ' V')
        print ('Number of strings: ' + str(N_string))# + ' V')
        print ('Heater current: ' + str(int(I_pipe)))# + ' A')
        print ('Step-down ratio: ' + str(round(k_max,2)))# + ' V')
        print ('Pipe length: ' + str(round(L_pipe*N_string,2)))# + ' m')
        print ('Transporation Pipe length: ' + str(int(L_pipe_trans)))
        print ('Cable cost: ' + str(int(C_cable/P_heater_tot*1000)))# + ' ' + str(int(C_cu_thick_2)))
        print ('Pipe cost: ' + str(int(C_pipe/P_heater_tot*1000)))# + ' USD')
        print ('Transporation pipe cost: ' + str(int(C_pipe_trans/P_heater_tot*1000)))# + ' USD')
        print ('busbar cost: ' + str(int(C_busbar/P_heater_tot*1000)))# + ' ' + str(int(C_cu_thick_2)))
        print ('converter cost: ' + str(int(C_converter/P_heater_tot*1000)))# + ' USD')
        #print ('Total cost: ' + str(int(C_tot)))# + ' USD')
        #print ('Unit cost: ' + str(round(P_unit,4)))# + ' USD/W')
        print ()

        
    def current_carrying_cables(self,I,material='Cu'):
        '''
        Parameters
        ----------
        I : float
            Current, A
        material : string
            Cable material

        Returns
        -------
        d_cable : float
            cable diameter, mm
        '''
        A_cable_Cu = np.array([1.00, 1.50, 2.50, 4.00, 6.00, 10.00, 16.00, 25.00, 35.00, 50.00, 70.00, 95.00, 120.00, 150.00, 185.00, 240.00, 300.00, 400.00, 500.00, 630.00])
        A_cable_Al = np.array([1.00, 1.50, 2.50, 4.00, 6.00, 10.00, 16.00, 25.00, 35.00, 50.00, 70.00, 95.00, 120.00, 150.00, 185.00, 240.00, 300.00, 400.00, 500.00, 630.00])
        I_cable_Cu = np.array([16, 21, 30, 40, 51, 69, 92, 124, 153, 187, 238, 295, 344, 395, 459, 549, 636, 744, 867, 1014])
        I_cable_Al = np.array([10, 15, 20, 28, 35, 50, 72, 96, 119,145,184,229,267,307,357,427,495,583,685,808])
        if material=='Cu':
            I_cable = I_cable_Cu
            A_cable = A_cable_Cu
        elif material=='Al':
            I_cable = I_cable_Al
            A_cable = A_cable_Al
        try:
            A = A_cable[I-I_cable<0][0]
            d_cable = np.sqrt(A*4/np.pi)
        except:
            d_cable = np.sqrt(A_cable[-1]*4/np.pi)
            print ('Error, out of cable current range')
        return d_cable
    
    def cal_cable_cost(self,d_cable,material='Cu'):
        '''
        Parameters
        ----------
        d_cable : float
            cable diameter, mm
            
        Returns
        -------
        C_cable : float
            unit length cable cost, mm
        '''
        A_cable = np.pi/4*d_cable**2
        if material=='Cu':
            C_cable = A_cable*0.26583518712379+2.37665205443601
        elif material=='Al':
            C_material = 2700*A_cable/1000000*2.2775
            C_labour = (0.2465*A_cable+2.2039)*0.6
            C_cable = C_material+C_labour
        return C_cable
    
    def cal_busbar_cost(self,A_busbar,material='Cu'):
        '''
        Parameters
        ----------
        A_busbar : float
            busbar cross-sectional area, mm2
            
        Returns
        -------
        C_busbar : float
            unit length busbar cost, mm
        '''
        if material=='Cu':
            mass = A_busbar/1e6*8960 # kg/m
            C_material = mass*8.114 # USD/m
            C_busbar = C_material*2.2204+12.933
        elif material=='Al':
            mass = A_busbar/1e6*2700 # kg/m
            C_material = mass*2.2775 # USD/m
            C_busbar = C_material*8.535+9.1451
            
        return C_busbar
    
    def current_carrying_busbars(self,I,material='Cu'):
        '''
        Parameters
        ----------
        I : float
            Current, A
        material : string
            busbar material
        
        Returns
        -------
        A : float
            busbar cross-sectional area, mm2
        shape: string
            busbar shape, width x thickness, mm x mm
        '''
        Width = np.array([10,10,12.5,16,20, 13, 16, 20, 25, 30, 20, 25, 30, 40, 50, 
                        40, 50, 60, 80, 100, 120, 160, 200, 250, 300])
        Th = np.array([1.6,2,2,2,2, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 16, 
                        16, 16, 16, 16, 16, 16, 16, 16, 16])
        Width_str = Width.astype(str)
        Th_str = Th.astype(str)
        Shape = np.char.add(np.char.add(Width_str, " × "), Th_str)
        
        A_busbar = Width*Th
        I_busbar_Cu = np.array([98,111,134,166,201, 238.45, 281.89, 338.54, 407.83, 475.86, 516.36, 615.97, 713.46, 903.86, 1089.91, 1199.69, 1437.28, 1670.35, 2126.77, 2573.82, 3014.01, 3879.31, 4729.75, 5777.50, 6812.21])
        I_busbar_Al = np.array([80,90,109,135,164,193, 229, 275, 331, 386, 419, 500, 580, 734, 886, 975, 1168, 1358, 1729, 2092, 2450, 3154, 3845, 4697, 5538])
        if material=='Cu':
            I_busbar = I_busbar_Cu
        elif material=='Al':
            I_busbar = I_busbar_Al
        A = A_busbar[I-I_busbar<0][0]
        shape = Shape[I-I_busbar<0][0]
        
        return A,shape
    
    def cal_pipe_heating_cost(self,d_in,th,option='1'):
        '''
        Parameters
        ----------
        d_in : float
            inner pipe diameter, m
        th : float
            pipe thickness, m

        Returns
        -------
        C_pipe_u : float
            unit length pipe cost, USD/m

        '''
        d_out = d_in + 2*th
        A = np.pi/4*(d_out**2-d_in**2)
        density = 8970 # kg/m3
        mass = density * A
        unit_price = 10.77 #USD/kg
        C_material = unit_price*mass
       
        C_installation = C_material*3.77198143088368+10.1962961808475 # linear relationship
        C_insulation = 0.17*(C_material+C_installation)
        return C_material+C_installation+C_insulation
    
    def cal_pipe_trans_cost(self,m_dot,N_panel,N_string,l_panel,w_panel,GCR,T_ref_K,N_vert,N_horiz,X,Y):
        def cal_pipe(m,limit = 3):
            V_dot = m/material_f.rho(T_ref_K)
            V_limit = 3 #m/s
            d_min = np.sqrt(4*V_dot/V_limit/np.pi)
            d_pipe,th = self.select_pipe(d_min)#_mass(m)
            return d_pipe,th
        limit = 3
        ratio = 0.15 # assume cheaper cold pipe
        material_f = Na() 
        
        #print (X,Y,N_horiz,N_vert,N_horiz*X,N_vert*Y)
        C_pipe_trans = 0
        C_pipe_trans_branch = 0
        L_pipe_trans = 0
        for i in range(N_horiz-1):
            m = m_dot*(i+1)
            d_pipe,th = cal_pipe(m,limit)
            C_pipe_trans_branch += self.cal_pipe_heating_cost(d_pipe,th)*X
            L_pipe_trans += X
        C_pipe_trans_branch = C_pipe_trans_branch*(1+ratio)*N_vert/(N_vert*N_horiz)
        
        m_tot = m_dot*N_vert*N_horiz
        
        d_pipe,th = cal_pipe(m_tot,limit)
        C_pipe_trans_cold1=self.cal_pipe_heating_cost(d_pipe,th)*X*N_horiz*ratio/(N_vert*N_horiz)
        d_pipe,th = cal_pipe(m_tot/2,limit)
        C_pipe_trans_cold2=self.cal_pipe_heating_cost(d_pipe,th)*Y*N_vert*ratio/(N_vert*N_horiz)
        
        C_pipe_trans_hot = 0
        for i in range(2,int(N_vert/2)):
            m = m_dot*N_horiz*i
            d_pipe,th = cal_pipe(m,limit)
            C_pipe_trans_hot+=self.cal_pipe_heating_cost(d_pipe,th)*Y*(1+ratio)*2/(N_vert*N_horiz)
        
        C_pipe_trans = C_pipe_trans_branch + C_pipe_trans_cold1 + C_pipe_trans_cold2 + C_pipe_trans_hot
        
        #print (C_pipe_trans_branch/81480,C_pipe_trans_cold1/81480,C_pipe_trans_cold2/81480,C_pipe_trans_hot/81480)
        
        L_pipe_trans = L_pipe_trans*N_vert*2
        
        #print (C_pipe_trans/N_vert/N_horiz)
        return C_pipe_trans,L_pipe_trans
    
if __name__=='__main__':
    GCR = 0.3
    N_panel = 10 # number of panels per string
    N_panel_tot = 120
    N_string = int(N_panel_tot/N_panel)
    T_in = 200
    T_out = 600
    eta_c = 0.97
    
    # technical parameters for the PV module
    panel_parameters = {
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
        'P0': 670,
        'w_panel': 2.278,
        'l_panel': 1.134,
    }
    
    DEH_parameters = {
        'N_panel': N_panel,
        'N_string': N_string,
        'N_panel_tot':N_panel_tot,
        'T_in':T_in,
        'T_out':T_out,
        'eta_c':eta_c,
        'GCR':GCR,
        'P_tot': 100e6,
        'k_max': 5
    }
    
    DEH_system = DEH_MPPT(folder = os.getcwd())
    #DEH_system.annual_simulation(GCR)
    DEH_parameters = DEH_system.DEH(DEH_parameters,panel_parameters)
    #DEH_system.DEH_cost(DEH_parameters,panel_parameters)
    #DEH_system.MPPT(DEH_parameters,panel_parameters)
    
    