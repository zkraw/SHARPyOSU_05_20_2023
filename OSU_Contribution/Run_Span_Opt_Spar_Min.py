from audioop import reverse
from cProfile import label
from tkinter import font
import h5py as h5
import numpy as np
import configparser
import os
import configobj
import json
import sys
from pathlib import Path
import pandas as pd
from pyrsistent import b
from scipy.interpolate import interp1d
import sympy as sp
from scipy.optimize import minimize
import subprocess

 
# sys.path is a list of absolute path strings
direct = __file__.split('/')
del direct[-1]; del direct[-1]
direct =  '/'.join(direct)
sys.path.append(direct)

import sharpy.utils.algebra as algebra
from OSU_Contribution.Main_Paramaterized import Paramaterized
import matplotlib.pyplot as plt
import sharpy.sharpy_main as sharpy_main
from OSU_Contribution.utils.SHARPy_Input_Writes import generate_dyn_file, generate_dyn_file_2
from OSU_Contribution.utils.SHARPy_Input_Writes import clean_test_files, write_aero, write_fem
from OSU_Contribution.utils.Rotation_Matrices import local_for_forces, nodal_a_for_2_b_for 


home = str(Path.home())
route = os.path.dirname(os.path.realpath(__file__)) 

# -----------------------------------------------------------------------------------------------------
#----------------------INPUTS--------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
file = 'T_Tail_Flat_High_Payload.xlsx'
global case_name
case_name = '46_2G'
excel = True 

input_file = os.path.dirname(os.path.abspath(__file__)) + '/Input Files/' + file 

global spar_cap_dict

from OSU_Contribution.utils.Dyn_Min_Spar_Func.user_setting import t_tail
set, case_name, spar_cap_dict = t_tail(case_name)

set['flow'] = ['BeamLoader',
        'AerogridLoader',
        #'NonLinearStatic',
        #'StaticUvlm',
        #'BeamOpt',
        'StaticTrim',
        #'StaticCoupled',
        'BeamLoads',
        'AerogridPlot',
        'BeamPlot',
        #'DynamicCoupled',
        #'Modal',
        #'LinearAssember',
        # 'AsymptoticStability',
         ]

# saving values to new file
#Paramaterized.excel_write(file, 'new_file.xlsx', ['start_Y'], [[3, 4, 5, 6, 7, 8]], [[123456789, 100, 100, 100, 100, 100]])

# generate FEM and Aero from input file
run = Paramaterized(case_name, input_file, set['xfoil'], excel)

set['nodes_app'] = run.app_nodes # where the thrust nodes are 
set['tail_cs_index'] = run.tail_index # the control surface index that will be used 
  
# OPTIMIZATION SETTINGS
save_data = True
set['beam_opt_nodes'] = run.beam_opt_nodes
set['beam_opt_elems'] = run.beam_opt_elems
set['beam_opt_sym'] = run.beam_opt_sys

set['alpha'] = 10.3161907367023*np.pi/180
set['cs_deflection'] = -9.14764443148759*np.pi/180
set['initial_thrust'] = 88.6765895305575

spar_cap_dict['core'] = .005

# Output File
beam_save_file = '49_32m2_higher_load_MS_1.xlsx'

# Span Mass Spar Creation
global span, surface_area, guess_x_loc, max_deflection, mag, control_points, min_ms
surface_area = 32
span = 46
min_ms = .8 

# used to scale max height
chord = surface_area/span

max_deflection = np.tan(20*np.pi/180)*span / 2

# cap_width = [.05]*6
# spar_width = [.1]*6 
# spar_height =[.1]*6 

# cap_width = [0.05]
# spar_width = [0.06705243, 0.06568396, 0.06673288, 0.06676578, 0.06688279, 0.06848254]
# spar_height = [0.04486957, 0.04486957, 0.03610573, 0.04116199, 0.04197148, 0.04382518]

# # 49 MS 1
# cap_width = [0.04030612, 0.03356816, 0.02821992, 0.02757655, 0.02757655, 0.02757655]
# spar_width = [0.06213276, 0.06211253, 0.06211253, 0.06211253, 0.062104, 0.062104]
# spar_height = [0.04212245, 0.04212245, 0.03983358, 0.03340894, 0.03145168, 0.03145168]

# # 46 MS 1
cap_width = [0.04456522, 0.02193532, 0.02152951, 0.02130852, 0.02130852, 0.02130852]
spar_width = [0.0672121, 0.06614813, 0.06614813, 0.06614813, 0.06614813, 0.06614813]
spar_height = [0.04486957, 0.04151363, 0.04083853, 0.04042375, 0.04030649, 0.04030649]

# 42 MS 1
# cap_width = [0.03479575, 0.0226716, 0.0226716, 0.0226716, 0.0226716, 0.0226716]
# spar_width =  [0.07220695, 0.07013747, 0.07013747, 0.07013747, 0.07013747, 0.07013747]
# spar_height = [0.04750399, 0.04281908, 0.04281908, 0.04281908, 0.04281908, 0.04281908]

# # 39 MS 1
# cap_width = [0.03171846, 0.01881422, 0.00923001, 0.00923001, 0.00923001, 0.00923001]
# spar_width = [0.07695528, 0.07444593, 0.07444593, 0.07444593, 0.07444593, 0.07444593]
# spar_height = [0.05051348, 0.04898298, 0.04459268, 0.04459268, 0.04459268, 0.04459268]

control_points = len(cap_width)

guess = cap_width + spar_width + spar_height 
guess_x_loc = np.linspace(0, 1, control_points)  

buffer = .025

# WEB IS ASSUMED TO BE 6 mm
min_cap_width = spar_cap_dict['core'] + spar_cap_dict['thickness']

# max spar height and width found from area optimization script
max_spar_width = (.2*chord)/2
min_spar_width = min_cap_width + buffer

max_cap_width = max_spar_width - buffer 

max_spar_height = (.129*chord)/(2)
min_spar_height = min_cap_width + buffer

min_bnd = [min_cap_width]*control_points + [min_spar_width]*control_points + [min_spar_height]*control_points
max_bnd = [max_cap_width]*control_points + [max_spar_width]*control_points + [max_spar_height]*control_points

# -----------------------------------------------------------------------------------------------------
#----------------------END OF USER INPUT--------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

#getting guess into the correct magnitudes
mag = 1000
guess = np.array(guess)*mag

# checking if output file write name already exists
count = 1
check = True
while check == True:
    if Path(beam_save_file).exists():
        save_elem = beam_save_file.split('.')
        if count == 1:
            save_elem[0] += '(' + str(count) + ')'
            beam_save_file = save_elem[0] + ('.') + save_elem[1]
            count += 1
        else:
            save_elem[0] += '(' + str(count) + ')'     
            beam_save_file = save_elem[0] + ('.') + save_elem[1]
            count += 1
    else:
        check = False


# -----------------------------------------------------------------------------------------------------
#----------------------LOG FILES--------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

def logs_def():
    logs = {}
    logs['guess_log'] = []
    logs['return_log'] = []
    logs['spar_mass'] = []
    logs['mass_log'] = []
    logs['max_strain_log'] = []
    logs['wing_tip_deflection_log'] = []
    logs['max_deflection_val_log'] = []
    logs['wing_deflection_diff_log'] = []
    logs['power_log'] = []
    logs['chord_log'] = []
    logs['penalties_log'] = []

    logs['ei_log'] = []
    logs['gj_log'] = []

    logs['alpha_log'] = []
    logs['elevator_log'] = []

    logs['ms_tu_log'] = []
    logs['ms_su_log'] = []
    logs['ms_cu_log'] = [] 
    logs['ms_tu_web_log'] = []
    logs['thrust_log'] = []
    return logs
logs = logs_def()

generate_dyn_file_2([set['dt'], set['n_tstep'], route, case_name, set['amplitude'], set['period'], set['free_flight']], run.num_nodes)

from OSU_Contribution.utils.Dyn_Min_Spar_Func.init_sharpy_settings import static_trim_viscous
static_trim_viscous(route, case_name, set)

def opt(guess):
    
    global logs
    guess = guess / mag
    guess_final = list(np.array(guess))
    
    
    f = []

    cutoff = len(guess)/control_points
    for i in range(1, int(cutoff)+1):
        if i == 1:
            start = 0
        else:
            start = control_points*(i-1)
        end = control_points*(i)
        
        f.append(interp1d(guess_x_loc, guess_final[start:end], kind= 'linear'))
    
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.h5_file_read_and_chord_vals import h5_file_read
    aero, fem = h5_file_read('OSU_Contribution/' + case_name)

#---------------------------------------------------------------------------------------------------------------------
#--------------SPAN_OPT SCALING AND AERO CONSIDERATIONS--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------
    
    # scaling wingspan 
    span_init =  fem['coordinates'][set['beam_opt_nodes'][-1]][1]
    scale = span / (2*span_init)
    
    # COORDINATES scaling update
    for i in range(len(set['beam_opt_nodes'])):
        fem['coordinates'][i][1] = scale * fem['coordinates'][i][1]
    
    if set['beam_opt_sym'] == True:
       for i in range(len(set['beam_opt_nodes'])):
           fem['coordinates'][set['beam_opt_nodes'][-1] + i + 1][1] =  scale * fem['coordinates'][set['beam_opt_nodes'][-1] + i + 1][1]

    # Scaling constant chord values to maintain area
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.h5_file_read_and_chord_vals import scale_chord_values
    chord_vals = scale_chord_values(aero, set['beam_opt_elems'], set['beam_opt_sym'], scale)
    
    if logs['alpha_log'] == []: 
    #UPDATING REYNOLDS NUMBER DUE TO CHANGE IN CHORD
        chord_ave = sum(chord_vals) / len(chord_vals)
        # reynolds number calculation
        set['xfoil']['Re'] = (set['xfoil']['rho'] * set['xfoil']['u_inf'] * chord_ave) / set['xfoil']['u']
        
        XFOIL_PATH = route + "/XFOIL6.99/xfoil.exe"
        
        # Removing previous polar_file.txt file so new one can be calculated
        if os.path.exists("polar_file.txt"):
            os.remove("polar_file.txt")
    
        # Calling Xfoil and updating drag polar 
        from OSU_Contribution.utils.Dyn_Min_Spar_Func.XFOIL_call import XFOIL_call_constant_chord
        XFOIL_call_constant_chord(XFOIL_PATH, set['xfoil'])
    
    # Updating the drag polar for the change in chord
    nu_pol = np.loadtxt("polar_file.txt", skiprows=12)
    new_polar = []
    for i in range(len(nu_pol)):
        new_polar.append([nu_pol[i][0]*np.pi / 180, nu_pol[i][1], nu_pol[i][2], nu_pol[i][4]])
    
    pol_ref = aero['airfoil_distribution'][set['beam_opt_elems'][0]][0]
    aero['polars'][pol_ref] = new_polar

#---------------------------------------------------------------------------------------------------------------------
#--------------UPDATING STIFFNESS & MASS PROPERTIES--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------
    
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.stiffness_and_mass import t_tail_three_variable
    ei_final, gj_final, mass_final, spar_geom, m_span = t_tail_three_variable(set['beam_opt_nodes'], fem, spar_cap_dict, f, chord_vals, scale, max_spar_height)
    
    # calculating mass of spar for minimzation
    m_spar = 0
    for i in range(len(m_span)-1):
        m_spar += m_span[i]*(fem['coordinates'][i+1][1] - fem['coordinates'][i][1]) 
    
    logs['spar_mass'].append(m_spar) 
    
    # adding final stiffness and mass values to input dictionary 
    for i in range(len(set['beam_opt_elems'])):
        fem['mass_db'][i][0][0] =  mass_final[i]
        fem['mass_db'][i][1][1] =  mass_final[i]
        fem['mass_db'][i][2][2] =  mass_final[i]
        
        fem['stiffness_db'][i][3][3] =  gj_final[i]
        fem['stiffness_db'][i][4][4] =  ei_final[i]

    # adding symmetric values
    if set['beam_opt_sym'] == True:
        for i in range(len(set['beam_opt_elems'])):
            fem['mass_db'][set['beam_opt_elems'][-1] + i + 1][0][0] =  mass_final[i]
            fem['mass_db'][set['beam_opt_elems'][-1] + i + 1][1][1] =  mass_final[i]
            fem['mass_db'][set['beam_opt_elems'][-1] + i + 1][2][2] =  mass_final[i] 
            
            fem['stiffness_db'][set['beam_opt_elems'][-1] + i + 1][3][3] =  gj_final[i]
            fem['stiffness_db'][set['beam_opt_elems'][-1] + i + 1][4][4] =  ei_final[i]
                
    #deletes old input files
    clean_test_files(case_name)
    # creates new aero file with updated twist
    write_fem(fem, case_name)
    # creates new aero file with updated twist
    write_aero(aero, case_name)
    
#---------------------------------------------------------------------------------------------------------------------
#-------------------------2G MAX LOADING--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------
    
    # 2G trim condition 
    # changes the settings in the SHARPy file to reflect 2G static trim run case 
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.sharpy_run_settings import two_g_trim
    g_load = two_g_trim(route, case_name, logs, set)
    
    # 2G Static Trim Information
    max_load = sharpy_main.main(['', route + '/' + case_name + '.sharpy'])
    
    if max_load.structure.div_flag.value != 1:
        logs['alpha_log'].append(max_load.alpha) # trim_calc.alpha
        logs['elevator_log'].append(max_load.elevator) # trim_calc.elevator
        logs['thrust_log'].append(max_load.thrust) # [trim_calc.thrust]    
    else:
        logs['alpha_log'].append(0) # trim_calc.alpha
        logs['elevator_log'].append(0) # trim_calc.elevator
        logs['thrust_log'].append(15) # [trim_calc.thrust]    
        
    
    df1 = pd.DataFrame([ logs['alpha_log'], logs['elevator_log'], logs['thrust_log']]) #, ei_log, gj_log, mass_log])
    
#---------------------------------------------------------------------------------------------------------------------
#-------------------------MASS CALCULATION--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------   
    
    # calculating the mass of the configuration
    grav_forces = max_load.structure.timestep_info[0].gravity_forces
    grav = 0
    for i in range(len(grav_forces)):
        grav += np.sqrt(grav_forces[i][0]**2 + grav_forces[i][1]**2 + grav_forces[i][2]**2)
    
    masss_temp = abs(grav/(g_load*9.81)) 
    
#---------------------------------------------------------------------------------------------------------------------
#-------------------------MS CALCULATIONS--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------
    
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.margin_of_safety import ms_OG_STATIC
    m_s_fcu, m_s_ftu, m_s_su, m_s_tu_web = ms_OG_STATIC(max_load, spar_cap_dict, set, spar_geom, chord_vals)

#---------------------------------------------------------------------------------------------------------------------
#-------------------------PENALTIES--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------
        
    penalties = 0
    if max_load.structure.div_flag.value != 1:
        # penalties for being under a margin or safety of .2 
        for i in range(len(m_s_ftu)): 
            penalties += max(0, (min_ms - m_s_ftu[i])/m_s_ftu[i])
            penalties += max(0, (min_ms - m_s_fcu[i])/m_s_fcu[i])
            penalties += max(0, (min_ms - m_s_su[i])/m_s_su[i])
            penalties += max(0, (min_ms - m_s_tu_web[i])/m_s_tu_web[i])
        
        
        max_deflection_val = max_load.structure.elements[set['beam_opt_elems'][-1]].coordinates_ini[1][2] + max_deflection
        
        # max deflection across the gust encounter
        deflection = []
        for i in range(max_load.ts+1):
            deflection.append(max_load.structure.timestep_info[i].pos[set['beam_opt_nodes'][-1]][2])
        
        deflection = max(deflection)
        
        # penalties for being over max wingtip deflection 
        penalties += (max(0, 50*(deflection-max_deflection_val)/max_deflection_val)**2)
        
        cost = m_spar/(logs['spar_mass'][0])

    else:
        # for i in range(len(m_s_ftu)): 
        #     penalties += max(0, (min_ms - m_s_ftu[i])/m_s_ftu[i])
        #     penalties += max(0, (min_ms - m_s_fcu[i])/m_s_fcu[i])
        #     penalties += max(0, (min_ms - m_s_su[i])/m_s_su[i])
        deflection = 0
        max_deflection_val = max_load.structure.elements[set['beam_opt_elems'][-1]].coordinates_ini[1][2] + max_deflection 
        penalties = 0
        if logs['return_log'] != []:
            cost = max(logs['return_log']) + .1
        else:
            cost = 1.5

#---------------------------------------------------------------------------------------------------------------------
#-------------------------FINAL LOGGING--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------     
    
    logs['guess_log'].append(guess)
    logs['mass_log'].append(masss_temp)
    logs['return_log'].append(cost + penalties)
    logs['penalties_log'].append(penalties)
    
    logs['ms_tu_log'].append(m_s_ftu)
    logs['ms_su_log'].append(m_s_su) 
    logs['ms_cu_log'].append(m_s_fcu)  
    logs['ms_tu_web_log'].append(m_s_tu_web)  
    logs['wing_tip_deflection_log'].append(deflection)
    logs['max_deflection_val_log'].append(max_deflection_val)    
    
    logs['power_log'].append(logs['thrust_log'][-1]*set['u_inf'])
    df2 = pd.DataFrame([logs['guess_log'], logs['spar_mass'], logs['mass_log'], logs['return_log'], logs['penalties_log'], logs['ms_tu_log'], logs['ms_su_log'], 
                        logs['ms_cu_log'], logs['ms_tu_web_log'], logs['wing_tip_deflection_log'], logs['max_deflection_val_log'], logs['power_log']]) #, ei_log, gj_log, mass_log])

    # writing to output file
    xlwriter = pd.ExcelWriter(beam_save_file)
    df1.to_excel( xlwriter, header= True,  sheet_name= 'Trim_Log')
    df2.to_excel( xlwriter, header= True,  sheet_name= 'Minimization_Log')
    xlwriter.close()
        
    return cost + penalties

# -----------------------------------------------------------------------------------------------------
#----------------------BONDS & CONS FOR MINIMZATION----------------------------------------------------
# -----------------------------------------------------------------------------------------------------

cutoff = len(guess)/control_points
# bounds on minimization function
bnds = []

for i in range(len(guess)):
    bnds.append(np.array([min_bnd[i], max_bnd[i]])*mag)

# constraints 
class non_increase():
    def __init__(self, id_num):
        self.id_num = id_num

    # preceding variable value must be greater than or equal to the next one
    def decrease_along_span(self, guess):
        first = guess[self.id_num]
        second = guess[self.id_num+1]
        return first - second

class width_or_height():
    def __init__(self, interest):
        self.first = interest[0]
        self.second = interest[1] 

    def min_buffer(self, guess):
        cap = guess[self.first]
        interest = guess[self.second]
        return interest - (cap + buffer) 

cons = []
count = 1
# non-increasing constraint
for i in range(len(guess)):
        if i != (control_points*count) -1:
            cons.append({'type': 'ineq', 'fun': non_increase(i).decrease_along_span})
        else:
            count += 1

# minimum web width
for i in range(control_points):
    cap = i
    width = i + control_points 
    cons.append({'type': 'ineq', 'fun': width_or_height([cap, width]).min_buffer})

opt(guess)
# results = minimize(opt, guess, method = 'SLSQP', bounds = bnds, constraints = cons)
# print(results)






                        
        