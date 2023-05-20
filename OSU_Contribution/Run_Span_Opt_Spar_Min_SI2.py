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
file = 'Scaled Solar Impulse Final.xlsx'
global case_name
case_name = 'Span_Opt'
excel = True 

# Output File
beam_save_file = 'Nu_Span_Opt.xlsx'

input_file = os.path.dirname(os.path.abspath(__file__)) + '/Input Files/' + file 

global spar_cap_dict, span

span = 58.5

from OSU_Contribution.utils.Dyn_Min_Spar_Func.user_setting import scaled_SI2
set, case_name, spar_cap_dict = scaled_SI2(case_name, span)

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
        'DynamicCoupled',
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

# Span Mass Spar Creation
global surface_area, mass_desired, guess_x_loc, max_deflection, mag
surface_area = 162
mass_desired = 890

max_deflection = np.tan(15*np.pi/180)*span / 2

# guess = [0.001]*6 # spar thickness
guess = [0.00208604, 0.0016987, 0.0007724, 0.00073755, 0.00050012, 0.00050005]
guess_x_loc = np.linspace(0, 1, 6)  
max_thick = .01
min_thick = .0003
mag = 1000

# bounds on minimization function
bnds = []
for i in range(len(guess_x_loc)):
    bnds.append(np.array([min_thick, max_thick])*mag)


# -----------------------------------------------------------------------------------------------------
#----------------------END OF USER INPUT--------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------
#----------------------LOG FILES--------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

logs = {}
logs['guess_log'] = []
logs['return_log'] = []
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
logs['thrust_log'] = []

# Mach = set['u_inf']/a

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

generate_dyn_file_2([set['dt'], set['n_tstep'], route, case_name, set['amplitude'], set['period'], set['free_flight']], run.num_nodes)

from OSU_Contribution.utils.Dyn_Min_Spar_Func.init_sharpy_settings import scaled_si2_dynamic
scaled_si2_dynamic(route, case_name, set)

def opt(guess):
    
    global logs
    guess = guess / mag
    guess_final = list(np.array(guess))
    
    f = interp1d(guess_x_loc, guess_final, kind= 'linear')
    
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
    XFOIL_call_constant_chord(XFOIL_PATH, set['xfoil'], set, aero)
    

#---------------------------------------------------------------------------------------------------------------------
#--------------UPDATING STIFFNESS & MASS PROPERTIES--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------
    
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.stiffness_and_mass import scaled_SI2_OG
    ei_final, gj_final, mass_final, spar_geom = scaled_SI2_OG(set['beam_opt_nodes'], fem, spar_cap_dict, f, chord_vals)
  
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
#-------------------------1G TRIM--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------
    
    # 1G trim condition 
    # needed as initial condition for dynamic coupled
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.sharpy_settings import one_g_trim
    g_load =one_g_trim(route, case_name, logs, set)
    
    # 1G Static Trim Information
    trim_calc = sharpy_main.main(['', route + '/' + case_name + '.sharpy'])

    logs['alpha_log'].append(trim_calc.alpha) # trim_calc.alpha
    logs['elevator_log'].append(trim_calc.elevator) # trim_calc.elevator
    logs['thrust_log'].append(trim_calc.thrust) # [trim_calc.thrust]  
    
    df1 = pd.DataFrame([ logs['alpha_log'], logs['elevator_log'], logs['thrust_log']]) #, ei_log, gj_log, mass_log])

#---------------------------------------------------------------------------------------------------------------------
#-------------------------GUST RESPONSE--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.sharpy_settings import turbulent_gust
    aero, fem = turbulent_gust(route, case_name, logs, set, fem)
    
    #deletes old input files
    clean_test_files(case_name)
    # creates aero input with trimmed elevator value
    write_aero(aero, case_name)
    # writing trimed condition before gust
    write_fem(fem, case_name)
    
    # GUST RESPONSE
    gust_calc = sharpy_main.main(['', route + '/' + case_name + '.sharpy'])
    
#---------------------------------------------------------------------------------------------------------------------
#-------------------------MASS CALCULATION--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------   
    
    # calculating the mass of the configuration
    grav_forces = gust_calc.structure.timestep_info[0].gravity_forces
    grav = 0
    for i in range(len(grav_forces)):
        grav += np.sqrt(grav_forces[i][0]**2 + grav_forces[i][1]**2 + grav_forces[i][2]**2)
    
    masss_temp = abs(grav/(g_load*9.81)) 
    
#---------------------------------------------------------------------------------------------------------------------
#-------------------------MS CALCULATIONS--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------
    
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.margin_of_safety import ms_OG_DYN
    m_s_fcu, m_s_ftu, m_s_su = ms_OG_DYN(gust_calc, spar_cap_dict, set, spar_geom, chord_vals)

#---------------------------------------------------------------------------------------------------------------------
#-------------------------PENALTIES--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------
        
    penalties = 0
    
    # penalties for being under a margin or safety of 1 
    for i in range(len(m_s_ftu)):
        for j in range(len(m_s_ftu[0])):
            penalties += max(0, (1 - m_s_ftu[i][j])/m_s_ftu[i][j])  
            penalties += max(0, (1 - m_s_fcu[i][j])/m_s_fcu[i][j])    
            penalties += max(0, (1 - m_s_su[i][j])/m_s_su[i][j])  
    
    
    max_deflection_val = gust_calc.structure.elements[set['beam_opt_elems'][-1]].coordinates_ini[1][2] + max_deflection
    # max deflection across the gust encounter
    deflection = []
    for i in range(gust_calc.ts):
        deflection.append(gust_calc.structure.timestep_info[i].pos[set['beam_opt_nodes'][-1]][2])
    
    deflection = max(deflection)
    
    # penalties for being over max wingtip deflection 
    penalties += max(0, 50*(deflection-max_deflection_val)/max_deflection_val)**2 
    
    # penalties if thickness is not decreasing
    # for i in range(1, len(guess)):
    #     penalties += max(0, (guess[i] - guess[i-1])/guess[i-1]) 
     
    cost = (masss_temp / mass_desired)

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
    logs['wing_tip_deflection_log'].append(deflection)
    logs['max_deflection_val_log'].append(max_deflection_val)    
    df2 = pd.DataFrame([logs['guess_log'], logs['mass_log'], logs['return_log'], logs['penalties_log'], logs['ms_tu_log'], logs['ms_su_log'], 
                        logs['ms_cu_log'], logs['wing_tip_deflection_log'], logs['max_deflection_val_log']]) #, ei_log, gj_log, mass_log])

    # writing to output file
    xlwriter = pd.ExcelWriter(beam_save_file)
    df1.to_excel( xlwriter, header= False,  sheet_name= 'Trim_Log')
    df2.to_excel( xlwriter, header= False,  sheet_name= 'Minimization_Log')
    xlwriter.close()
        
    return cost + penalties

# data = sharpy_main.main(['', route + '/' + case_name + '.sharpy'])
results = minimize(opt, guess, bounds = bnds)
print(results)

# thick = [
# [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005], # 25 m 
# [0.00071649, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005], # 30 m 
# [0.00084737, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005], # 32 m
# [0.00100547, 0.00080657, 0.00066043, 0.0005, 0.0005, 0.0005], # 34.14 m
# [.00167506181, .00122905654, .00079361084, 0.0005, 0.0005, .0005], # 36 m 
# [0.00284516, 0.00133893, 0.00123705, 0.00123378, 0.00123262, 0.00119995] # 39 m 
# ]

# spans = [34.14, 36, 39]

# for i in range(len(thick)):
#     span = spans[i]
#     span_opt(thick[i])

# opt(guess)


