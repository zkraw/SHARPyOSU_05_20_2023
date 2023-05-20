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
file = 'SI_2 No Aileron Def.xlsx'
global case_name
case_name = 'Span_Opt'
excel = True 

# Output File
beam_save_file = 'warren_check_new_torsion.xlsx'

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

set['alpha'] = 9.3634*np.pi/180
set['initial_thrust'] = 499.2316 
set['cs_deflection'] =  -5.3938*np.pi/180
set['g_loading'] = 1.8

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


y = [0, .507, 1.621, 2.735, 4.536, 7.239, 9.940, 12.363, 14.897, 16.219, 18.355, 20.980, 
     23.543, 24.824, 26.108, 28.675, 30.114, 31.553, 33.094, 34.634, 35.404, 36.174]

eiy = [2.81E+13, 2.72E+13, 2.68E+13, 2.63E+13, 2.36E+13, 2.23E+13, 2.21E+13, 1.86E+13, #Nmm2 
       1.66E+13, 1.26E+13, 1.00E+13, 8.53E+12, 6.34E+12, 5.95E+12, 5.14E+12, 3.34E+12, 2.54E+12, 1.75E+12, 1.54E+12, 1.27E+12, 1.03E+12, 8.19E+11] 
eiy = np.array(eiy)*(1/10**6)

gj = [2.70E+13, 1.98E+13, 1.61E+13, 1.23E+13, 1.23E+13, 1.23E+13, 1.23E+13, 1.23E+13, 1.23E+13, #Nmm2
      1.23E+13, 1.08E+13, 1.08E+13, 8.83E+12, 8.83E+12, 7.74E+12, 5.82E+12, 4.88E+12, 3.33E+12, 2.68E+12, 1.55E+12, 1.36E+12, 1.19E+12]
gj = np.array(gj)*(1/10**6)

w_a_cross = [5.4E+05, 5.4E+05, 5.4E+05, 5.4E+05, 5.4E+05, 5.4E+05, 5.4E+05, 5.4E+05, 5.4E+05, 5.4E+05, 5.4E+05, 5.4E+05, #mm2 
                5.4E+05, .4E+05, 4.9E+05, 4.0E+05, 3.6E+05, .1E+05, 2.7E+05, 2.3E+05, 2.1E+05, 1.9E+05]
w_a_cross = np.array(w_a_cross)*(1/10**6)

w_t_fwdaft = [1.59, 1.17,  0.95,  0.73,  0.73,  0.73,  0.73,  0.73,  0.73,  0.73, 0.64,  0.64,  0.52,  0.52, # mm
                0.52, 0.52 ,  0.52,  0.43,  0.43,  0.31,  0.31,  0.31]
w_t_fwdaft = np.array(w_t_fwdaft)*(1/10**3)

w_a_fwdaft = [865.3, 635.1, 515.7, 396.3, 396.3, 396.3, 396.3, 396.3, 396.3, 396.3, 346.3, 346.3, 283.4, # mm2
              283.4, 274.0, 255.1, 244.4, 192.4, 182.8, 126.5, 122.9, 119.4] 
w_a_fwdaft = np.array(w_a_fwdaft)*(1/10**6)

w_h_spar = [542.8, 542.8, 542.8, 542.8, 542.8, 542.8, 542.8, 542.8, 542.8, 542.8, 542.8, 542.8, 542.8, # mm nedd to divide by 2
            542.8, 524.9, 488.6, 468.1, 447.5, 425.2, 402.8, 391.5, 380.2] 
w_h_spar = np.array(w_h_spar)*(1/(2*10**3))

control_points = len(y)

interp = np.concatenate((eiy, gj, w_a_cross, w_t_fwdaft, w_a_fwdaft, w_h_spar), axis=0)
# -----------------------------------------------------------------------------------------------------
#----------------------END OF USER INPUT--------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

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

    logs['shear'] = []
    logs['moment'] = []
    return logs
logs = logs_def()

# Mach = set['u_inf']/a


generate_dyn_file_2([set['dt'], set['n_tstep'], route, case_name, set['amplitude'], set['period'], set['free_flight']], run.num_nodes)

from OSU_Contribution.utils.Dyn_Min_Spar_Func.init_sharpy_settings import static_trim_viscous
static_trim_viscous(route, case_name, set)

def opt(guess):
    
    y_w = [0, .507, 1.621, 2.735, 4.536, 7.239, 9.940, 12.363, 14.897, 16.219, 18.355, 20.980, 
          23.543, 24.824, 26.108, 28.675, 30.114, 31.553, 33.094, 34.634, 35.404, 36.174]

    w_shear = [6389, 6289, 5923, 5393, 4667, 3555, 2271, 5075, 3983, 7060, 6316, # N
                    5288, 4226, 3459, 2970, 2247, 1568, 1120, 692, 322, 110, 6] 

    w_moment = [129401518, 126213187, 119615051, 113607051, 105201493, 95592771, 89458690,
                    77161384, 67069693, 57735972, 44244206, 30362717, 19531953, 15100730, 11287237, 5520281, 3263296, 1651379, 585723, 89465, 4530, 0]
    w_moment = np.array(w_moment)*(1/10**3)

    w_torsion = [-9757832, -9791463, -9898694, -10044720, -10233222, -10516090, -10836618, -4863442,
                    -5117491, 818655, 665717, 475589, 301268, 200930, 145063, 67499, 8852, -21429, -38861, -40007, -27018, -12967]
    w_torsion = np.array(w_torsion)*(1/10**3)

    w_ms_tu = [4.8, 4.7, 4.9, 5.1, 4.9, 5.2, 5.6, 5.4, 5.6, 4.8, 5.0, 6.4, 7.6, 9.4, 11.5, 16.8, 
                    22.9, 33.0, 88.2, 507.9, 8352.6, 4000]

    w_ms_cu = [1.1, 1.1, 1.1, 1.2, 1.1, 1.2, 1.4, 1.3, 1.4, 1.1, 1.2, 1.7, 2.1, 2.8, 3.5, 5.4,
                    7.6, 11.2, 31.1, 181.9, 3000.5, 4000]

    w_ms_su = [23.99,12.27, 10.08, 7.88, 8.45, 9.48, 11.00, 11.35, 13.17, 11.46, 11.24, 13.73,
                    14.21, 17.71, 20.13, 25.21, 35.49, 38.84, 58.54, 82.26, 212.63, 1504.50]

    w_ms_tu_web = [3.5, 3.43, 3.60, 3.76, 3.61, 3.79, 4.08, 3.95, 4.09, 3.47, 3.64, 4.77, 5.67, 7.09, 8.68, 12.82, 17.56, 25.38, 68.14, 393.39, 6472.87, 7000]

    
    # Interpolation do stiffness and geomtries 
    f = []

    # might need to change this
    cutoff = len(guess)/control_points
    for i in range(1, int(cutoff)+1):
        if i == 1:
            start = 0
        else:
            start = control_points*(i-1)
        end = control_points*(i)
        
        f.append(interp1d(y, guess[start:end], kind= 'linear'))
    
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.h5_file_read_and_chord_vals import h5_file_read
    aero, fem = h5_file_read('OSU_Contribution/' + case_name)


    scale = 1
    # Scaling constant chord values to maintain area
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.h5_file_read_and_chord_vals import scale_chord_values
    chord_vals = scale_chord_values(aero, set['beam_opt_elems'], set['beam_opt_sym'], scale)
    

    span_init =  fem['coordinates'][set['beam_opt_nodes'][-1]][1]
    max_deflection = np.tan(25*np.pi/180)*span_init 

    y_s = []
    for i in range(len(set['beam_opt_nodes'])):
        y_s.append(fem['coordinates'][i][1])
#---------------------------------------------------------------------------------------------------------------------
#--------------UPDATING STIFFNESS & MASS PROPERTIES--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------
    
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.stiffness_and_mass import warren_interp
    ei_final, gj_final, spar_geom = warren_interp(set['beam_opt_nodes'], fem, spar_cap_dict, f, chord_vals, scale)
  
    # adding final stiffness and mass values to input dictionary 
    for i in range(len(set['beam_opt_elems'])):
        fem['stiffness_db'][i][3][3] =  gj_final[i]
        fem['stiffness_db'][i][4][4] =  ei_final[i]

    # adding symmetric values
    if set['beam_opt_sym'] == True:
        for i in range(len(set['beam_opt_elems'])):
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
    
    # Trim condition 
    # Check value of trim in utils settings 
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.sharpy_run_settings import any_g_trim
    g_load = any_g_trim(route, case_name, logs, set)
    
    print('G Loading for this case is = ' + str(g_load) + "\n\n Change g_loading in utils.user_setting.py if needed") 

    # Static Trim Information
    max_load = sharpy_main.main(['', route + '/' + case_name + '.sharpy'])
    
    if max_load.structure.div_flag.value != 1:
        logs['alpha_log'].append(max_load.alpha) # trim_calc.alpha
        logs['elevator_log'].append(max_load.elevator) # trim_calc.elevator
        logs['thrust_log'].append(max_load.thrust) # [trim_calc.thrust]    
    else:
        logs['alpha_log'].append(0) # trim_calc.alpha
        logs['elevator_log'].append(0) # trim_calc.elevator
        logs['thrust_log'].append(15) # [trim_calc.thrust]    
        
    
     #, ei_log, gj_log, mass_log])
    
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
    
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.margin_of_safety import warren_check
    m_s_fcu, m_s_ftu, m_s_su, m_s_tu_web, shear, moment, torsion, restor_moment_trans, restor_moment_copy, moment_trans, moment_copy, m_grav_trans, m_grav_copy = warren_check(max_load, spar_cap_dict, set, spar_geom, chord_vals)

#---------------------------------------------------------------------------------------------------------------------
#-------------------------FINAL LOGGING--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------     
    max_deflection_val = max_load.structure.elements[set['beam_opt_elems'][-1]].coordinates_ini[1][2] + max_deflection
    
    # max deflection across the gust encounter
    deflection = []
    for i in range(max_load.ts+1):
        deflection.append(max_load.structure.timestep_info[i].pos[set['beam_opt_nodes'][-1]][2])
    
    deflection = max(deflection)
    
    logs['mass_log'].append(masss_temp)
    
    logs['ms_tu_log'].append(m_s_ftu)
    logs['ms_su_log'].append(m_s_su) 
    logs['ms_cu_log'].append(m_s_fcu)  
    logs['ms_tu_web_log'].append(m_s_tu_web)  

    logs['shear'].append(shear) 
    logs['moment'].append(moment) 
    
    logs['wing_tip_deflection_log'].append(deflection)
    logs['max_deflection_val_log'].append(max_deflection_val)    

    df1 = pd.DataFrame([masss_temp, deflection, max_deflection_val, logs['alpha_log'], logs['elevator_log'], logs['thrust_log']])
    df2 = pd.DataFrame([y_s, shear.tolist(), moment.tolist(), torsion, m_s_ftu, m_s_fcu, m_s_su, m_s_tu_web, restor_moment_trans, restor_moment_copy, moment_trans, moment_copy, m_grav_trans, m_grav_copy])
    df3 = pd.DataFrame([y_w, w_shear, w_moment, w_torsion, w_ms_tu, w_ms_cu, w_ms_su, w_ms_tu_web])

    # writing to output file
    xlwriter = pd.ExcelWriter(beam_save_file)
    df1.to_excel(xlwriter, header= False,  sheet_name= 'Trim_Log')
    df2.to_excel(xlwriter, header= False,  sheet_name= 'Minimization_Log')
    df3.to_excel(xlwriter, header= False,  sheet_name= 'Warren Values')
    xlwriter.close()
        
    return

opt(interp)



                        
        