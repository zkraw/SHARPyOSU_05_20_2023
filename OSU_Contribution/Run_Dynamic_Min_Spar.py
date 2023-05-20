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
file = 'Scaled_Solar_Impulse_Updated_Debug.xlsx'
global case_name
case_name = 'Solar_Impulse'
excel = True 

input_file = os.path.dirname(os.path.abspath(__file__)) + '/Input Files/' + file 

flow = ['BeamLoader',
        'AerogridLoader',
        # 'NonLinearStatic',
        # 'StaticUvlm',
        #'BeamOpt',
        'StaticTrim',
        'StaticCoupled',
        'BeamLoads',
        'AerogridPlot',
        'BeamPlot',
        'DynamicCoupled',
        # 'Modal',
        # 'LinearAssembler',
        # 'AsymptoticStability',
         ]


# FLIGHT CONDITIONS for AERO SOLVER
u_inf = 23.7 #14 for XHALE, 10 for T_Tail
rho = 0.1736 # 1.225 

# default elevator is -0.020708598 radians 3
trim_set = {}
trim_set['alpha'] = 4.7319*np.pi/180 # 2.2953
trim_set['roll'] = 0
trim_set['beta'] = 0
trim_set['cs_deflection'] = -0.6246*np.pi/180 # 2.1124
trim_set['initial_thrust'] = 67.3539 # 37.07

gravity = 'on'
g_loading = 1 # acceleration of gravity (9.81) * factor 

# gust settings
gust_intensity =  0.18*u_inf
gust_length = u_inf / .307
gust_offset = 0.1*u_inf

#DYNAMIC FLIGHT SETTINGS
free_flight = True
amplitude = 0
period = 0
if not free_flight:
    case_name += '_prescribed'
    amplitude = 0*np.pi/180
    period = 3
    case_name += '_amp_' + str(amplitude).replace('.', '') + '_period_' + str(period)

# Time for dynamic solver  
physical_time = .1
dt = .05
# dt = 1/m/u_inf*tstep_factor
n_tstep = round(physical_time/dt)

m = 8
tstep_factor = 1

num_cores= 4
# n_step = 100

#numerics
n_step = 5
fsi_tolerance = 1e-3
structural_relaxation_factor = 0.3
relaxation_factor = 0.5
tolerance = 1e-2

num_cores = 4

# XFIOL INPUTS

# speed of sound 
a = 343         #speed of sound (m/s)
u = 1.802e-05   #dynamic viscosity (kg/(m*s))

xfoil = {
    'Mach': u_inf/a,
    'n_iter': 100,
    'alfa_1': -16,
    'alfa_2': 16,
    'alfa_step': .5, 
    'u_inf': u_inf,
    'rho': rho,
    'u': u,
}

# saving values to new file
#Paramaterized.excel_write(file, 'new_file.xlsx', ['start_Y'], [[3, 4, 5, 6, 7, 8]], [[123456789, 100, 100, 100, 100, 100]])

# generate FEM and Aero from input file
run = Paramaterized(case_name, input_file, xfoil, excel)

trim_set['nodes_app'] = run.app_nodes # where the thrust nodes are 
trim_set['tail_cs_index'] = run.tail_index # the control surface index that will be used 
  

# OPTIMIZATION SETTINGS

save_data = True
trim_set['beam_opt_nodes'] = run.beam_opt_nodes
trim_set['beam_opt_elems'] = run.beam_opt_elems
trim_set['beam_opt_sym'] = run.beam_opt_sys
beam_save_file = 'Solar_Impulse_Min.xlsx'

# Span Mass Spar Creation
global span, surface_area, mass_desired, guess_x_loc, max_deflection, mag
surface_area = 162
span = 58.5
mass_desired = 890

max_deflection = np.tan(15*np.pi/180)*span / 2

# guess = [0.001]*6 # spar thickness
guess = [0.00175577, 0.00132001, 0.00080799, 0.00074273, 0.00044357, 0.0002]
guess_x_loc = np.linspace(0, 1, 6)  
max_thick = [.002, .002, .002, .0015, .00075, .00075]
min_thick = [.0002]*6
mag = 1000

# bounds on minimization function
bnds = []
for i in range(len(guess_x_loc)):
    bnds.append(np.array([min_thick[i], max_thick[i]])*mag)


# -----------------------------------------------------------------------------------------------------
#----------------------SPAR LOOKUP VALUES--------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

global spar_cap_dict

spar_cap_dict = {}
spar_cap_dict['E_cap'] = 1.783e11
spar_cap_dict['E_lam'] = 1.212e10
spar_cap_dict['G'] = 4.51e10
spar_cap_dict['Misc_Factor'] = .59

spar_cap_dict['x_taper_start'] = span * 0.686 / 2

# height and width of spar before taper
spar_cap_dict['h_spar_1'] = 0.388
spar_cap_dict['b_spar_1'] = 0.703
#height and widrh of spar after taper
spar_cap_dict['h_spar_2'] = 0.357
spar_cap_dict['b_spar_2'] = 0.477

spar_cap_dict['density'] = 1540 

# margin of safety calculations
spar_cap_dict['F_tu'] = 1.283e9
spar_cap_dict['F_cu'] = 4.61e8
spar_cap_dict['F_su'] = 2.35e8



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

Mach = u_inf/a

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

generate_dyn_file_2([dt, n_tstep, route, case_name, amplitude, period, free_flight], run.num_nodes)

def generate_solver_file():
    file_name = route + '/' + case_name + '.sharpy'
    settings = dict()
    settings['SHARPy'] = {'case': case_name,
                        'route': route,
                        'flow': flow,
                        'write_screen': 'on',
                        'write_log': 'on',
                        'log_folder': route + '/output',
                        'log_file': case_name + '.log'}

    settings['BeamLoader'] = {'unsteady': 'on',
                            'orientation': algebra.euler2quat(np.array([trim_set['roll'],
                                                                        trim_set['alpha'],
                                                                        trim_set['beta']]))}
    settings['AerogridLoader'] = {'unsteady': 'on',
                                'aligned_grid': 'on',
                                'mstar': int(20/tstep_factor),
                                'freestream_dir': ['1', '0', '0'],
                                'wake_shape_generator': 'StraightWake',
                                'wake_shape_generator_input': {'u_inf': u_inf,
                                                               'u_inf_direction': ['1', '0', '0'],
                                                               'dt': dt}}

    settings['NonLinearStatic'] = {'print_info': 'off',
                                'max_iterations': 150,
                                'num_load_steps': 1,
                                'delta_curved': 1e-1,
                                'min_delta': tolerance,
                                'gravity_on': gravity,
                                'gravity': 9.81*g_loading}

    settings['StaticUvlm'] = {'print_info': 'on', #previously off
                            'horseshoe': 'off', # previously on
                            'num_cores': num_cores,
                            'n_rollup': 0,
                            'rollup_dt': dt,
                            'rollup_aic_refresh': 1,
                            'rollup_tolerance': 1e-4,
                            'velocity_field_generator': 'SteadyVelocityField',
                            'velocity_field_input': {'u_inf': u_inf,
                                                    'u_inf_direction': [1., 0, 0]},
                            'rho': rho}


    settings['StaticCoupled'] = {'print_info': 'off',
                                'structural_solver': 'NonLinearStatic',
                                'structural_solver_settings': settings['NonLinearStatic'],
                                'aero_solver': 'StaticUvlm',
                                'aero_solver_settings': settings['StaticUvlm'],
                                'max_iter': 100,
                                'n_load_steps': n_step,
                                'tolerance': fsi_tolerance,
                                'relaxation_factor': structural_relaxation_factor}
                                # 'correct_forces_method': 'PolarCorrection',
                                # 'correct_forces_settings': {'cd_from_cl': 'off',  # recommended settings (default)
                                #                                 'correct_lift': 'off',
                                #                                 'moment_from_polar': 'off'}}




    settings['StaticTrim'] = {'solver': 'StaticCoupled',
                            'solver_settings': settings['StaticCoupled'],
                            'initial_alpha': trim_set['alpha'],
                            'initial_deflection': trim_set['cs_deflection'],
                            'tail_cs_index': trim_set['tail_cs_index'],
                            'thrust_nodes': trim_set['nodes_app'],
                            'initial_thrust': trim_set['initial_thrust'],
                            'fz_tolerance': .1,
                            'fx_tolerance': .1,
                            'm_tolerance': .1,
                            'max_iter': 250}
    
    settings['NonLinearDynamicCoupledStep'] = {'print_info': 'off',
                                            'max_iterations': 100,
                                            'delta_curved': 1e-1,
                                            'min_delta': tolerance,
                                            'newmark_damp': 5e-3,
                                            'gravity_on': gravity,
                                            'gravity': 9.81*g_loading,
                                            'num_steps': n_tstep,
                                            'dt': dt,
                                            'initial_velocity': u_inf}

    settings['NonLinearDynamicPrescribedStep'] = {'print_info': 'off',
                                        'max_iterations': 100,
                                        'delta_curved': 1e-1,
                                        'min_delta': tolerance,
                                        'newmark_damp': 5e-3,
                                        'gravity_on': gravity,
                                        'gravity': 9.81*g_loading,
                                        'num_steps': n_tstep,
                                        'dt': dt,
                                        'initial_velocity': u_inf*int(free_flight)}

    relative_motion = 'off'
    if not free_flight:
        relative_motion = 'on'

    settings['StepUvlm'] = {'print_info': 'off',
                            'num_cores': num_cores,
                            'convection_scheme': 2,
                            'gamma_dot_filtering': 3,
                            'velocity_field_generator': 'GustVelocityField',
                            'velocity_field_input': {'u_inf': int(not free_flight) * u_inf,
                                                     'u_inf_direction': [1., 0, 0],
                                                     'gust_shape': '1-cos',
                                                     'gust_parameters': {'gust_length': gust_length,
                                                                         'gust_intensity': gust_intensity},
                                                     'offset': gust_offset,
                                                     'relative_motion': relative_motion},
                            'rho': rho,
                            'n_time_steps': n_tstep,
                            'dt': dt}

    if free_flight:
        solver = 'NonLinearDynamicCoupledStep'
    else:
        solver = 'NonLinearDynamicPrescribedStep'

    settings['DynamicCoupled'] = {'structural_solver': solver,
                                'structural_solver_settings': settings[solver],
                                'aero_solver': 'StepUvlm',
                                'aero_solver_settings': settings['StepUvlm'],
                                'fsi_substeps': 200,
                                'fsi_tolerance': fsi_tolerance,
                                'relaxation_factor': relaxation_factor,
                                'minimum_steps': 1,
                                'relaxation_steps': 150,
                                'final_relaxation_factor': 0.5,
                                'n_time_steps': n_tstep,
                                'dt': dt,
                                'include_unsteady_force_contribution': 'on', # originally on
                                'postprocessors': ['BeamLoads', 'BeamPlot', 'AerogridPlot'],
                                'postprocessors_settings': {'BeamLoads': {'csv_output': 'off'},
                                                            'BeamPlot': {'include_rbm': 'on',
                                                                        'include_applied_forces': 'on'},
                                                            'AerogridPlot': {
                                                                'include_rbm': 'on',
                                                                'include_applied_forces': 'on',
                                                                'minus_m_star': 0},
                                                            }}

    settings['BeamLoads'] = {'csv_output': 'off'}

    settings['BeamPlot'] = {'include_rbm': 'on',
                            'include_applied_forces': 'on'}
                        

    settings['AerogridPlot'] = {'include_rbm': 'on',
                                'include_forward_motion': 'off',
                                'include_applied_forces': 'on',
                                'minus_m_star': 0,
                                'u_inf': u_inf,
                                'dt': dt}

    settings['Modal'] = {'print_info': True,
                    'use_undamped_modes': True,
                    'NumLambda': 30,
                    'rigid_body_modes': True,
                    'write_modes_vtk': 'on',
                    'print_matrices': 'on',
                    #  'write_data': 'on',
                    'continuous_eigenvalues': 'off',
                    'dt': dt,
                    'plot_eigenvalues': True}

    settings['LinearAssembler'] = {'linear_system': 'LinearAeroelastic',
                                    'linear_system_settings': {
                                        'beam_settings': {'modal_projection': False,
                                                        'inout_coords': 'nodes',
                                                        'discrete_time': True,
                                                        'newmark_damp': 0.05,
                                                        'discr_method': 'newmark',
                                                        'dt': dt,
                                                        'proj_modes': 'undamped',
                                                        'use_euler': 'off',
                                                        'num_modes': 40,
                                                        'print_info': 'on',
                                                        'gravity': 'on',
                                                        'remove_dofs': []},
                                        'aero_settings': {'dt': dt,
                                                        'integr_order': 2,
                                                        'density': rho,
                                                        'remove_predictor': False,
                                                        'use_sparse': True,
                                                        'rigid_body_motion': free_flight,
                                                        'use_euler': False,
                                                        'remove_inputs': ['u_gust']},
                                        'rigid_body_motion': free_flight}}

    settings['AsymptoticStability'] = {'sys_id': 'LinearAeroelastic',
                                        'print_info': 'on',
                                        'modes_to_plot': [],
                                        'display_root_locus': 'off',
                                        'frequency_cutoff': 0,
                                        'export_eigenvalues': 'off',
                                        'num_evals': 40,}
                                    

    config = configobj.ConfigObj()
    config.filename = file_name
    for k, v in settings.items():
        config[k] = v
    config.write()

generate_solver_file()

def opt(guess):
    
    global logs
    guess = guess / mag
    guess_final = list(np.array(guess))
    
    f = interp1d(guess_x_loc, guess_final, kind= 'linear')
    
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.h5_file_read_and_chord_vals import h5_file_read
    aero, fem = h5_file_read('OSU_Contribution/' + case_name)

    from OSU_Contribution.utils.Dyn_Min_Spar_Func.h5_file_read_and_chord_vals import find_chord_vals
    chord_vals = find_chord_vals(aero, trim_set['beam_opt_elems'], trim_set['beam_opt_sym'])

#---------------------------------------------------------------------------------------------------------------------
#--------------UPDATING STIFFNESS & MASS PROPERTIES--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------
    
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.stiffness_and_mass import scaled_SI2_OG
    ei_final, gj_final, mass_final, spar_geom = scaled_SI2_OG(trim_set['beam_opt_nodes'], fem, spar_cap_dict, f, chord_vals)
  
    # adding final stiffness and mass values to input dictionary 
    for i in range(len(trim_set['beam_opt_elems'])):
        fem['mass_db'][i][0][0] =  mass_final[i]
        fem['mass_db'][i][1][1] =  mass_final[i]
        fem['mass_db'][i][2][2] =  mass_final[i]
        
        fem['stiffness_db'][i][3][3] =  gj_final[i]
        fem['stiffness_db'][i][4][4] =  ei_final[i]

    # adding symmetric values
    if trim_set['beam_opt_sym'] == True:
        for i in range(len(trim_set['beam_opt_elems'])):
            fem['mass_db'][trim_set['beam_opt_elems'][-1] + i + 1][0][0] =  mass_final[i]
            fem['mass_db'][trim_set['beam_opt_elems'][-1] + i + 1][1][1] =  mass_final[i]
            fem['mass_db'][trim_set['beam_opt_elems'][-1] + i + 1][2][2] =  mass_final[i] 
            
            fem['stiffness_db'][trim_set['beam_opt_elems'][-1] + i + 1][3][3] =  gj_final[i]
            fem['stiffness_db'][trim_set['beam_opt_elems'][-1] + i + 1][4][4] =  ei_final[i]
                
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
    g_load =one_g_trim(route, case_name, logs, trim_set)
    
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
    aero, fem = turbulent_gust(route, case_name, logs, trim_set, fem)
    
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
#------------------------- MS CALCULATIONS--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------
    
    from OSU_Contribution.utils.Dyn_Min_Spar_Func.margin_of_safety import ms_OG_DYN
    m_s_fcu, m_s_ftu, m_s_su = ms_OG_DYN(gust_calc, spar_cap_dict, trim_set, spar_geom, chord_vals)
    
    # finding the minimum margin of safety for along the span for all time steps
    m_s_su_final = []
    m_s_fcu_final = []
    m_s_ftu_final = []
    
    for i in range(len(m_s_su[0])):
        m_s_su_final.append(min(m_s_su[:, i]))
        m_s_fcu_final.append(min(m_s_fcu[:, i]))
        m_s_ftu_final.append(min(m_s_ftu[:, i]))

#---------------------------------------------------------------------------------------------------------------------
#-------------------------PENALTIES--------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------
        
    penalties = 0
    # penalties for being under a margin or safety of 1 
    for i in range(len(m_s_ftu_final)): 
        penalties += max(0, (1 - m_s_ftu_final[i])/m_s_ftu_final[i])  
        penalties += max(0, (1 - m_s_fcu_final[i])/m_s_fcu_final[i])    
        penalties += max(0, (1 - m_s_su_final[i])/m_s_su_final[i])  
    
    max_deflection_val = gust_calc.structure.elements[trim_set['beam_opt_elems'][-1]].coordinates_ini[1][2] + max_deflection
    
    # max deflection across the gust encounter
    deflection = []
    for i in range(gust_calc.ts):
        deflection.append(gust_calc.structure.timestep_info[i].pos[trim_set['beam_opt_nodes'][-1]][2])
    
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
    
    logs['ms_tu_log'].append(m_s_ftu_final)
    logs['ms_su_log'].append(m_s_su_final) 
    logs['ms_cu_log'].append(m_s_fcu_final)  
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
# results = minimize(opt, guess, bounds = bnds)
# print(results)

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

opt(guess)


