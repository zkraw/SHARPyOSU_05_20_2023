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

# import sharpy.sharpy_main as sharpy_main
import sharpy.utils.algebra as algebra
from OSU_Contribution.Main_Paramaterized import Paramaterized
import matplotlib.pyplot as plt
import sharpy.sharpy_main as sharpy_main
from OSU_Contribution.utils.SHARPy_Input_Writes import generate_dyn_file
from OSU_Contribution.utils.SHARPy_Input_Writes import clean_test_files, write_aero, write_fem
from OSU_Contribution.utils.Rotation_Matrices import local_for_forces, nodal_a_for_2_b_for 


home = str(Path.home())
route = os.path.dirname(os.path.realpath(__file__)) 

# -----------------------------------------------------------------------------------------------------
#----------------------INPUTS--------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
file = 'Scaled_Solar_Impulse_Updated_Debug.xlsx'
global case_name
case_name = 'Initial Velocity Set'
excel = True 

input_file = os.path.dirname(os.path.abspath(__file__)) + '/Input Files/' + file

flow = ['BeamLoader',
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

# FLIGHT CONDITIONS for AERO SOLVER
# the simulation is set such that the aircraft flies at a u_inf velocity while
# the air is calm.
u_inf = 23.7 #14 for XHALE, 10 for T_Tail
rho = 0.1736 # 1.225 

# default elevator is -0.020708598 radians 3
alpha = 4.2*np.pi/180 # 1G 2.0251 alpha, .0533 N thrust, 3.7940 degrees elevator / 2G 4.287 degrees alpha, .2173 N thrust, 6.9285 degrees elevator 
beta = 0
roll = 0
gravity = 'on'
g_loading = 1 # acceleration of gravity (9.81) * factor 

# gust settings
gust_length = 1*u_inf
gust_offset = 0.5*u_inf
gust_intensity = 0.20

#DYNAMIC FLIGHT SETTINGS
free_flight = True
amplitude = 0
period = 0
if not free_flight:
    case_name += '_prescribed'
    amplitude = 0*np.pi/180
    period = 3
    case_name += '_amp_' + str(amplitude).replace('.', '') + '_period_' + str(period)
    
m = 5
physical_time = .1
tstep_factor = .2
# dt = 4./m/u_inf*tstep_factor
dt = physical_time
n_tstep = round(physical_time/dt)

span_main = 6

num_cores= 4
# n_step = 100

#numerics
n_step = 5
fsi_tolerance = 1e-4
structural_relaxation_factor = 0.3

relaxation_factor = 0.5
tolerance = 1e-2

num_cores = 4

# saving values to new file
#Paramaterized.excel_write(file, 'new_file.xlsx', ['start_Y'], [[3, 4, 5, 6, 7, 8]], [[123456789, 100, 100, 100, 100, 100]])

# XFIOL INPUTS

# speed of sound
a = 294.9        #speed of sound (m/s)
u = 1.4508e-05   #dynamic viscosity (kg/(m*s))

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

# generate FEM and Aero from input file
run = Paramaterized(case_name, input_file, xfoil, excel)


# Static Trim Settings
nodes_app = run.app_nodes # where the thrust nodes are 
cs_deflection = 2.9402*np.pi/180 # initial control surface deflection 1G .06621779 / 2G 0.1209251 
tail_cs_index = run.tail_index # the control surface index that will be used 
thrust = 1.332 # 1G .0533 N  / 2G .2173 N  

# OPTIMIZATION SETTINGS
global beam_opt_nodes
global beam_opt_sym

save_data = True
beam_opt_nodes = run.beam_opt_nodes
beam_opt_elems = run.beam_opt_elems
beam_opt_sym = run.beam_opt_sys
beam_save_file = 'Solar_Impulse_Min.xlsx'

# Span Mass Spar Creation
global span, surface_area, mass_desired, guess_x_loc, max_deflection, mag
surface_area = 162
span = 58.5
mass_desired = 890

max_deflection = np.tan(15*np.pi/180)*span / 2

# guess = [0.001]*6 # spar thickness
guess = [0.005]*6
guess_x_loc = np.linspace(0, 1, 6)  
max_thick = .05
min_thick = .0005
mag = 1000

x_taper_start = span * 0.686 / 2

global spar_cap_dict

spar_cap_dict = {}
spar_cap_dict['E_cap'] = 1.783e11
spar_cap_dict['E_lam'] = 1.212e10
spar_cap_dict['G'] = 4.51e10

h_spar_1 = 0.138
b_spar_1 = 0.25
h_spar_2 = 0.141
b_spar_2 = 0.188

spar_cap_dict['density'] = 1540 #163.934

# margin of safety calculations
spar_cap_dict['F_tu'] = 1.283e9
spar_cap_dict['F_cu'] = 4.61e8
spar_cap_dict['F_su'] = 2.35e8


global guess_log, diff_sum_out, masss, wing_tip_deflection_log, max_deflection_val, wing_deflection_diff_log, power_log, chord_log, penalties_log
guess_log = []
return_log = []
masss = []
max_strain_log = []
wing_tip_deflection_log = []
max_deflection_val_log = []
wing_deflection_diff_log = []
power_log = []
chord_log = []
penalties_log = []
global alpha_log, elevator_log, ei_log, gj_log, mass_log
alpha_log = []
elevator_log = []
ei_log = []
gj_log = []
mass_log = []
global ms_tu_log, ms_su_log, ms_cu_log
ms_tu_log = []
ms_su_log = []
ms_cu_log = [] 

global thrust_log
thrust_log = []

# -----------------------------------------------------------------------------------------------------
#----------------------End of Input--------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
Mach = u_inf/a

guess = np.array(guess)*mag
# bounds on minimization function
bnds = []
for i in range(len(guess_x_loc)):
    bnds.append(np.array([min_thick, max_thick])*mag)

# checking output file write
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

generate_dyn_file([dt, n_tstep, route, case_name, amplitude, period, free_flight], run.num_nodes)

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
                            'orientation': algebra.euler2quat(np.array([roll,
                                                                        alpha,
                                                                        beta]))}
    settings['AerogridLoader'] = {'unsteady': 'on',
                                'aligned_grid': 'on',
                                'mstar': int(10/tstep_factor),
                                'freestream_dir': ['1', '0', '0']}

    settings['NonLinearStatic'] = {'print_info': 'off',
                                'max_iterations': 100,
                                'num_load_steps': 1,
                                'delta_curved': 1e-1,
                                'min_delta': tolerance,
                                'gravity_on': gravity,
                                'gravity': 9.81*g_loading}

    settings['StaticUvlm'] = {'print_info': 'off',
                            'horseshoe': 'on',
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
                                'relaxation_factor': structural_relaxation_factor,
                                'correct_forces_method': 'PolarCorrection',
                                'correct_forces_settings': {'cd_from_cl': 'off',  # recommended settings (default)
                                                                'correct_lift': 'off',
                                                                'moment_from_polar': 'off'}}




    settings['StaticTrim'] = {'solver': 'StaticCoupled',
                            'solver_settings': settings['StaticCoupled'],
                            'initial_alpha': alpha,
                            'initial_deflection': cs_deflection,
                            'tail_cs_index': tail_cs_index,
                            'thrust_nodes': nodes_app,
                            'initial_thrust': thrust}

    
    settings['NonLinearDynamicCoupledStep'] = {'print_info': 'off',
                                            'max_iterations': 500,
                                            'delta_curved': 1e-1,
                                            'min_delta': tolerance,
                                            'newmark_damp': 5e-3,
                                            'gravity_on': gravity,
                                            'gravity': 9.81*g_loading,
                                            'num_steps': n_tstep,
                                            'dt': dt,
                                            'initial_velocity': 0}

    settings['NonLinearDynamicPrescribedStep'] = {'print_info': 'off',
                                        'max_iterations': 950,
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
                            #'horseshoe': 'on',
                            'num_cores': num_cores,
                            #'n_rollup': 0,
                            'convection_scheme': 2,
                            #'rollup_dt': dt,
                            #'rollup_aic_refresh': 1,
                            #'rollup_tolerance': 1e-4,
                            'gamma_dot_filtering': 6,
                            # 'velocity_field_generator': 'GustVelocityField',
                            # 'velocity_field_input': {'u_inf': int(not free_flight)*u_inf,
                            #                          'u_inf_direction': [1., 0, 0],
                            #                          'gust_shape': '1-cos',
                            #                          'gust_length': gust_length,
                            #                          'gust_intensity': gust_intensity*u_inf,
                            #                          'offset': gust_offset,
                            #                          'span': span_main,
                            #                          'relative_motion': relative_motion},
                            'velocity_field_generator': 'SteadyVelocityField',
                            'velocity_field_input': {'u_inf': u_inf,
                                                    'u_inf_direction': [1., 0, 0]},
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
                                'fsi_substeps': 100,
                                'fsi_tolerance': fsi_tolerance,
                                'relaxation_factor': relaxation_factor,
                                'minimum_steps': 1,
                                'relaxation_steps': 150,
                                'final_relaxation_factor': 0.5,
                                'n_time_steps': n_tstep,
                                'dt': dt,
                                'include_unsteady_force_contribution': 'on',
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
    
    guess = guess / mag

    # getting aero data
    aero = {}
    with h5.File('OSU_Contribution/' + case_name + '.aero.h5', 'r') as hdf:
        ls = list(hdf.keys())
        for i in ls:
            if i == 'airfoils':
                air = list(hdf.get(i).items())
                final = []
                for j in range(len(air)):
                    temp = np.array(air[j][1].value)
                    final.append(temp)
                aero[i] = final
                
            elif i == 'polars':
                pol = list(hdf.get(i).items())
                final_pol = []
                for j in range(len(pol)):
                    temp = np.array(pol[j][1].value)
                    final_pol.append(temp)
                aero[i] = final_pol
            
            elif i == 'm_distribution':
                aero[i] = 'uniform'
                
            else: 
                aero[i] = np.array(hdf.get(i)) 

    # getting fem data
    fem = {}
    with h5.File('OSU_Contribution/' + case_name + '.fem.h5', 'r') as hdf:
        ls = list(hdf.keys())
        for i in ls:
            fem[i] = np.array(hdf.get(i)) 
    
    # scaling wingspan 
    span_init =  fem['coordinates'][beam_opt_nodes[-1]][1]
    scale = span / (2*span_init)
    
    # COORDINATES update
    for i in range(len(beam_opt_nodes)):
        fem['coordinates'][i][1] = scale * fem['coordinates'][i][1]
    
    if beam_opt_sym == True:
       for i in range(len(beam_opt_nodes)):
           fem['coordinates'][beam_opt_nodes[-1] + i + 1][1] =  scale * fem['coordinates'][beam_opt_nodes[-1] + i + 1][1]
    
    # CHORD update
    chord_vals = []
    for i in range(len(beam_opt_elems)):
        aero['chord'][i] = (1/scale) * aero['chord'][i]
        for j in range(3):
            chord_vals.append(aero['chord'][i][j])
    
    if beam_opt_sym == True:
        for i in range(len(beam_opt_elems)):
            aero['chord'][beam_opt_elems[-1] + i + 1] = (1/scale) * aero['chord'][beam_opt_elems[-1] + i + 1]
            for j in range(3):
                chord_vals.append(aero['chord'][beam_opt_elems[-1] + i + 1][j])
 
    
    guess_final = list(np.array(guess))
       
    f = interp1d(guess_x_loc, guess_final, kind= 'linear')

    max_span =  fem['coordinates'][beam_opt_nodes[-1]][1] 
    
    thickness = []
    global span_loc_log
    span_loc_log = []
    h_spar = []
    b_spar = []
    # finding ei along the span
    for i in range(len(beam_opt_nodes)):
        span_loc_log.append(fem['coordinates'][i][1])
        thickness.append(f(fem['coordinates'][i][1]/max_span))
        if fem['coordinates'][i][1] <= x_taper_start:
            h_spar.append((1/scale) * h_spar_1)
            b_spar.append((1/scale) * b_spar_1)
        else:
            h_spar.append((1/scale)*((h_spar_2 - h_spar_1)*(fem['coordinates'][i][1] - x_taper_start)/(span/2 - x_taper_start) + h_spar_1))
            b_spar.append((1/scale)*((b_spar_2 - b_spar_1)*(fem['coordinates'][i][1] - x_taper_start)/(span/2 - x_taper_start) + b_spar_1))
        
    
    # stiffness and mass value updates
    ei_span = []
    gj_span = []
    m_span = []
    # torsion ms values
    a_cross = []
    a_fwdaft = []
    a_topbot = []
    # constant values of the geometry 
    for i in range(len(thickness)):
        thick = thickness[i]
        # new width values 
        b_cap = .15*b_spar[i]
        b_lam_up_low = b_spar[i] - 2 * b_cap
       
        # EI Calculation
        t_fwdaft = thick
        h_fwdaft = h_spar[i] - 2 * thick
        EI_fwdaft = spar_cap_dict['E_lam'] * t_fwdaft * h_fwdaft ** 3 / 12
        a_fwdaft.append(t_fwdaft*h_fwdaft)
        a_topbot.append(thick * b_lam_up_low)
        H_couple = (h_spar[i]) 
        EI_topbot = spar_cap_dict['E_lam'] * a_topbot[-1] * (H_couple / 2) ** 2
        a_cap = thick * b_cap
        EI_cap = spar_cap_dict['E_cap'] * a_cap * (H_couple / 2) ** 2
        ei_span.append(2 * EI_topbot + 2 * EI_fwdaft + 4 * EI_cap)
        # GJ calculation
        s = 2 * (h_spar[i] - 2*thick) + 2 * (b_lam_up_low)
        a_cross.append(h_spar[i] * b_spar[i])
        gj_span.append(4 * a_cross[-1] ** 2 / (s / (spar_cap_dict['G'] * thick)))
        # m_bar calculation
        m_span.append((spar_cap_dict['density']*(4*a_cap + 2*a_topbot[-1] + 2*a_fwdaft[-1])))
 
    ei_final = []
    #taking the average of the nodal ei
    for i in range(0, len(ei_span) - 2, 2):
        ei_final.append(abs((ei_span[i] + ei_span[i + 1] + ei_span[i + 2]) / 3)) 
        
      # adding final stiffness values to input dictionary 
    for i in range(len(beam_opt_elems)):
        fem['stiffness_db'][i][4][4] =  ei_final[i]

    # adding symmetric values
    if beam_opt_sym == True:
        for i in range(len(beam_opt_elems)):
            fem['stiffness_db'][beam_opt_elems[-1] + i + 1][4][4] =  ei_final[i]

    gj_final = []
    #taking the average of the nodal gj
    for i in range(0, len(gj_span) - 2, 2):
        gj_final.append(abs((gj_span[i] + gj_span[i + 1] + gj_span[i + 2]) / 3))
        
    # adding final stiffness values to input dictionary 
    for i in range(len(beam_opt_elems)):
        fem['stiffness_db'][i][3][3] =  gj_final[i]

    # adding symmetric values
    if beam_opt_sym == True:
        for i in range(len(beam_opt_elems)):
            fem['stiffness_db'][beam_opt_elems[-1] + i + 1][3][3] =  gj_final[i]
            
    # updating mass values 
    mass_final = []
    # mass from solar panels and other components of wing
    for i in range(0, len(m_span)-2, 2):
        m_solar = 0.437 * chord_vals[i]
        m_lete = .1043 * chord_vals[i]**2 * 3.722 
        mass_final.append((m_span[i] + m_span[i+1] + m_span[i+2])/3 + m_solar + m_lete)

    for i in range(len(beam_opt_elems)):
        fem['mass_db'][i][0][0] =  mass_final[i]
        fem['mass_db'][i][1][1] =  mass_final[i]
        fem['mass_db'][i][2][2] =  mass_final[i]

    if beam_opt_sym == True:
        for i in range(len(beam_opt_elems)):
            fem['mass_db'][beam_opt_elems[-1] + i + 1][0][0] =  mass_final[i]
            fem['mass_db'][beam_opt_elems[-1] + i + 1][1][1] =  mass_final[i]
            fem['mass_db'][beam_opt_elems[-1] + i + 1][2][2] =  mass_final[i] 
        
    #deletes old input files
    clean_test_files(case_name)
    # creates new aero file with updated twist
    write_fem(fem, case_name)
    # creates new aero file with updated twist
    write_aero(aero, case_name)
    
    # 2G max loading condition 
    # modifiying input to have 2G condition
    file = open(route + '/' + case_name + '.sharpy')
    lines = file.readlines()
    nu_lines = []
    for i in lines:
        current = i.replace(' ', '').split('=')
        
        if current[0] == 'gravity':
            current[1] = str(1*9.81) + '\n'
            nu_lines.append(' = '.join(current))
            continue

        if current[0] == 'initial_alpha': 
            if alpha_log != []:
                current[1] = str(alpha_log[-1]) + '\n'
                nu_lines.append(' = '.join(current))
            continue

        if current[0] == 'initial_deflection':
            if elevator_log != []:
                current[1] = str(elevator_log[-1]) + '\n'
                nu_lines.append(' = '.join(current))
            continue

        if current[0] == 'initial_thrust':
            current[1] = '2.3116\n'
            nu_lines.append(' = '.join(current))
            continue
        
        nu_lines.append(i) 

    os.remove(route + '/' + case_name + '.sharpy')
    
    with open(route + '/' + case_name + '.sharpy', 'w') as f:
        for line in nu_lines:
            f.write(line)

    # 2G Static Trim case
    data2 = sharpy_main.main(['', route + '/' + case_name + '.sharpy'])
    
    if data2.structure.div_flag.value != 1:
        node_forces, local_applied_forces_log = local_for_forces(data2) 
        
        nodes = beam_opt_nodes
        shear = np.zeros(len(nodes))
        for i in range(len(nodes)-1, -1, -1):
            if i == len(nodes)-1:
                shear[i] = node_forces[i][2]
            elif i == 0:
                shear[i] = .5 * node_forces[i][2] + shear[i+1]
            else:
                shear[i] = node_forces[i][2] + shear[i+1]

        moment =  np.zeros(len(nodes))
        for i in range(len(nodes) - 1, -1, -1):
            if i == len(nodes) - 1:
                moment[i] = 0
            else:
                back = data2.structure.timestep_info[-1].pos[i+1]
                current = data2.structure.timestep_info[-1].pos[i]
                moment[i] = shear[i+1] * (back[1] - current[1]) + moment[i+1]

        strain = np.zeros(len(nodes)) 
        elem_stiff = data2.structure.elem_stiffness
        stiff = data2.structure.stiffness_db
        elem_stiff_counter = int(0)
        m_s_ftu = []
        m_s_fcu = []
        
        # cap structural calculations
        for i in range(len(nodes)):
            strain[i] = abs(((moment[i] * H_couple * .5) / (stiff[int(elem_stiff[elem_stiff_counter])][4][4]))) # ei other taken out to find strain on spar
            if strain[i] != 0:
                m_s_ftu.append(spar_cap_dict['F_tu']/(spar_cap_dict['E_cap']*strain[i]) - 1)
                m_s_fcu.append(spar_cap_dict['F_cu']/(spar_cap_dict['E_cap']*strain[i]) - 1)
            else:
                m_s_ftu.append(1000)
                m_s_fcu.append(1000)
            if ( i % 2) == 0 and i != 0:
                elem_stiff_counter +=1
        
        m_s_su = []
        fwd_shear_log = []
        aft_shear_log = []
        bending_shear_log = []
        msfwd_log = []
        msaft_log = []
        
        elastic_axis = []
        row_count = 0
        order_count = 0
        
        # laminant shear stress calculations
        for i in range(len(nodes)):
            if i == 0:
                order = [0, 2, 1]
            else:
                order = [2, 1]    
            if len(order) == order_count:
                row_count += 1
                order_count = 0
            elastic_axis.append(data2.aero.aero_dict['elastic_axis'][row_count][order[order_count]])
            order_count += 1
            
        torsion = []
        for i in range(len(nodes)-1, -1, -1):
            if i == len(nodes)-1:
                torsion.append(local_applied_forces_log[i][2] * (chord_vals[i]*(elastic_axis[i]-.25)))
            else:
                torsion.append(local_applied_forces_log[i][2] * (chord_vals[i]*(elastic_axis[i]-.25)) + torsion[-1])
        
        torsion.reverse()
        for i in range(len(nodes)):
            # torsional shear    
            shear_flow_torsion = torsion[i] / (2*a_cross[i])
            
            #up_lwr_shear = shear_flow_torsion / self.up_low_thick[i]
            fwd_shear_tors = shear_flow_torsion / t_fwdaft
            aft_shear_tors = -shear_flow_torsion / t_fwdaft

            # bending shear
            if i == 0:
                bending_shear = .5 * shear[i] / a_fwdaft[i]
            else:
                bending_shear = shear[i] / a_fwdaft[i]

            fwd_shear_log.append(fwd_shear_tors)
            aft_shear_log.append(aft_shear_tors)
            bending_shear_log.append(bending_shear)
            
            fwd_shear = fwd_shear_tors + bending_shear
            aft_shear = aft_shear_tors + bending_shear
           
            M_S_fwd = (spar_cap_dict['F_su'] / abs(fwd_shear)) - 1
            M_S_aft = (spar_cap_dict['F_su'] / abs(aft_shear)) - 1
            m_s_su.append(min([M_S_fwd, M_S_aft]))
            msfwd_log.append(M_S_fwd)
            msaft_log.append(M_S_aft)
            
        alpha_log.append(data2.alpha)
        elevator_log.append(data2.elevator)
        
    # if diverged, has inputs for required variables
    else:
        m_s_fcu = [.05]
        m_s_ftu = [.05]
        m_s_su = [.05]
        alpha_log.append(0)
        elevator_log.append(0)
    
    max_deflection_val = data2.structure.elements[beam_opt_elems[-1]].coordinates_ini[1][2] + max_deflection
    deflection = data2.structure.timestep_info[-1].pos[beam_opt_nodes[-1]][2]  
    
    penalties = 0
    # penalties for being under a margin or safety of 1 
    for i in range(len(m_s_ftu)): 
        penalties += max(0, (1 - m_s_ftu[i])/m_s_ftu[i])  
        penalties += max(0, (1 - m_s_fcu[i])/m_s_fcu[i])    
        penalties += max(0, (1 - m_s_su[i])/m_s_su[i])  
      
    # penalties for being over max wingtip deflection 
    penalties += max(0, 50*(deflection-max_deflection_val)/max_deflection_val)**2 
    
    # penalties if thickness is not decreasing
    for i in range(1, len(guess)):
        penalties += max(0, (guess[i] - guess[i-1])/guess[i-1]) 

    # calculating the mass of the configuration
    grav_forces = data2.structure.timestep_info[0].gravity_forces
    grav = 0
    for i in range(len(grav_forces)):
        # grav += grav_forces[i][2]
        grav += np.sqrt(grav_forces[i][0]**2 + grav_forces[i][1]**2 + grav_forces[i][2]**2)
    
    masss_temp = abs(grav/(1*9.81)) 
     
    cost = (masss_temp / mass_desired)
            
    guess_log.append(guess)
    masss.append(masss_temp)
    return_log.append(cost + penalties)
    penalties_log.append(penalties)
    ms_tu_log.append(m_s_ftu)
    ms_su_log.append(m_s_su) 
    ms_cu_log.append(m_s_fcu)  
    chord_log.append(chord_vals[0])
    wing_tip_deflection_log.append(data2.structure.timestep_info[-1].pos[beam_opt_nodes][-1][2])
    max_deflection_val_log.append(max_deflection_val)
    thrust_log.append([data2.thrust*u_inf])     
    df = pd.DataFrame([guess_log, masss, return_log, penalties_log, ms_tu_log, ms_su_log, ms_cu_log, chord_log, alpha_log, elevator_log, wing_tip_deflection_log, max_deflection_val_log, thrust_log]) #, ei_log, gj_log, mass_log])
    df.to_excel( beam_save_file, header= False,  sheet_name= 'Minimization_Logs')
        
    return cost + penalties

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


