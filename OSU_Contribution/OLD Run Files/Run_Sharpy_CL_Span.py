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
import h5py
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
from Main_Paramaterized import Paramaterized
import matplotlib.pyplot as plt
import sharpy.sharpy_main as sharpy_main

home = str(Path.home())
route = os.path.dirname(os.path.realpath(__file__)) 

# -----------------------------------------------------------------------------------------------------
#----------------------INPUTS--------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
file = 'T_Tail_Flat2.xlsx'
global case_name
case_name = 'test_Beam'
excel = True 

input_file = os.path.dirname(os.path.abspath(__file__)) + '/' + file

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

free_flight = True
if not free_flight:
    case_name += '_prescribed'
    amplitude = 0*np.pi/180
    period = 3
    case_name += '_amp_' + str(amplitude).replace('.', '') + '_period_' + str(period)


# FLIGHT CONDITIONS for AERO SOLVER
# the simulation is set such that the aircraft flies at a u_inf velocity while
# the air is calm.
global u_inf, rho
u_inf = 12 #14 for XHALE, 10 for T_Tail
rho = 1.225 # 1.225 

# default elevator is -0.020708598 radians 3
alpha = 5*np.pi/180 # 1G 2.0251 alpha, .0533 N thrust, 3.7940 degrees elevator / 2G 4.287 degrees alpha, .2173 N thrust, 6.9285 degrees elevator 
beta = 0
roll = 0
gravity = 'on'
global g_loading
g_loading = 1 # acceleration of gravity (9.81) * factor 

# gust settings
gust_length = 1*u_inf
gust_offset = 0.5*u_inf
gust_intensity = 0.20

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

# temporal discretisation
m = 5
physical_time = .1
tstep_factor = .2
# dt = 4./m/u_inf*tstep_factor
dt = physical_time
n_tstep = round(physical_time/dt)

# saving values to new file
#Paramaterized.excel_write(file, 'new_file.xlsx', ['start_Y'], [[3, 4, 5, 6, 7, 8]], [[123456789, 100, 100, 100, 100, 100]])

# generate FEM and Aero from input file
run = Paramaterized(case_name, input_file, excel)


# Static Trim Settings
nodes_app = run.app_nodes # where the thrust nodes are 
cs_deflection = -3.51*np.pi/180 # initial control surface deflection 1G .06621779 / 2G 0.1209251 
tail_cs_index = run.tail_index # the control surface index that will be used 
thrust = 10.2102   

# beam Opt settings
global beam_opt_nodes
global beam_opt_sym

save_data = True
beam_opt_nodes = run.beam_opt_nodes
beam_opt_elems = run.beam_opt_elems
beam_opt_sym = run.beam_opt_sys

global max_strain, beam_save_file
global min_thick
global bb
global h, e, g

max_strain = .002
max_deflection = 5 # 4 meters in z deflection
beam_save_file = '32m_Stiff.xlsx'
min_thick = .0001
max_thick = .003
low_twist = -10
high_twist = 10

bb =  .11999999999
h =  .13235760
e = 155e9
pois = .309 # poissons ratio
g = e / (2*(1+ pois))

global ei_other, gj_other, mass_other, density
ei_other = 6000
gj_other = 3000
mass_other = .375
# 70% intermediate carbon density and 30% epoxy 
density = 1790 * .7 + 1200 * .3

global guess_log, diff_sum_out, masss, wing_tip_deflection_log, max_deflection_val, wing_deflection_diff_log, power_log, chord_log, penalties_log
guess_log = []
diff_sum_out = []
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

global span, surface_area
surface_area = 32
span = 25
guess = [span]
max_span = 80
min_span = 10

starting_mass = 150
starting_power = 250

global spar_cap_mode, spar_cap_dict
spar_cap_mode = False

spar_cap_dict = {}
spar_cap_dict['E_cap'] = 1.783e11
spar_cap_dict['E_lam'] = 1.212e10
spar_cap_dict['G'] = 4.51e10

spar_cap_dict['h_spar'] = 0.129
spar_cap_dict['b_spar'] = 0.2

spar_cap_dict['h_cap'] = 0.012
spar_cap_dict['b_cap'] = 0.018
spar_cap_dict['density'] = 163.934

# XFOIL Inputs
global airfoil_name, Mach, u, n_iter, alfa_1, alfa_2, alfa_step
airfoil_name = 'naca0015'
a = 343         #speed of sound (m/s)
u = 1.802e-05   #dynamic viscosity (kg/(m*s))
n_iter = 100

# alfa steps for polars in XFOIL
alfa_1 = -16
alfa_2 = 16
alfa_step = .5

# -----------------------------------------------------------------------------------------------------
#----------------------End of Input--------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
Mach = u_inf/a

# bounds on minimization function
bnds = [[min_span, max_span]]

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

def generate_dyn_file():
    global dt
    global n_tstep
    global route
    global case_name
    global num_elem
    global num_node_elem
    num_node = run.num_nodes
    global amplitude
    global period
    global free_flight

    # num_node_main = 3

    dynamic_forces_time = None
    with_dynamic_forces = True
    with_forced_vel = False
    if not free_flight:
        with_forced_vel = True

    if with_dynamic_forces:
        f1 = 0
        dynamic_forces = np.zeros((num_node, 6))
        app_node = [int(16)]
        dynamic_forces[app_node, 2] = f1
        force_time = np.zeros((n_tstep, ))
        limit = round(0.05/dt)
        force_time[0:3] = 1

        dynamic_forces_time = np.zeros((n_tstep, num_node, 6))
        for it in range(n_tstep):
            dynamic_forces_time[it, :, :] = force_time[it]*dynamic_forces

    forced_for_vel = None
    if with_forced_vel:
        forced_for_vel = np.zeros((n_tstep, 6))
        forced_for_acc = np.zeros((n_tstep, 6))
        for it in range(n_tstep):
            # if dt*it < period:
            # forced_for_vel[it, 2] = 2*np.pi/period*amplitude*np.sin(2*np.pi*dt*it/period)
            # forced_for_acc[it, 2] = (2*np.pi/period)**2*amplitude*np.cos(2*np.pi*dt*it/period)

            forced_for_vel[it, 3] = 2*np.pi/period*amplitude*np.sin(2*np.pi*dt*it/period)
            forced_for_acc[it, 3] = (2*np.pi/period)**2*amplitude*np.cos(2*np.pi*dt*it/period)

    if with_dynamic_forces or with_forced_vel:
        with h5.File(route + '/' + case_name + '.dyn.h5', 'a') as h5file:
            if with_dynamic_forces:
                h5file.create_dataset(
                    'dynamic_forces', data=dynamic_forces_time)
            if with_forced_vel:
                h5file.create_dataset(
                    'for_vel', data=forced_for_vel)
                h5file.create_dataset(
                    'for_acc', data=forced_for_acc)
            h5file.create_dataset(
                'num_steps', data=n_tstep)


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

    settings['BeamOpt'] = {'solver': 'StaticCoupled',
                        'solver_settings': settings['StaticCoupled'],
                        'Max Strain': max_strain,
                        'Beam_Opt_Nodes': beam_opt_nodes,
                        'Beam_Opt_Elems' : beam_opt_elems,
                        'Sym': beam_opt_sym,
                        'initial_deflection': cs_deflection}

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
generate_dyn_file()

def clean_test_files():
    fem_file_name = 'Zack_Beam' + '/' + case_name + '.fem.h5'
    if os.path.isfile(fem_file_name):
        os.remove(fem_file_name)

    aero_file_name = 'Zack_Beam' + '/' + case_name + '.aero.h5'
    if os.path.isfile(aero_file_name):
        os.remove(aero_file_name)
        
def write_aero(aero):
    
    with h5.File('Zack_Beam' + '/' + case_name + '.aero.h5', 'a') as h5file:
        airfoils_group = h5file.create_group('airfoils')
        polars_group = h5file.create_group('polars')

        for i in range(len(aero['airfoils'])):  # len(air_data)
            airfoil = airfoils_group.create_dataset(str(i), data= aero['airfoils'][i])
            polars_input = polars_group.create_dataset(str(i), data=aero['polars'][i])
        
        # chord
        chord_input = h5file.create_dataset('chord', data= aero['chord'])
        dim_attr = chord_input.attrs['units'] = 'm'

        # twist
        twist_input = h5file.create_dataset('twist', data=aero['twist'])
        dim_attr = twist_input.attrs['units'] = 'rad'

        # sweep
        sweep_input = h5file.create_dataset('sweep', data=aero['sweep'])
        dim_attr = sweep_input.attrs['units'] = 'rad'

        # airfoil distribution
        airfoil_distribution_input = h5file.create_dataset('airfoil_distribution', data=aero['airfoil_distribution'])
        surface_distribution_input = h5file.create_dataset('surface_distribution', data=aero['surface_distribution'])
        surface_m_input = h5file.create_dataset('surface_m', data=aero['surface_m'])
        m_distribution_input = h5file.create_dataset('m_distribution', data=aero['m_distribution'].encode('ascii', 'ignore'))
        aero_node_input = h5file.create_dataset('aero_node', data=aero['aero_node'])
        elastic_axis_input = h5file.create_dataset('elastic_axis', data=aero['elastic_axis'])
        control_surface_input = h5file.create_dataset('control_surface', data=aero['control_surface'])
        control_surface_deflection_input = h5file.create_dataset('control_surface_deflection', data=aero['control_surface_deflection'])
        control_surface_chord_input = h5file.create_dataset('control_surface_chord', data=aero['control_surface_chord'])
        control_surface_hinge_coord_input = h5file.create_dataset('control_surface_hinge_coord', data=aero['control_surface_hinge_coord'])
        control_surface_types_input = h5file.create_dataset('control_surface_type', data=aero['control_surface_type'])

def write_fem(fem):
    
    with h5.File('Zack_Beam' + '/' + case_name + '.fem.h5', 'a') as h5file:
        coordinates = h5file.create_dataset('coordinates', data = fem['coordinates'])
        conectivities = h5file.create_dataset('connectivities', data = fem['connectivities'])
        num_nodes_elem_handle = h5file.create_dataset('num_node_elem', data = fem['num_node_elem'])
        num_nodes_handle = h5file.create_dataset('num_node', data = fem['num_node'])
        num_elem_handle = h5file.create_dataset('num_elem', data = fem['num_elem'])
        stiffness_db_handle = h5file.create_dataset('stiffness_db', data = fem['stiffness_db'])
        stiffness_handle = h5file.create_dataset('elem_stiffness', data = fem['elem_stiffness'])
        mass_db_handle = h5file.create_dataset('mass_db', data = fem['mass_db'])
        mass_handle = h5file.create_dataset('elem_mass', data = fem['elem_mass'])
        frame_of_reference_delta_handle = h5file.create_dataset('frame_of_reference_delta', data=fem['frame_of_reference_delta'])
        structural_twist_handle = h5file.create_dataset('structural_twist', data = fem['structural_twist'])
        bocos_handle = h5file.create_dataset('boundary_conditions', data=fem['boundary_conditions'])
        beam_handle = h5file.create_dataset('beam_number', data=fem['beam_number'])
        app_forces_handle = h5file.create_dataset('app_forces', data=fem['app_forces'])
        lumped_mass_nodes_handle = h5file.create_dataset('lumped_mass_nodes', data=fem['lumped_mass_nodes'])
        lumped_mass_handle = h5file.create_dataset('lumped_mass', data=fem['lumped_mass'])
        lumped_mass_inertia_handle = h5file.create_dataset('lumped_mass_inertia', data=fem['lumped_mass_inertia'])
        lumped_mass_position_handle = h5file.create_dataset('lumped_mass_position', data=fem['lumped_mass_position'])
        
# converts from aircraft frame of reference to local frame of reference
def nodal_a_for_2_b_for(data, nodal, tstep, filter=np.array([True] * 6)):
    nodal_a = nodal.copy(order='F')
    for i_node in range(data.structure.num_node):
        # get master elem and i_local_node
        i_master_elem, i_local_node = data.structure.node_master_elem[i_node, :]
        crv = tstep.psi[i_master_elem, i_local_node, :]
        cab = algebra.crv2rotation(crv)
        temp = np.zeros((6,))
        temp[0:3] = np.dot(cab, nodal[i_node, 0:3])
        temp[3:6] = np.dot(cab, nodal[i_node, 3:6])
        for i in range(6):
            if filter[i]:
                nodal_a[i_node, i] = temp[i]
    return nodal_a

def local_for_forces(data, tstep=None):
    if tstep is None:
        tstep = data.structure.timestep_info[data.ts]
    applied_forces = data.structure.nodal_b_for_2_a_for(tstep.steady_applied_forces,
                                                                tstep)

    applied_forces_copy = applied_forces.copy()
    global lift
    lift = applied_forces_copy
    #lift[0] = lift[0]*2 
    gravity_forces_copy = tstep.gravity_forces.copy()
    
    
    for i_node in range(data.structure.num_node):
        applied_forces_copy[i_node, 3:6] += np.cross(tstep.pos[i_node, :],
                                                        applied_forces_copy[i_node, 0:3])
        gravity_forces_copy[i_node, 3:6] += np.cross(tstep.pos[i_node, :],
                                                        gravity_forces_copy[i_node, 0:3])


    # forces are in frame of reference a
    together = np.zeros((len(applied_forces_copy), 6))
    for i in range(len(applied_forces_copy)):
        if i == 0:
            together[i, ...] = applied_forces_copy[i] + gravity_forces_copy[i]
        else:
            together[i, ...] = applied_forces_copy[i] + gravity_forces_copy[i]
            
    # # only need half at the beginning of span
    # applied_forces_copy[0] = applied_forces_copy[0] * .5
    # gravity_forces_copy[0] = gravity_forces_copy[0] * .5

    # for the logger 
    global local_grav_forces_log, local_applied_forces_log
    local_grav_forces_log = nodal_a_for_2_b_for(data, gravity_forces_copy, tstep)
    local_applied_forces_log = nodal_a_for_2_b_for(data, applied_forces_copy, tstep)
    
    return nodal_a_for_2_b_for(data, together, tstep)
    
def span_opt(guess):

    guess_log.append(guess[0])

    # getting aero data
    aero = {}
    with h5py.File('Zack_Beam/' + case_name + '.aero.h5', 'r') as hdf:
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
    with h5py.File('Zack_Beam/' + case_name + '.fem.h5', 'r') as hdf:
        ls = list(hdf.keys())
        for i in ls:
            fem[i] = np.array(hdf.get(i)) 
    
    span_init =  fem['coordinates'][beam_opt_nodes[-1]][1]
    scale = guess[0] / (2*span_init)
    
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

    chord_ave = sum(chord_vals) / len(chord_vals)
    
    # reynolds number calculation
    Re = (rho * u_inf * chord_ave) / u
    
    # Removing previous polar_file.txt file
    if os.path.exists("polar_file.txt"):
        os.remove("polar_file.txt")
    
    XFOIL_PATH = str(os.getcwd()) + '/Zack_Beam/XFOIL6.99/xfoil.exe'

    xfoil = subprocess.Popen(XFOIL_PATH, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         universal_newlines=True)
    
    actions = ["NACA 0015 \n ",
               #["LOAD {0}.dat \n ".format(airfoil_name),
                                      # "{0} \n ".format(airfoil_name),
                                       "PANE \n ",
                                       "OPER \n ",
                                       "Visc {0} \n ".format(Re),
                                       "Mach {0} \n ".format(Mach),
                                       "PACC \n ",
                                       "polar_file.txt \n\n ",
                                       "ITER {0} \n ".format(n_iter),
                                       "aseq {0} {1} {2} \n ". format(alfa_1, alfa_2, alfa_step),
                                       "quit"]
    command = ''.join(actions)
        
    xfoil.communicate(input=command)
    
    # Updating the drag polar for the change in chord
    nu_pol = np.loadtxt("polar_file.txt", skiprows=12)
    new_polar = []
    for i in range(len(nu_pol)):
        new_polar.append([nu_pol[i][0]*np.pi / 180, nu_pol[i][1], nu_pol[i][2], nu_pol[i][4]])
    
    pol_ref = aero['airfoil_distribution'][beam_opt_elems[0]][0]
    aero['polars'][pol_ref] = new_polar
    
    # Calculate m_bar to keep CL constant
    m_nu = (((.5*rho*u_inf**2*surface_area*CL) / 9.81) - 54.2) / guess[0]
    for i in range(len(beam_opt_elems)):
            fem['mass_db'][i][0][0] =  m_nu
            fem['mass_db'][i][1][1] =  m_nu
            fem['mass_db'][i][2][2] =  m_nu

    if beam_opt_sym == True:
        for i in range(len(beam_opt_elems)):
            fem['mass_db'][beam_opt_elems[-1] + i + 1][0][0] =  m_nu
            fem['mass_db'][beam_opt_elems[-1] + i + 1][1][1] =  m_nu
            fem['mass_db'][beam_opt_elems[-1] + i + 1][2][2] =  m_nu
    
    # CAP SPAR MODE
    if spar_cap_mode == True:
        # scaling the geometric values according to the change in chord
        spar_cap_dict['h_spar'] = (1/scale) * spar_cap_dict['h_spar']
        spar_cap_dict['b_spar'] = (1/scale) * spar_cap_dict['b_spar']
        spar_cap_dict['h_cap'] = (1/scale) * spar_cap_dict['h_cap']
        spar_cap_dict['b_cap'] = (1/scale) * spar_cap_dict['b_cap']
        
        t_fwdaft = spar_cap_dict['h_cap']
        h_fwdaft = spar_cap_dict['h_spar'] - 2 * spar_cap_dict['h_cap']
        
        EI_fwdaft = spar_cap_dict['E_lam'] * t_fwdaft * h_fwdaft ** 3 / 12

        A_topbot = spar_cap_dict['h_cap'] * (spar_cap_dict['b_spar'] - 2 * spar_cap_dict['b_cap'])
        H_couple = (spar_cap_dict['h_spar'] - spar_cap_dict['h_cap']) / 2

        EI_topbot = spar_cap_dict['E_lam'] * A_topbot * (H_couple / 2) ** 2

        A_cap = spar_cap_dict['h_cap'] * spar_cap_dict['b_cap']

        EI_cap = spar_cap_dict['E_cap'] * A_cap * (H_couple / 2) ** 2

        EI_y = 2 * EI_topbot + 2 * EI_fwdaft + 4 * EI_cap

        t_spar = spar_cap_dict['h_cap']
        s = 2 * (spar_cap_dict['h_spar'] - 2*spar_cap_dict['h_cap']) + 2 * (spar_cap_dict['b_spar'] - 2*spar_cap_dict['b_cap'])

        A_cross = spar_cap_dict['h_spar'] * spar_cap_dict['b_spar']

        GJ = 4 * A_cross ** 2 / (s / (spar_cap_dict['G'] * t_spar))
        
         # adding final stiffness values to input dictionary 
        for i in range(len(beam_opt_elems)):
            fem['stiffness_db'][i][4][4] =  EI_y

        # adding symmetric values
        if beam_opt_sym == True:
            for i in range(len(beam_opt_elems)):
                fem['stiffness_db'][beam_opt_elems[-1] + i + 1][4][4] =  EI_y
            
        # adding final stiffness values to input dictionary 
        for i in range(len(beam_opt_elems)):
            fem['stiffness_db'][i][3][3] =  GJ

        # adding symmetric values
        if beam_opt_sym == True:
            for i in range(len(beam_opt_elems)):
                fem['stiffness_db'][beam_opt_elems[-1] + i + 1][3][3] =  GJ

        m_nu = spar_cap_dict['density'] * (spar_cap_dict['h_spar']*spar_cap_dict['b_spar'] - ((spar_cap_dict['h_spar']-2*t_spar)*(spar_cap_dict['b_spar']-2*t_spar)))
        
        for i in range(len(beam_opt_elems)):
            fem['mass_db'][i][0][0] =  m_nu
            fem['mass_db'][i][1][1] =  m_nu
            fem['mass_db'][i][2][2] =  m_nu

        if beam_opt_sym == True:
            for i in range(len(beam_opt_elems)):
                fem['mass_db'][beam_opt_elems[-1] + i + 1][0][0] =  m_nu
                fem['mass_db'][beam_opt_elems[-1] + i + 1][1][1] =  m_nu
                fem['mass_db'][beam_opt_elems[-1] + i + 1][2][2] =  m_nu
        
    # deletes previous fem file 
    clean_test_files()

    # creates new aero file with updated twist
    write_fem(fem)
    # creates new aero file with updated twist
    write_aero(aero)
    
    # 1G static trim at cruise  
    # modifiying input to have 1G condition
    file = open(route + '/' + case_name + '.sharpy')
    lines = file.readlines()
    nu_lines = []
    for i in lines:
        current = i.replace(' ', '').split('=')
        
        if current[0] == 'gravity':
            current[1] = '9.81\n'
            nu_lines.append(' = '.join(current))
            continue

        if current[0] == 'initial_alpha': 
            current[1] = str(2.8197*np.pi/180) + '\n'
            nu_lines.append(' = '.join(current))
            continue

        if current[0] == 'initial_deflection':
            current[1] = str(-1*np.pi/180) + '\n'
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

    # 1G Static Trim case
    data2 = sharpy_main.main(['', route + '/' + case_name + '.sharpy'])
    
    grav_forces = data2.structure.timestep_info[0].gravity_forces
    grav = 0
    for i in range(len(grav_forces)):
        # grav += grav_forces[i][2]
        grav += np.sqrt(grav_forces[i][0]**2 + grav_forces[i][1]**2 + grav_forces[i][2]**2)
    
    masss_temp = abs(grav/(9.81*g_loading)) 
    masss.append(masss_temp) 
       
    wing_tip_deflection_log.append(data2.structure.timestep_info[-1].pos[beam_opt_nodes][-1][2])
    
    max_deflection_val = data2.structure.elements[beam_opt_elems[-1]].coordinates_ini[1][2] + max_deflection
    max_deflection_val_log.append(max_deflection_val)
    deflection = data2.structure.timestep_info[-1].pos[beam_opt_nodes[-1]][2]  
    
    wing_deflection_diff = (data2.structure.elements[beam_opt_elems[-1]].coordinates_ini[1][2] + max_deflection) - data2.structure.timestep_info[-1].pos[beam_opt_nodes[-1]][2]
    wing_deflection_diff_log.append(wing_deflection_diff)
    
    
    if data2.structure.div_flag.value != 1:
        node_forces2 = local_for_forces(data2) 
        power = data2.thrust * u_inf
        alpha_log.append(data2.alpha)
        elevator_log.append(data2.elevator)

    # if there was divergence of 1g Trim
    else:  
        power = starting_power*5 # handled by increase in mass 
        alpha_log.append(0)
        elevator_log.append(0)
    
            
    power_log.append(power)
    cost = (power / starting_power) 
    chord_log.append(chord_vals[0])        
    diff_sum_out.append(cost)
    ei_log.append(data2.structure.stiffness_db[0][4][4])
    gj_log.append(data2.structure.stiffness_db[0][3][3])
    mass_log.append(data2.structure.mass_db[0][0][0])
    
     
    df = pd.DataFrame([guess_log, power_log, masss, chord_log, alpha_log, elevator_log, wing_tip_deflection_log, max_deflection_val_log, ei_log, gj_log, mass_log])
    df.to_excel( 'temp_' + beam_save_file, header= False,  sheet_name= 'CL = {0}'.format(CL))
    if CL == CL_range[0]:
        df.to_excel(beam_save_file, header= False,  sheet_name= 'CL = {0}'.format(CL))
    else:
        if guess[0] == spans[-1]:
            with pd.ExcelWriter(beam_save_file, engine='openpyxl', mode='a') as writer:  
                df.to_excel(writer, header = False, sheet_name= 'CL = {0}'.format(CL))
        
    return cost 

# results = minimize(span_opt, guess, bounds = bnds)
# print(results)
# beam_opt(results)

global CL, CL_range, spans
CL_range = [.2, .3, .4, .5, .6, .7, .8, .9]
spans = np.linspace(10, 80, 30)
# spans = [32]
for i in CL_range:
    CL = i
    guess_log = []
    diff_sum_out = []
    masss = []
    max_strain_log = []
    wing_tip_deflection_log = []
    max_deflection_val_log = []
    wing_deflection_diff_log = []
    power_log = []
    chord_log = []
    penalties_log = []
    alpha_log = []
    elevator_log = []
    ei_log = []
    gj_log = []
    mass_log = []
    for j in spans:
        span_opt([j])

script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, 'results/')

if not os.path.isdir(results_dir):
    os.makedirs(results_dir)

# fig = plt.figure(figsize= (8, 8))
# ax1 = fig.add_subplot(211)
# ax1.plot(span_loc_log, thickness_log, 'r') 
# ax1.set_title('Thickness vs Span', size = 16)
# ax1.set_xlabel('Span (m)', size = 12)
# ax1.set_ylabel('Thickness (m)', size = 12)

# ax2 = fig.add_subplot(212)
# ax2.plot(span_loc_log, area_spar_log, 'b')
# ax2.set_title('Area vs Span', size = 16)
# ax2.set_xlabel('Span (m)', size = 12)
# ax2.set_ylabel('Area (m^2)', size = 12)

# fig.tight_layout()
# fig.savefig(results_dir +'thickness_area_' + beam_save_file+ str(g_loading) + 'G.png')

# #----------------------------------------------------------------

# fig = plt.figure(figsize= (8, 8))
# ax1 = fig.add_subplot(211)
# ax1.plot(span_loc_log, ei_span_log, 'r') 
# ax1.set_title('EI_yy vs Span', size = 16)
# ax1.set_xlabel('Span (m)', size = 12)
# ax1.set_ylabel('EI_yy (kg-m^2)', size = 12)

# ax2 = fig.add_subplot(212)
# ax2.plot(span_loc_log, gj_span_log, 'b')
# ax2.set_title('GJ vs Span', size = 16)
# ax2.set_xlabel('Span (m)', size = 12)
# ax2.set_ylabel('GJ (kg-m^2)', size = 12)

# fig.tight_layout()
# fig.savefig(results_dir +'EI_GJ_' + beam_save_file + str(g_loading) + 'G.png')

# #---------------------------------------------------------------

# fig = plt.figure(figsize= (7, 8))

# ax1 = fig.add_subplot(311)
# ax1.plot(span_loc_log, shear_log, 'r') 
# ax1.set_title('Shear vs Span', size = 16)
# ax1.set_xlabel('Span (m)', size = 12)
# ax1.set_ylabel('Shear (N)', size = 12)

# ax2 = fig.add_subplot(312)
# ax2.plot(span_loc_log, moment_log, 'g')
# ax2.set_title('Moment vs Span', size = 16)
# ax2.set_xlabel('Span (m)', size = 12)
# ax2.set_ylabel('Moment (N-m)', size = 12)

# ax3 = fig.add_subplot(313)
# ax3.plot(span_loc_log, strain_log, 'b')
# ax3.set_title('Strain vs Span', size = 16)
# ax3.set_xlabel('Span (m)', size = 12)
# ax3.set_ylabel('Strain', size = 12)

# fig.tight_layout()
# fig.savefig(results_dir +'Shaer_Moment_Strain_' + beam_save_file + str(g_loading) + 'G.png')

# #---------------------------------------------------------------
# plt.figure(figsize = (12,6))
# plt.rc('axes', titlesize = 16)
# plt.rc('axes', labelsize = 12)
# plt.rc('legend', fontsize = 10)
# plt.plot(span_loc_log, node_forces_z_log, 'r', label = 'Net-Force')
# plt.plot(span_loc_log, local_grav_forces_log[0:len(beam_opt_nodes), 2], 'g', label = 'Grav-Force')
# plt.plot(span_loc_log, local_applied_forces_log[0:len(beam_opt_nodes), 2], 'b', label = 'Lift-Force')
# plt.title('Aircraft Loads')
# plt.xlabel('Span (m)')
# plt.ylabel('Load (N)')
# plt.legend(loc = 'lower right')
# plt.savefig(results_dir +'Loads_' + beam_save_file + str(g_loading) + 'G.png')

# #---------------------------------------------------------------
# plt.figure(figsize = (12,6))
# plt.rc('axes', titlesize = 16)
# plt.rc('axes', labelsize = 12)
# plt.rc('legend', fontsize = 10)
# plt.plot(span_loc_log, lift[0:len(beam_opt_nodes), 2], 'r', label = 'Lift-Force')
# plt.title('Lift Distribution')
# plt.xlabel('Span (m)')
# plt.ylabel('Lift (N)')
# plt.savefig(results_dir +'Lift_Distribution_' + beam_save_file + '.png')

#df = pd.DataFrame([span_loc_log, thickness_log, area_spar_log, node_forces_z_log, local_grav_forces_log[0:len(beam_opt_nodes), 2], local_applied_forces_log[0:len(beam_opt_nodes), 2], ei_span_log, gj_span_log, shear_log, moment_log, strain_log])

# with pd.ExcelWriter(beam_save_file, engine='openpyxl', mode='a') as writer:  
#     df.to_excel(writer, header= False, sheet_name= 'Results From Opt')


# # tx1 = 'Span (m)'
# # tx2 = 'Thickness (m)'
# # tx3 = 'Area of Spar (m^2)'j
# # tx4 = 'Local Node Forces (N)'
# # tx5 = 'Local Grav_Forces_Z (N)'
# # tx6 = 'Local App_Forces_Z (N)'
# # tx7 = 'EI_y (N-m^2)'
# # tx8 = 'GJ (N-m^2)'
# # tx9 = 'Shear (N)'
# # tx10 = 'Moment (N-m)'
# # tx11 = 'Strain'

# # writer = pd.ExcelWriter(beam_save_file, engine="openoyxl")
# # act = writer.sheets['Results From Opt']
# # act.write(0,0, tx1)
# # act.write(1,0, tx2)
# # act.write(2,0, tx3)
# # act.write(3,0, tx4)
# # act.write(4,0, tx5)
# # act.write(5,0, tx6)
# # act.write(6,0, tx7)
# # act.write(7,0, tx8)
# # act.write(8,0, tx9)
# # act.write(9,0, tx10)
# # act.write(10,0, tx11)



# # # data = sharpy_main.main(['', route + '/' + case_name + '.sharpy'])


