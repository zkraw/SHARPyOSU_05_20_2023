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
import scipy
 
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
file = 'T_Tail.xlsx'
global case_name
case_name = 'test_Beam'
excel = True 

input_file = os.path.dirname(os.path.abspath(__file__)) + '/' + file

flow = ['BeamLoader',
        'AerogridLoader',
        #'NonLinearStatic',
        #'StaticUvlm',
        'BeamOpt',
        #'StaticTrim',
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
u_inf = 10 #14 for XHALE, 10 for T_Tail
rho = 1.225 # 1.225 
3
# default elevator is -0.020708598 radians 3
alpha = 2.2793*np.pi/180 # 1G 2.0251 alpha, .0533 N thrust, 3.7940 degrees elevator / 2G 4.287 degrees alpha, .2173 N thrust, 6.9285 degrees elevator 
beta = 0
roll = 0
gravity = 'on'
g_loading = 2 # acceleration of gravity (9.81) * factor 

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
cs_deflection = 2.9402*np.pi/180 # initial control surface deflection 1G .06621779 / 2G 0.1209251 
tail_cs_index = run.tail_index # the control surface index that will be used 
thrust = 1.332 # 1G .0533 N  / 2G .2173 N  

global beam_opt_nodes
global beam_opt_sym
#Beam Opt Settings
max_strain = .0001
beam_save_file = 'T_Tail_1G.xlsx'
save_data = True
beam_opt_nodes = run.beam_opt_nodes
beam_opt_elems = run.beam_opt_elems
beam_opt_sym = run.beam_opt_sys
max_iter = 5


# -----------------------------------------------------------------------------------------------------
#----------------------End of Input--------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

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
                                'relaxation_factor': structural_relaxation_factor}




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
                        'Max_iter' : max_iter,
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
    fem_file_name = 'Zack_Beam' + '/' + case_name + '.aero.h5'
    if os.path.isfile(fem_file_name):
        os.remove(fem_file_name)
        
def write_aero(aero):
    
    with h5.File('Zack_Beam' + '/' + case_name + '.aero.h5', 'a') as h5file:
        airfoils_group = h5file.create_group('airfoils')

        for i in range(len(aero['airfoils'])):  # len(air_data)
            airfoil = airfoils_group.create_dataset(str(i), data= aero['airfoils'][i])
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

        
def twist_opt(twist):
    span_locs = np.linspace(0, 1, len(twist))
    interp = scipy.interpolate.interp1d(x=span_locs,
                                     y=twist,
                                     kind='quadratic')
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
    
    max_span =  fem['coordinates'][beam_opt_nodes[-1]][1] 

    # changing twist values for optimization
    for i in range(len(beam_opt_nodes)):
        for j in range(len(fem['connectivities'])):
            for k in [0, 2, 1]:
                if beam_opt_nodes[i] == fem['connectivities'][j][k] and aero['chord'][j][k] != 0:
                    aero['twist'][j][k] = interp(fem['coordinates'][fem['connectivities'][j][k]][1] / max_span)
    
    
    #symmetric data 
    if beam_opt_sym == True:
        for i in range(1, len(beam_opt_nodes)):
            for j in range(len(fem['connectivities'])):
                for k in [0, 2, 1]:
                    if beam_opt_nodes[-1] + beam_opt_nodes[i] == fem['connectivities'][j][k] and aero['chord'][j][k] != 0:
                        aero['twist'][j][k] = interp(abs(fem['coordinates'][fem['connectivities'][j][k]][1] / max_span))
    
    # deletes previous aero file 
    clean_test_files()

    # creates new aero file with updated twist
    write_aero(aero)
    
    # runs beam_opt
    data = sharpy_main.main(['', route + '/' + case_name + '.sharpy'])
    
    
    return # difference from max strain + drag/power 
    
# twist_opt([2., 1, 0, 0])


data = sharpy_main.main(['', route + '/' + case_name + '.sharpy'])


