import h5py as h5
import numpy as np
import configparser
import os
import configobj
import json
import sys
from pathlib import Path
import pandas as pd 
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
case_name = 'test_Beam'
excel = True 

input_file = os.path.dirname(os.path.abspath(__file__)) + '/' + file

flow = ['BeamLoader',
        'AerogridLoader',
        'NonLinearStatic',
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
rho = 1.225 # DO NOT CHANGE THIS

# default elevator is -0.020708598 radians 
alpha = 1.4641*np.pi/180 # 1G 2.0251 alpha, .0533 N thrust, 3.7940 degrees elevator / 2G 4.287 degrees alpha, .2173 N thrust, 6.9285 degrees elevator 
beta = 0
roll = 0
gravity = 'on'
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
dt = .1
n_tstep = round(physical_time/dt)

# saving values to new file
#Paramaterized.excel_write(file, 'new_file.xlsx', ['start_Y'], [[3, 4, 5, 6, 7, 8]], [[123456789, 100, 100, 100, 100, 100]])

# generate FEM and Aero from input file
run = Paramaterized(case_name, input_file, excel)


# Static Trim Settings
nodes_app = run.app_nodes # where the thrust nodes are 
cs_deflection = 0 # initial control surface deflection 1G .06621779 / 2G 0.1209251 
tail_cs_index = run.tail_index # the control surface index that will be used 
thrust = 5 # 1G .0533 N  / 2G .2173 N  

#Beam Opt Settings
max_strain_space = np.linspace(.00001, .001 , 50)
beam_save_file = 'T_Tail_1G.xlsx'
mass_strain = []

beam_opt_nodes = run.beam_opt_nodes
beam_opt_elems = run.beam_opt_elems
beam_opt_sym = run.beam_opt_sys
max_iter = 10

# -----------------------------------------------------------------------------------------------------
#----------------------End of Input--------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

# making sure there is not a save file with the name already 
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
            save_elem[0][-2] = str(count)     
            count += 1
    else:
        check = False

for i in range(len(max_strain_space)):
    #Beam Opt Settings
    max_strain = max_strain_space[i]

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
        with_dynamic_forces = False
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


        settings['BeamOpt'] = {'solver': 'StaticCoupled',
                            'solver_settings': settings['StaticCoupled'],
                            'Max Strain': max_strain,
                            'Spar couple H' :.86,
                            'Beam_Opt_Nodes': beam_opt_nodes,
                            'Beam_Opt_Elems' : beam_opt_elems,
                            'Sym': beam_opt_sym,
                            'Max_iter' : max_iter,
                            'initial_deflection': cs_deflection}


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
    # generate_dyn_file()

    try:
        data = sharpy_main.main(['', route + '/' + case_name + '.sharpy'])

        grav_forces = data.structure.timestep_info[0].gravity_forces
        grav_z = 0
        for i in range(len(grav_forces)):
            grav_z += grav_forces[i][2]

        # CG calculation
        cg = []
        mass_bar = []
        elem_start = []
        elem_end = []

        # outputs
        r = []
        m = []
        span = []
        mr = []

        for i in range(len(data.structure.mass_db)):
            mass_db = data.structure.mass_db[i]
            mass_bar.append(mass_db[0][0])
            elem_start.append(data.structure.elements[i].coordinates_ini[0])
            elem_end.append(data.structure.elements[i].coordinates_ini[1])
            cg.append([mass_db[1][5]/mass_db[0][0], mass_db[2][3]/mass_db[0][0], mass_db[0][4]/mass_db[0][0]])
            
            # right wing 
            if run.ref_lookup[i] == 'right': 
                r.append((elem_start[i] + elem_end[i])/2 + [-1*cg[i][1], cg[i][0], cg[i][2]])
                m.append(abs(elem_start[i][1] - elem_end[i][1])*mass_bar[i])
                mr.append(m[i]*r[i])
            
            # left wing
            if run.ref_lookup[i] == 'left': 
                r.append((elem_start[i] + elem_end[i])/2 + [1*cg[i][1], cg[i][0], cg[i][2]])
                m.append(abs(elem_start[i][1] - elem_end[i][1])*mass_bar[i])
                mr.append(m[i]*r[i])
            
            # v_fin 
            if run.ref_lookup[i] == 'v_fin': 
                r.append((elem_start[i] + elem_end[i])/2 + [-1*cg[i][1], cg[i][0], cg[i][2]])
                m.append(abs(elem_start[i][2] - elem_end[i][2])*mass_bar[i])
                mr.append(m[i]*r[i])
            
            # fuselage assumed that mass is located on the Elastic Axis
            if run.ref_lookup[i] == 'fuselage': 
                r.append((elem_start[i] + elem_end[i])/2)
                m.append(abs(elem_start[i][0] - elem_end[i][0])*mass_bar[i])
                mr.append(m[i]*r[i])

            # fuselage_neg assumed that mass is located on the Elastic Axis
            if run.ref_lookup[i] == 'fuselage_neg': 
                r.append((elem_start[i] + elem_end[i])/2) 
                m.append(abs(elem_start[i][0] - elem_end[i][0])*mass_bar[i])
                mr.append(m[i]*r[i])

        for i in range(len(run.lumped_loc)):
            mr.append(data.structure.lumped_mass[i]*run.lumped_loc[i])

        # total mass with the addition of the lumped masses
        try:
            m_total = sum(m) + sum(data.structure.lumped_mass)    
        except:
            m_total = sum(m)
            
        mr_total = sum(mr) 
        cg_total = mr_total / m_total
        print('\nMass From Applied Forces: ' + str(abs(grav_z/(9.81*g_loading))+ sum(data.structure.lumped_mass)) + '\n')
        print('Initial Span Mass: ' + str(m_total[0]) + '\n')
        print('x_cg = ' + str(cg_total[0]) + ' y_cg = ' + str(cg_total[1]) + ' z_cg = ' + str(cg_total[2]) + '\n')
        
        massss = abs(grav_z/(9.81*g_loading)) + sum(data.structure.lumped_mass) 
        mass_strain.append(massss[0])
        df = pd.DataFrame([max_strain_space, mass_strain])
        df.to_excel(beam_save_file, index=['Strain', 'Mass'], header= False) 
    except:
        mass_strain.append('0')
        df = pd.DataFrame([max_strain_space, mass_strain])
        df.to_excel(beam_save_file, index=['Strain', 'Mass'], header= False) 
      