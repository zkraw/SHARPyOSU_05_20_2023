import configobj
import numpy as np
import sharpy.utils.algebra as algebra
import os


def t_tail_static(route, case_name, set):
    # initial sharpy settings for t_tail HALE span_opt paper
    flow = set['flow']
    roll = set['roll']
    alpha = set['alpha']
    beta = set['beta']
    
    tstep_factor = set['tstep_factor']
    tolerance = set['tolerance']
    gravity = set['gravity']
    g_loading = set['g_loading']
    num_cores = set['num_cores']
    dt = set['dt']
    u_inf = set['u_inf']
    rho = set['rho']
    n_step = set['n_step']
    
    fsi_tolerance = set['fsi_tolerance']
    structural_relaxation_factor = set['structural_relaxation_factor']
    n_step = set['n_step']
    relaxation_factor = set['relaxation_factor']

    cs_deflection = set['cs_deflection']
    tail_cs_index = set['tail_cs_index']
    nodes_app = set['nodes_app']
    thrust = set['initial_thrust']
    
    n_tstep = set['n_tstep']

    free_flight = set['free_flight']

    file_name = os.path.join(route, '{}.sharpy'.format(case_name))
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
                                'max_iter': 200,
                                'n_load_steps': n_step,
                                'tolerance': fsi_tolerance,
                                'relaxation_factor': structural_relaxation_factor}


    settings['StaticTrim'] = {'solver': 'StaticCoupled',
                            'solver_settings': settings['StaticCoupled'],
                            'initial_alpha': alpha,
                            'initial_deflection': cs_deflection,
                            'tail_cs_index': tail_cs_index,
                            'thrust_nodes': nodes_app,
                            'initial_thrust': thrust,
                            'max_iter': 200,
                            'fz_tolerance': 0.1,
                            'fx_tolerance': 0.1,
                            'm_tolerance': 0.1}

    # settings['BeamOpt'] = {'solver': 'StaticCoupled',
    #                     'solver_settings': settings['StaticCoupled'],
    #                     'Max Strain': 0,
    #                     'Beam_Opt_Nodes': beam_opt_nodes,
    #                     'Beam_Opt_Elems' : beam_opt_elems,
    #                     'Sym': beam_opt_sym,
    #                     'initial_deflection': cs_deflection}

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

def steady_dyn(route, case_name, set):
    # steady velocity field generation to see what happens after running a configuration through StaticTrim
    flow = set['flow']
    # if static trim is ran before dynamiccoupled, velocity field is already defined. 
    # Adding another steady velocity field will essentially double the velocity
    if 'StaticTrim' in flow:
        trim_check = 0
    else:
        trim_check = 1
    
    roll = set['roll']
    alpha = set['alpha']
    beta = set['beta']
    
    tstep_factor = set['tstep_factor']
    tolerance = set['tolerance']
    gravity = set['gravity']
    g_loading = set['g_loading']
    num_cores = set['num_cores']
    dt = set['dt']
    u_inf = set['u_inf']
    rho = set['rho']
    n_step = set['n_step']
    
    fsi_tolerance = set['fsi_tolerance']
    structural_relaxation_factor = set['structural_relaxation_factor']
    n_step = set['n_step']
    relaxation_factor = set['relaxation_factor']

    cs_deflection = set['cs_deflection']
    tail_cs_index = set['tail_cs_index']
    nodes_app = set['nodes_app']
    initial_thrust = set['initial_thrust']
    
    n_tstep = set['n_tstep']

    free_flight = set['free_flight']
                                                     
    gust_length = set['gust_length']
    gust_intensity = set['gust_intensity']
    gust_offset = set['gust_offset']

    file_name = os.path.join(route, '{}.sharpy'.format(case_name))
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
                            'initial_alpha': alpha,
                            'initial_deflection': cs_deflection,
                            'tail_cs_index': tail_cs_index,
                            'thrust_nodes': nodes_app,
                            'initial_thrust': initial_thrust,
                            'fz_tolerance': .00001,
                            'fx_tolerance': .00001,
                            'm_tolerance': .00001,
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
                            'velocity_field_generator': 'SteadyVelocityField',
                            'velocity_field_input': {'u_inf': trim_check*u_inf,
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

def linearassembler(route, case_name, set):

    flow = ['BeamLoader',
        'AerogridLoader',
        'StaticTrim',
        'BeamPlot',
        'AerogridPlot',
        'AeroForcesCalculator',
        'DynamicCoupled',
        'Modal',
        'LinearAssembler',
        'AsymptoticStability',
        ]

    roll = set['roll']
    alpha = set['alpha']
    beta = set['beta']
    
    tstep_factor = set['tstep_factor']
    tolerance = set['tolerance']
    gravity = set['gravity']
    g_loading = set['g_loading']
    num_cores = set['num_cores']
    dt = set['dt']
    u_inf = set['u_inf']
    rho = set['rho']
    
    fsi_tolerance = set['fsi_tolerance']
    structural_relaxation_factor = set['structural_relaxation_factor']
    n_step = set['n_step']
    relaxation_factor = set['relaxation_factor']

    cs_deflection = set['cs_deflection']
    tail_cs_index = set['tail_cs_index']
    nodes_app = set['nodes_app']
    initial_thrust = set['initial_thrust']
    
    n_tstep = set['n_tstep']

    free_flight = set['free_flight']
                                                     
    gust_length = set['gust_length']
    gust_intensity = set['gust_intensity']
    gust_offset = set['gust_offset']

    # numerics default for T-Tail
    n_step = 5
    structural_relaxation_factor = 0.6
    relaxation_factor = 0.35
    tolerance = 1e-6
    fsi_tolerance = 1e-4
    num_cores = 2

    m_star = set['m_star']
    file_name = route + '/' + case_name + '.sharpy'
    settings = dict()
    settings['SHARPy'] = {'case': case_name,
                        'route': route,
                        'flow': flow,
                        'write_screen': 'on',
                        'write_log': 'on',
                        'log_folder': route + '/output',
                        'log_file': case_name + '.log'}

    settings['BeamLoader'] = {'unsteady': 'off',
                            'orientation': algebra.euler2quat(np.array([roll,
                                                                        alpha,
                                                                        beta]))}
    settings['AerogridLoader'] = {'unsteady': 'off',
                                'aligned_grid': 'on',
                                'mstar': int(4*5),
                                'freestream_dir': ['1', '0', '0'],
                                'control_surface_deflection': [''],
                                'wake_shape_generator': 'StraightWake',
                                'wake_shape_generator_input': {'u_inf': u_inf,
                                    'u_inf_direction': ['1', '0', '0'],
                                    'dt': dt}}

    settings['NonLinearStatic'] = {'print_info': 'off',
                                'max_iterations': 100,
                                'num_load_steps': 1,
                                'delta_curved': 1e-1,
                                'min_delta': tolerance,
                                'gravity_on': gravity,
                                'gravity': 9.81*g_loading}

    settings['StaticUvlm'] = {'print_info': 'off',
                            'horseshoe': 'off',
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
                             'structural_solver_settings': {'print_info': 'off',
                                                            'max_iterations': 200,
                                                            'num_load_steps': 1,
                                                            'delta_curved': 1e-5,
                                                            'min_delta': tolerance,
                                                            'gravity_on': 'on',
                                                            'gravity': 9.81},
                             'aero_solver': 'StaticUvlm',
                             'aero_solver_settings': {'print_info': 'on',
                                                      'horseshoe': False,
                                                      'num_cores': 4,
                                                      'n_rollup': int(0),
                                                      'rollup_dt': dt,
                                                      'rollup_aic_refresh': 1,
                                                      'rollup_tolerance': 1e-4,
                                                      'velocity_field_generator': 'SteadyVelocityField',
                                                      'velocity_field_input': {'u_inf': u_inf,
                                                                               'u_inf_direction': [1., 0, 0]},
                                                      'rho': rho},
                             'max_iter': 200,
                             'n_load_steps': 1,
                             'tolerance': tolerance,
                             'relaxation_factor': 0.2}

    settings['StaticTrim'] = {'solver': 'StaticCoupled',
                            'solver_settings': settings['StaticCoupled'],
                            'initial_alpha': alpha,
                            'initial_deflection': cs_deflection,
                            'tail_cs_index': tail_cs_index,
                            'thrust_nodes': nodes_app,
                            'initial_thrust': initial_thrust,
                            'fz_tolerance': 1e-2,
                          'fx_tolerance': 1e-2,
                          'm_tolerance': 1e-2}


    settings['NonLinearDynamicCoupledStep'] = {'print_info': 'off',
                                            'max_iterations': 500,
                                            'delta_curved': 1e-1,
                                            'min_delta': tolerance,
                                            'newmark_damp': 5e-3,
                                            'gravity_on': gravity,
                                            'gravity': 9.81*g_loading,
                                            'num_steps': n_tstep,
                                            'dt': dt,
                                            'initial_velocity': u_inf}

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
                            # 'horseshoe': 'off',
                            'num_cores': 4,
                            'convection_scheme': 2,
                            'gamma_dot_filtering': 3,
                            'velocity_field_generator': 'SteadyVelocityField',
                            'velocity_field_input': {'u_inf': u_inf*0,
                                               'u_inf_direction': [1., 0., 0.]},
                            # 'velocity_field_input': {'u_inf': int(not free_flight) * u_inf,
                            #                          'u_inf_direction': [1., 0, 0],
                            #                          'gust_shape': '1-cos',
                            #                          'gust_parameters': {'gust_length': gust_length,
                            #                                              'gust_intensity': gust_intensity * u_inf},
                            #                          'offset': gust_offset,
                            #                          'relative_motion': relative_motion},
                            'rho': rho,
                            'n_time_steps': n_tstep,
                            'dt': dt}

    if free_flight:
        solver = 'NonLinearDynamicCoupledStep'
    else:
        solver = 'NonLinearDynamicPrescribedStep'
        
    struct_solver_settings = {'print_info': 'off',
                          'initial_velocity_direction': [-1., 0., 0.],
                          'max_iterations': 950,
                          'delta_curved': 1e-6,
                          'min_delta': tolerance,
                          'newmark_damp': 5e-3,
                          'gravity_on': True,
                          'gravity': 9.81,
                          'num_steps': n_tstep,
                          'dt': dt,
                          'initial_velocity': u_inf * 1}

    settings['DynamicCoupled'] = {'structural_solver': solver,
                                'structural_solver_settings': struct_solver_settings,
                                'aero_solver': 'StepUvlm',
                                'aero_solver_settings': settings['StepUvlm'],
                                'fsi_substeps': 200,
                                'fsi_tolerance': fsi_tolerance,
                                'relaxation_factor': relaxation_factor,
                                'minimum_steps': 1,
                                # 'relaxation_steps': 150,
                                # 'final_relaxation_factor': 0.5,
                                'n_time_steps': n_tstep,
                                'dt': dt,
                                'include_unsteady_force_contribution': 'off'
                                # 'postprocessors': ['BeamLoads', 'BeamPlot', 'AerogridPlot'],
                                # 'postprocessors_settings': {'BeamLoads': {'csv_output': 'off'},
                                #                             'BeamPlot': {'include_rbm': 'on',
                                #                                         'include_applied_forces': 'on'},
                                #                             'AerogridPlot': {
                                #                                 'include_rbm': 'on',
                                #                                 'include_applied_forces': 'on',
                                #                                 'minus_m_star': 0},
                                                            }

    settings['BeamLoads'] = {'csv_output': 'off'}

    settings['BeamPlot'] = {'include_rbm': 'on',
                            'include_applied_forces': 'on',
                            'include_FoR': 'on'}
                        

    settings['AerogridPlot'] = {'include_rbm': 'off',
                                'include_applied_forces': 'on',
                                'minus_m_star': 0,
                                'u_inf': u_inf
                                }
    
    settings['AeroForcesCalculator'] = {
                                    'write_text_file': 'off',
                                    'text_file_name': case_name + '_aeroforces.csv',
                                    'screen_output': 'on',
                                    'coefficients': True,
                                    'q_ref': 0.5 * rho * u_inf ** 2,
                                    'S_ref': 162,
                                    }

    settings['Modal'] = {'print_info': True,
                    'use_undamped_modes': True,
                    'NumLambda': 30,
                    'rigid_body_modes': True,
                    'write_modes_vtk': 'on',
                    'print_matrices': 'off',
                    #  'write_data': 'on',
                    'continuous_eigenvalues': 'off',
                    'dt': dt,
                    'plot_eigenvalues': False,
                    'rigid_modes_cg': False}

    settings['LinearAssembler'] = {'linear_system': 'LinearAeroelastic',
                                    'linear_system_settings': {
                                        'beam_settings': {'modal_projection': True,
                                                        'inout_coords': 'modes',
                                                        'discrete_time': True,
                                                        'newmark_damp': 0.05,
                                                        'discr_method': 'newmark',
                                                        'dt': dt,
                                                        'proj_modes': 'undamped',
                                                        # 'use_euler': 'off',
                                                        'num_modes': 9,
                                                        'print_info': 'on',
                                                        'gravity': 'on',
                                                        'remove_dofs': []},
                                        'aero_settings': {'dt': dt,
                                                        'integr_order': 2,
                                                        'density': rho,
                                                        'remove_predictor': False,
                                                        'use_sparse': False,
                                                        # 'rigid_body_motion': free_flight,
                                                        # 'use_euler': False,
                                                        'remove_inputs': ['u_gust']},
                                        'track_body': True,
                                        'use_euler': True,
                                        
                                        }}
                                        # 'rigid_body_motion': free_flight}}

    settings['AsymptoticStability'] = { 'print_info': 'on',
                                        # 'modes_to_plot': [],
                                        # 'display_root_locus': 'off',
                                        'frequency_cutoff': 0,
                                        'export_eigenvalues': 'on',
                                        # 'sys_id': 'LinearAeroelastic',
                                        'num_evals': 500,}
                                    

    config = configobj.ConfigObj()
    config.filename = file_name
    for k, v in settings.items():
        config[k] = v
    config.write() 

def turbulent_velocity_feild(route, case_name, set):
    # initial sharpy settings for t_tail HALE span_opt paper
    flow = set['flow']
    roll = set['roll']
    alpha = set['alpha']
    beta = set['beta']
    
    tstep_factor = set['tstep_factor']
    tolerance = set['tolerance']
    gravity = set['gravity']
    g_loading = set['g_loading']
    num_cores = set['num_cores']
    dt = set['dt']
    u_inf = set['u_inf']
    rho = set['rho']
    
    fsi_tolerance = set['fsi_tolerance']
    structural_relaxation_factor = set['structural_relaxation_factor']
    n_step = set['n_step']
    relaxation_factor = set['relaxation_factor']

    cs_deflection = set['cs_deflection']
    tail_cs_index = set['tail_cs_index']
    nodes_app = set['nodes_app']
    initial_thrust = set['initial_thrust']
    
    n_tstep = set['n_tstep']

    free_flight = set['free_flight']
                                                     
    gust_length = set['gust_length']
    gust_intensity = set['gust_intensity']
    gust_offset = set['gust_offset']

    # numerics default for T-Tail
    n_step = 5
    structural_relaxation_factor = 0.6
    relaxation_factor = 0.35
    tolerance = 1e-6
    fsi_tolerance = 1e-4
    num_cores = 2

    m_star = set['m_star']

    file_name = os.path.join(route, '{}.sharpy'.format(case_name))
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
                                'mstar': int(m_star/tstep_factor),
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
                            'initial_alpha': alpha,
                            'initial_deflection': cs_deflection,
                            'tail_cs_index': tail_cs_index,
                            'thrust_nodes': nodes_app,
                            'initial_thrust': initial_thrust,
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
                            'velocity_field_generator': 'TurbVelocityField',
                            'velocity_field_input': {'turbulent_field': set['xmf'], # file location of xdmf
                                                     'offset': set['offset'], # spatial offset in the 3 dimensions
                                                     'store_field': False, # stores a snapshot of the xdmf in memory. saves from allocating GB's of memory to the velocity field
                                                     'centre_y': False, # Flat for changing the domain to ['-y_max/2', 'y_max/2'] CHECK to see how it affects results 
                                                     'periodicity': 'xy', # axes in which periodicity is enforced
                                                     'frozen': False}, # If 'True', the turbulent field will not be updtaed in time 
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

def turbulent_velocity_feild_controller(route, case_name, set):
    # initial sharpy settings for t_tail HALE span_opt paper
    flow = set['flow']
    roll = set['roll']
    alpha = set['alpha']
    beta = set['beta']
    
    tstep_factor = set['tstep_factor']
    tolerance = set['tolerance']
    gravity = set['gravity']
    g_loading = set['g_loading']
    num_cores = set['num_cores']
    dt = set['dt']
    u_inf = set['u_inf']
    rho = set['rho']
    
    fsi_tolerance = set['fsi_tolerance']
    structural_relaxation_factor = set['structural_relaxation_factor']
    n_step = set['n_step']
    relaxation_factor = set['relaxation_factor']

    cs_deflection = set['cs_deflection']
    tail_cs_index = set['tail_cs_index']
    nodes_app = set['nodes_app']
    initial_thrust = set['initial_thrust']
    
    n_tstep = set['n_tstep']

    free_flight = set['free_flight']
                                                     
    gust_length = set['gust_length']
    gust_intensity = set['gust_intensity']
    gust_offset = set['gust_offset']

    # numerics default for T-Tail
    n_step = 5
    structural_relaxation_factor = 0.6
    relaxation_factor = 0.35
    tolerance = 1e-6
    fsi_tolerance = 1e-4
    num_cores = 2

    m_star = set['m_star']

    file_name = os.path.join(route, '{}.sharpy'.format(case_name))
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
                                'mstar': int(m_star/tstep_factor),
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
                            'initial_alpha': alpha,
                            'initial_deflection': cs_deflection,
                            'tail_cs_index': tail_cs_index,
                            'thrust_nodes': nodes_app,
                            'initial_thrust': initial_thrust,
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
                            'velocity_field_generator': 'TurbVelocityField',
                            'velocity_field_input': {'turbulent_field': set['xmf'], # file location of xdmf
                                                     'offset': set['offset'], # spatial offset in the 3 dimensions
                                                     'store_field': False, # stores a snapshot of the xdmf in memory. saves from allocating GB's of memory to the velocity field
                                                     'centre_y': False, # Flat for changing the domain to ['-y_max/2', 'y_max/2'] CHECK to see how it affects results 
                                                     'periodicity': 'xy', # axes in which periodicity is enforced
                                                     'frozen': True}, # If 'True', the turbulent field will not be updtaed in time 
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
                                'fsi_substeps': 50,
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
                                                            },

                                'controller_id': {'roll': 'ControlSurfacePidController',
                                                #   'vel_hold': 'ControlSurfacePidController',
                                                  'altitude': 'ControlSurfacePidController'
                                                },
                                      
                                'controller_settings': {
                                #     'vel_hold':{'P': .7,
                                #                 'I': 0,
                                #                 'D': 0,
                                #                 'input_type': 'for_vel_x',
                                #                 'dt': dt,
                                #                 'controlled_surfaces': [0] , # note this does is not needed for this controller, but there is an error otherwise
                                #                 'controller_log_route': os.path.join(route, 'output', case_name),
                                #                 'time_history_input_file': os.path.join(route, 'x_vel.csv')
                                                                    # },
                                    'roll':{'P': 3.5,
                                                'I': 6,
                                                'D': 0,
                                                'input_type': 'roll',
                                                'dt': dt,
                                                'controlled_surfaces': [0, 1, 2, 3] , # list of integers that specify the control surface ID's to be used 
                                                'controlled_surfaces_coeff': [1, -1, 0, 0], # allows for antisymmetric deflection of aileron
                                                'controller_log_route': os.path.join(route, 'output', case_name),
                                                'time_history_input_file': os.path.join(route, 'roll.csv')
                                                                     },

                                    'altitude':{'P': .00875,
                                                'I': -.0005,
                                                'D': 0,
                                                'input_type': 'for_pos_z',
                                                'dt': dt,
                                                'controlled_surfaces': [0, 1, 2, 3] , # list of integers that specify the control surface ID's to be used 
                                                'controlled_surfaces_coeff': [0, 0, 0, -1], # allows for antisymmetric deflection of aileron
                                                'controller_log_route': os.path.join(route, 'output', case_name),
                                                'time_history_input_file': os.path.join(route, 'altitude.csv')
                                                                     }}

                                    # 'altitude':{'P': 6.2,
                                    #             'I': 3.2,
                                    #             'D': 0,
                                    #             'input_type': 'pitch',
                                    #             'dt': dt,
                                    #             'controlled_surfaces': [0, 1, 2, 3] , # list of integers that specify the control surface ID's to be used 
                                    #             'controlled_surfaces_coeff': [0, 0, 0, -1], # allows for antisymmetric deflection of aileron
                                    #             'controller_log_route': os.path.join(route, 'output', case_name),
                                    #             'time_history_input_file': os.path.join(route, 'pitch_angle.csv')
                                    #                                  }}
                                                                    }


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

def t_tail_1_cos_gust(route, case_name, set):
    # initial sharpy settings for t_tail HALE span_opt paper
    flow = set['flow']

    flow = ['BeamLoader',
            'AerogridLoader',
            # 'NonLinearStatic',
            # 'StaticUvlm',
            'StaticCoupled',
            'StaticTrim',
            'BeamLoads',
            'AerogridPlot',
            'BeamPlot',
            'DynamicCoupled',
            # 'Modal',
            # 'LinearAssember',
            # 'AsymptoticStability',
            ]

    roll = set['roll']
    alpha = set['alpha']
    beta = set['beta']
    
    tstep_factor = set['tstep_factor']
    tolerance = set['tolerance']
    gravity = set['gravity']
    g_loading = set['g_loading']
    num_cores = set['num_cores']
    dt = set['dt']
    u_inf = set['u_inf']
    rho = set['rho']
    n_step = set['n_step']
    

    fsi_tolerance = set['fsi_tolerance']
    structural_relaxation_factor = set['structural_relaxation_factor']
    n_step = set['n_step']
    relaxation_factor = set['relaxation_factor']

    cs_deflection = set['cs_deflection']
    tail_cs_index = set['tail_cs_index']
    nodes_app = set['nodes_app']
    initial_thrust = set['initial_thrust']
    
    n_tstep = set['n_tstep']

    free_flight = set['free_flight']
                                                     
    # gust_length = set['gust_length']
    # gust_intensity = set['gust_intensity']
    # gust_offset = set['gust_offset']

    # gust settings
    gust_intensity = 0.20
    gust_length = 1 * u_inf
    # gust_offset = 0.5 * u_inf
    gust_offset = 0 

    # numerics
    n_step = 5
    structural_relaxation_factor = 0.6
    relaxation_factor = 0.35
    tolerance = 1e-6
    fsi_tolerance = 1e-4

    num_cores = 2

    file_name = os.path.join(route, '{}.sharpy'.format(case_name))
    settings = dict()
    file_name = route + '/' + case_name + '.sharpy'
    settings = dict()
    settings['SHARPy'] = {'case': case_name,
                          'route': route,
                          'flow': flow,
                          'write_screen': 'on',
                          'write_log': 'on',
                          'log_folder': route + '/output/',
                          'log_file': case_name + '.log'}

    settings['BeamLoader'] = {'unsteady': 'on',
                              'orientation': algebra.euler2quat(np.array([roll,
                                                                          alpha,
                                                                          beta]))}
    settings['AerogridLoader'] = {'unsteady': 'on',
                                  'aligned_grid': 'on',
                                  'mstar': int(20 / tstep_factor),
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

    settings['StaticUvlm'] = {'print_info': 'on',
                              'horseshoe': 'off',
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
                            'initial_thrust': initial_thrust,
                            'fz_tolerance': .1,
                            'fx_tolerance': .1,
                            'm_tolerance': .1,
                            'max_iter': 250}

    settings['NonLinearDynamicCoupledStep'] = {'print_info': 'off',
                                               'max_iterations': 950,
                                               'delta_curved': 1e-1,
                                               'min_delta': tolerance,
                                               'newmark_damp': 5e-3,
                                               'gravity_on': gravity,
                                               'gravity': 9.81*g_loading,
                                               'num_steps': n_tstep,
                                               'dt': dt,
                                               'initial_velocity': u_inf}

    settings['NonLinearDynamicPrescribedStep'] = {'print_info': 'off',
                                                  'max_iterations': 950,
                                                  'delta_curved': 1e-1,
                                                  'min_delta': tolerance,
                                                  'newmark_damp': 5e-3,
                                                  'gravity_on': gravity,
                                                  'gravity': 9.81*g_loading,
                                                  'num_steps': n_tstep,
                                                  'dt': dt,
                                                  'initial_velocity': u_inf * int(free_flight)}

    relative_motion = 'off'
    if not free_flight:
        relative_motion = 'on'
    settings['StepUvlm'] = {'print_info': 'off',
                            'num_cores': num_cores,
                            'convection_scheme': 2,
                            'gamma_dot_filtering': 6,
                            'velocity_field_generator': 'GustVelocityField',
                            'velocity_field_input': {'u_inf': int(not free_flight) * u_inf,
                                                     'u_inf_direction': [1., 0, 0],
                                                     'gust_shape': '1-cos',
                                                     'gust_parameters': {'gust_length': gust_length,
                                                                         'gust_intensity': gust_intensity * u_inf},
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
                         'continuous_eigenvalues': 'off',
                         'dt': dt,
                         'plot_eigenvalues': False}

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
                                                         'remove_inputs': ['u_gust']}
                                   }}

    settings['AsymptoticStability'] = {'print_info': 'on',
                                       'modes_to_plot': [],
                                       'display_root_locus': 'off',
                                       'frequency_cutoff': 0,
                                       'export_eigenvalues': 'off',
                                       'num_evals': 40}

    import configobj
    config = configobj.ConfigObj()
    config.filename = file_name
    for k, v in settings.items():
        config[k] = v
    config.write()

def scaled_si2_dynamic(route, case_name, set):
    # initial sharpy settings for t_tail HALE span_opt paper
    flow = set['flow']
    roll = set['roll']
    alpha = set['alpha']
    beta = set['beta']
    
    tstep_factor = set['tstep_factor']
    tolerance = set['tolerance']
    gravity = set['gravity']
    g_loading = set['g_loading']
    num_cores = set['num_cores']
    dt = set['dt']
    u_inf = set['u_inf']
    rho = set['rho']
    n_step = set['n_step']
    
    fsi_tolerance = set['fsi_tolerance']
    structural_relaxation_factor = set['structural_relaxation_factor']
    n_step = set['n_step']
    relaxation_factor = set['relaxation_factor']

    cs_deflection = set['cs_deflection']
    tail_cs_index = set['tail_cs_index']
    nodes_app = set['nodes_app']
    initial_thrust = set['initial_thrust']
    
    n_tstep = set['n_tstep']

    free_flight = set['free_flight']
                                                     
    gust_length = set['gust_length']
    gust_intensity = set['gust_intensity']
    gust_offset = set['gust_offset']

    file_name = os.path.join(route, '{}.sharpy'.format(case_name))
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
                            'initial_alpha': alpha,
                            'initial_deflection': cs_deflection,
                            'tail_cs_index': tail_cs_index,
                            'thrust_nodes': nodes_app,
                            'initial_thrust': initial_thrust,
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

def static_trim_viscous(route, case_name, set):
    # initial sharpy settings for t_tail HALE span_opt paper
    flow = set['flow']
    roll = set['roll']
    alpha = set['alpha']
    beta = set['beta']
    
    tstep_factor = set['tstep_factor']
    tolerance = set['tolerance']
    gravity = set['gravity']
    g_loading = set['g_loading']
    num_cores = set['num_cores']
    dt = set['dt']
    u_inf = set['u_inf']
    rho = set['rho']
    n_step = set['n_step']
    
    fsi_tolerance = set['fsi_tolerance']
    structural_relaxation_factor = set['structural_relaxation_factor']
    n_step = set['n_step']
    relaxation_factor = set['relaxation_factor']

    cs_deflection = set['cs_deflection']
    tail_cs_index = set['tail_cs_index']
    nodes_app = set['nodes_app']
    thrust = set['initial_thrust']
    
    n_tstep = set['n_tstep']

    free_flight = set['free_flight']

    file_name = os.path.join(route, '{}.sharpy'.format(case_name))
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
                                'max_iter': 200,
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
                            'initial_thrust': thrust,
                            'max_iter': 200,
                            'fz_tolerance': 0.5,
                            'fx_tolerance': 0.5,
                            'm_tolerance': 0.5}

    # settings['BeamOpt'] = {'solver': 'StaticCoupled',
    #                     'solver_settings': settings['StaticCoupled'],
    #                     'Max Strain': 0,
    #                     'Beam_Opt_Nodes': beam_opt_nodes,
    #                     'Beam_Opt_Elems' : beam_opt_elems,
    #                     'Sym': beam_opt_sym,
    #                     'initial_deflection': cs_deflection}

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