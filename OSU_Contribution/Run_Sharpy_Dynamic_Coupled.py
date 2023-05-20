import pandas as pd
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
from OSU_Contribution.Main_Paramaterized import Paramaterized
import matplotlib.pyplot as plt
import sharpy.sharpy_main as sharpy_main
from OSU_Contribution.utils.SHARPy_Input_Writes import generate_dyn_file



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

# # Flow for Flight Dynamics Modes
# flow = ['BeamLoader',
#         'AerogridLoader',
#         'StaticTrim',
#         'BeamPlot',
#         'AerogridPlot',
#         'AeroForcesCalculator',
#         'DynamicCoupled',
#         'Modal',
#         'LinearAssembler',
#         'AsymptoticStability',
#         ]

# FLIGHT CONDITIONS for AERO SOLVER
u_inf = 23.7 #14 for XHALE, 10 for T_Tail
rho = 0.1736 # 1.225 

# default elevator is -0.020708598 radians 3
alpha =  2.5115*np.pi/180 # 2.2953
roll = 0
beta = 0
cs_deflection =  1.5124*np.pi/180 # 2.1124
thrust = 46.4832 # 37.07

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
    
m = 8
tstep_factor = 1
# dt = 1/m/u_inf*tstep_factor
dt = .05
physical_time = 5
n_tstep = round(physical_time/dt)

num_cores= 4

#numerics
n_step = 8
fsi_tolerance = 1e-3
structural_relaxation_factor = 0.5

relaxation_factor = 0.5
tolerance = 1e-3


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

# Static Trim Settings
nodes_app = run.app_nodes # where the thrust nodes are 
tail_cs_index = run.tail_index # the control surface index that will be used 

global beam_opt_nodes
global beam_opt_sym
#Beam Opt Settings
max_strain = .0001
beam_save_file = 'T_Tail.xlsx'
save_data = True
beam_opt_nodes = run.beam_opt_nodes
beam_opt_elems = run.beam_opt_elems
beam_opt_sym = run.beam_opt_sys
max_iter = 10

# -----------------------------------------------------------------------------------------------------
#----------------------End of Input--------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

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
                            'initial_thrust': thrust,
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

data = sharpy_main.main(['', route + '/' + case_name + '.sharpy'])

# eigenvalues_trim = np.loadtxt('./output/' + case_name + '/stability/eigenvalues.dat')

# fig = plt.figure()
# plt.scatter(eigenvalues_trim[:, 0], eigenvalues_trim[:, 1],
#            marker='x',
#            color='k')
# plt.xlim(-0.5, 0.5)
# plt.ylim(-0.5, 0.5)
# plt.grid()
# plt.xlabel('Real Part, $Re(\lambda)$ [rad/s]')
# plt.ylabel('Imaginary Part, $Im(\lambda)$ [rad/s]');

# grav_forces = data.structure.timestep_info[0].gravity_forces
# grav_z = 0
# for i in range(len(grav_forces)):
#     grav_z += grav_forces[i][2]

grav_forces = data.structure.timestep_info[0].gravity_forces
grav = 0
for i in range(len(grav_forces)):
    grav += np.sqrt(grav_forces[i][0]**2 + grav_forces[i][1]**2 + grav_forces[i][2]**2)
    
masss_temp = abs(grav/(g_loading*9.81)) 

def normal_accel_calc(data, mass):
    total_force_z = []
    acceleration = []
    for i in range(data.ts):
        total_force_z.append(data.structure.timestep_info[i].total_forces[2])
        acceleration.append(data.structure.timestep_info[i].total_forces[2]/mass)

    df1 = pd.DataFrame([total_force_z, acceleration]) 

    beam_save_file = 'Normal_Acceleration_cg.xlsx'

    xlwriter = pd.ExcelWriter(beam_save_file)
    df1.to_excel( xlwriter, header= False,  sheet_name= 'Normal_Accel')
    xlwriter.close()
    
    t_step = range(0, 100)
    plt(t_step, acceleration)
    plt.xlabel('Time Step')
    plt.ylabel('Normal Acceleration (ms^-2)')

normal_accel_calc(data, masss_temp)  

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
    if mass_db[0][0] != 0:
        cg.append([mass_db[1][5]/mass_db[0][0], mass_db[2][3]/mass_db[0][0], mass_db[0][4]/mass_db[0][0]])
    else: 
        cg.append([0, 0, 0])
    
    # right wing 
    if run.ref_lookup[i] == 'right': 
        r.append((elem_start[i] + elem_end[i])/2 + [-1*cg[i][1], cg[i][0], cg[i][2]])
        if str(r[-1][0]) == 'nan':
            print('hi')
        m.append(abs(elem_start[i][1] - elem_end[i][1])*mass_bar[i])
        # if str(m[i]*r[i]) != 'nan': 
        mr.append(m[i]*r[i])
    
    # left wing
    if run.ref_lookup[i] == 'left': 
        r.append((elem_start[i] + elem_end[i])/2 + [1*cg[i][1], cg[i][0], cg[i][2]])
        if str(r[-1][0]) == 'nan':
            print('hi')
        m.append(abs(elem_start[i][1] - elem_end[i][1])*mass_bar[i])
        # if str(m[i]*r[i]) != 'nan': 
        mr.append(m[i]*r[i])
    
    # v_fin 
    if run.ref_lookup[i] == 'v_fin': 
        r.append((elem_start[i] + elem_end[i])/2 + [-1*cg[i][1], cg[i][0], cg[i][2]])
        if str(r[-1][0]) == 'nan':
            print('hi')
        m.append(abs(elem_start[i][2] - elem_end[i][2])*mass_bar[i])
        # if str((m[i]*r[i])[0]) != 'nan': 
        mr.append(m[i]*r[i])
    
    # fuselage assumed that mass is located on the Elastic Axis
    if run.ref_lookup[i] == 'fuselage': 
        r.append((elem_start[i] + elem_end[i])/2)
        if str(r[-1][0]) == 'nan':
            print('hi')
        m.append(abs(elem_start[i][0] - elem_end[i][0])*mass_bar[i])
        # if str(m[i]*r[i]) != 'nan': 
        mr.append(m[i]*r[i])

    # fuselage_neg assumed that mass is located on the Elastic Axis
    if run.ref_lookup[i] == 'fuselage_neg': 
        r.append((elem_start[i] + elem_end[i])/2)
        if str(r[-1][0]) == 'nan':
            print('hi') 
        m.append(abs(elem_start[i][0] - elem_end[i][0])*mass_bar[i])
        # if str(m[i]*r[i]) != 'nan': 
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
# print('\nMass From Applied Forces: ' + str(abs(grav_z/(9.81*g_loading)))+ '\n')
print('Mass found by multiplying m_bar by length for each section: ' + str(m_total[0]) + '\n')
print('x_cg = ' + str(cg_total[0]) + ' y_cg = ' + str(cg_total[1]) + ' z_cg = ' + str(cg_total[2]) + '\n')

# LONGITUDINAL NEUTRAL POINT APPROXIMATION 
def long_neutral_approx():
    ar = 21.1 #AR wing
    ar_t = 7.56 # AR of H-Tail
    s_w = 162 # surface area of wing 
    b = 58.5 # wing span
    s_t = 12.6 # surface area of tail
    c_bar =  s_w/b # average chord 
    c_mac = 2.6 # mean aerodynamic chord of wing
    
    cl_alpha_0 = 6.01 # 2D lift vs alpha relation of wing. Obtained from drag polar info
    cl_alpha_0_t = 6.01 # 2D lift vs alpha relation of tail. Obtained from drag polar info 
    
    # WING TO TAIL INPUTS
    w_0 = 6.51 # x_cordinate elastic axis loction of wing
    ea_w = .383 # % chord relationship with elastic axis. Distance between EA and leading edge of wing
    c_w = 2.813 # chord of the wing
    
    t_0 = 17.65 # x_cordinate elastic axis loction of tail
    ea_t = .346 # % chord relationship with elastic axis. Distance between EA and leading edge of tail
    c_t = 1.3 # chord of the tail
    
    eta = 1 # HTail efficiency factor 
    
    def wing_correction(a0, ar):
        return a0/ (1 + a0/(np.pi*ar))

    def vh(s_t, l, s_w, c_mac):
        return (s_t*l)/(s_w*c_mac) 
    
    def wing_tail_dis(w_0, ea_w, c_w, t_0, ea_t, c_t):
        w_loc = w_0 - c_w*ea_w + .25*c_w
        t_loc = t_0 - c_t*ea_t + .25*c_t
        return t_loc - w_loc

    def neutral(c_bar, eta, vh, cl_alpha, cl_alpha_t, ar, w_0, c_w, ea_w):
        return (c_bar * (.25 - 0 + eta*vh*(cl_alpha_t/cl_alpha)*(1 - (2*cl_alpha)/(np.pi*ar)))) + (w_0 - c_w*ea_w)
    
    cl_alpha = wing_correction(cl_alpha_0, ar)
    cl_alpha_t = wing_correction(cl_alpha_0_t, ar_t)
    l = wing_tail_dis(w_0, ea_w, c_w, t_0, ea_t, c_t)
    vh_out = vh(s_t, l, s_w, c_mac)
    neutral_out = neutral(c_bar, eta, vh_out, cl_alpha, cl_alpha_t, ar, w_0, c_w, ea_w)

    print('Nuetral Point Estimation = ' + str(neutral_out) + '/n')

# long_neutral_approx()

# Writes save file for BeamOpt
if save_data == True:
    for i in range(len(flow)):
        if flow[i] == 'BeamOpt':
            massss = abs(grav_z/(9.81*g_loading)) + sum(data.structure.lumped_mass) 
            df = pd.DataFrame([max_strain, massss[0]])

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
                    
            df.to_excel(beam_save_file, index=['Strain', 'Mass'], header= False)
            break 



