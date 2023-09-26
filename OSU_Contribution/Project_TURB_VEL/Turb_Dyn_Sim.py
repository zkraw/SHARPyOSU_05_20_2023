import numpy as np
import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import json

 
# sys.path is a list of absolute path strings
current = os.path.abspath(__file__)
sharpy_directory = os.path.abspath(os.path.join(current, '..', '..', '..'))
sys.path.append(sharpy_directory)

import sharpy.utils.algebra as algebra
from OSU_Contribution.Main_Paramaterized import Paramaterized
import sharpy.sharpy_main as sharpy_main
from OSU_Contribution.utils.SHARPy_Input_Writes import generate_dyn_file, generate_dyn_file_2


home = str(Path.home())
route = os.path.dirname(os.path.realpath(__file__))
# xmf_file = os.path.join(route, "Velocity Field", "test.xmf") 
# xmf_file = os.path.join(route, "Velocity Field", "Python Creation", "1_cos.xdmf") 
# xmf_file = os.path.join(route, "Velocity Field", "Python Creation", "Reversed_Fine_h5_write", 'RUN_Reversed_Fine.xdmf') 
# xmf_file = os.path.join(route, "Velocity Field", "Python Creation", "Iterp_Validation_h5_write", 'RUN_Iterp_Validation.xdmf') 

# xmf_file = os.path.join(route, "Velocity Field", "Python Creation", "VTK_resampled_fine_h5_write", 'RUN_VTK_resampled_fine.xdmf') 

xmf_file = os.path.join(route, "Velocity Field", "Python Creation", "const_8_x_longer_h5_write", 'RUN_const_8_x_longer.xdmf') 
# xmf_file = os.path.join(route, "Velocity Field", "Python Creation", "const_3_z_h5_write", 'RUN_const_3_z_Assymmetric.xdmf') 
# xmf_file = os.path.join(route, "Velocity Field", "Python Creation", "const_3_z_h5_write", 'RUN_const_3_z_Symmetric.xdmf') 
# -----------------------------------------------------------------------------------------------------
#----------------------INPUTS--------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
plot_file_loc = os.path.join(route, "Sim_Plots_and_Results") 

# file = 'T_Tail_Stiff_Flat linassem.xlsx'
file = 'Urban Air Mobility.xlsx'

global case_name
case_name = 'new_case'
excel = True 

input_file = os.path.dirname(os.path.abspath(__file__)) + '/Input Files/' + file 

from OSU_Contribution.utils.OPT_functions.user_setting import t_tail
set, case_name, spar_cap_dict = t_tail(case_name)
set['xmf'] = xmf_file

set['flow'] = ['BeamLoader',
        'AerogridLoader',
        # 'NonLinearStatic',
        # 'StaticUvlm',
        'StaticCoupled',
        'StaticTrim',
        'BeamLoads',
        'AerogridPlot',
        'BeamPlot',
        # 'DynamicCoupled',
        # 'Modal',
        # 'LinearAssembler',
        # 'AsymptoticStability',
         ]

# saving values to new file
#Paramaterized.excel_write(file, 'new_file.xlsx', ['start_Y'], [[3, 4, 5, 6, 7, 8]], [[123456789, 100, 100, 100, 100, 100]])

# generate FEM and Aero from input file
run = Paramaterized(case_name, route, input_file, set['xfoil'], excel)

set['nodes_app'] = run.app_nodes # where the thrust nodes are 
set['tail_cs_index'] = run.tail_index # the control surface index that will be used 
  
# OPTIMIZATION SETTINGS
save_data = True
set['beam_opt_nodes'] = run.beam_opt_nodes
set['beam_opt_elems'] = run.beam_opt_elems
set['beam_opt_sym'] = run.beam_opt_sys

set['u_inf'] = 79.74 
set['rho'] = 1.225

set['alpha'] =  3.8865*np.pi/180
set['cs_deflection'] = 9.4579*np.pi/180
set['initial_thrust'] = 1035.0551

# T_Tail with dihedral
# set['alpha'] =  3.8732*np.pi/180
# set['cs_deflection'] = -.6596*np.pi/180
# set['initial_thrust'] = 4.6123 

set['m_star'] = 20

# Time for dynamic solver
m = 4 # chordwise panels  
set['tstep_factor'] = 1
# physical_time = 10 
physical_time = 30 
# set['dt'] = .2
set['dt'] = 1.0 / m / set['u_inf'] * set['tstep_factor']
set['n_tstep'] = round(physical_time/set['dt'])

# offset to move aircraft to proper positions in wind field
set['offset'] = np.array([0, 0, 0])

generate_dyn_file_2([set['dt'], set['n_tstep'], route, case_name, set['amplitude'], set['period'], set['free_flight']], run.num_nodes)

from OSU_Contribution.utils.OPT_functions.init_sharpy_settings import static_trim_viscous
static_trim_viscous(route, case_name, set) 

# from OSU_Contribution.utils.OPT_functions.init_sharpy_settings import turbulent_velocity_feild_controller
# turbulent_velocity_feild_controller(route, case_name, set) 

plot_sim_flag = False 

# from OSU_Contribution.utils.OPT_functions.init_sharpy_settings import linearassembler
# linearassembler(route, case_name, set) 

# from OSU_Contribution.utils.OPT_functions.init_sharpy_settings import scaled_si2_dynamic
# scaled_si2_dynamic(route, case_name, set) 

turb_sim = sharpy_main.main(['', route + '/' + case_name + '.sharpy'])
test = turb_sim.structure.timestep_info[0].cga()

def plot_sim(data):
        case_data = data
        n_tsteps = len(case_data.structure.timestep_info)
#-----------------------------------------------------------------------------------
#--------------------------TRAJECTORY--------------------------------------------------
#-----------------------------------------------------------------------------------
        fig, ax = plt.subplots(1, 1, figsize=(7, 4))
        # extract information
        xyz = np.zeros((n_tsteps, 3))
        for it in range(n_tsteps):
                xyz[it, 0] = -case_data.structure.timestep_info[it].for_pos[0] # the - is so that increasing time -> increasing x
                xyz[it, 1] = case_data.structure.timestep_info[it].for_pos[1]
                xyz[it, 2] = case_data.structure.timestep_info[it].for_pos[2]

        ax.plot(xyz[:, 0], xyz[:, 2])
        fig.suptitle('Longitudinal Trajectory')
        ax.set_xlabel('X [m]')
        ax.set_ylabel('Z [m]')

        # Set the limits for the y-axis
        val = 2 
        ax.set_ylim(np.min(xyz[:, 2])- val, np.max(xyz[:,2]) + val)


        plt.savefig(os.path.join(plot_file_loc, 'Long_Trajectory__' + case_name + '.png'))

        fig, ax = plt.subplots(1, 1, figsize=(7, 6), sharex=True)
        ax.plot(xyz[:, 1], xyz[:, 0])
        fig.suptitle('Lateral Trajectory')
        ax.set_xlabel('Y [m]')
        ax.set_ylabel('X [m]')

        # Set the limits for the y-axis
        val = 2 
        ax.set_xlim(np.min(xyz[:, 1])- val, np.max(xyz[:,1]) + val)


        plt.savefig(os.path.join(plot_file_loc, 'Lat_Trajectory__' + case_name + '.png'))
#-----------------------------------------------------------------------------------
#--------------------------VELOCITIES-----------------------------------------------
#-----------------------------------------------------------------------------------
        fig, ax = plt.subplots(3, 1, figsize=(7, 6), sharex=True)
        ylabels = ['Vx [m/s]', 'Vy [m/s]', 'Vz [m/s]']

        # extract information
        dt = case_data.settings['DynamicCoupled']['dt']
        time_vec = np.linspace(0, n_tsteps*dt, n_tsteps)
        for_vel_trans = np.zeros((n_tsteps, 3))
        for it in range(n_tsteps):
                for_vel_trans[it, 0:3] = case_data.structure.timestep_info[it].for_vel[0:3]

        for idim in range(3):
                ax[idim].plot(time_vec, for_vel_trans[:, idim])
                ax[idim].set_ylabel(ylabels[idim])

        ax[2].set_xlabel('time [s]')
        plt.subplots_adjust(hspace=0)
        fig.suptitle(case_name + 'Linear RBM Velocities');
        # ax.set_xlabel('X [m]')
        # ax.set_ylabel('Z [m]');
        # Set the limits for the y-axis
        val = 2 
        ax[0].set_ylim(np.min(for_vel_trans[:, 0])- val, np.max(for_vel_trans[:,0]) + val)
        ax[1].set_ylim(np.min(for_vel_trans[:, 1])- val, np.max(for_vel_trans[:,1]) + val)
        ax[2].set_ylim(np.min(for_vel_trans[:, 2])- val, np.max(for_vel_trans[:,2]) + val)
        plt.savefig(os.path.join(plot_file_loc, 'Velocities__' + case_name + '.png'))

#-----------------------------------------------------------------------------------
#--------------------------ANGULAR VELOCITIES---------------------------------------
#-----------------------------------------------------------------------------------
        fig, ax = plt.subplots(3, 1, figsize=(7, 6), sharex=True)
        ylabels = ['Roll rate [deg/s]', 'Pitch rate [deg/s]', 'Yaw rate [deg/s]']

        # extract information
        n_tsteps = len(case_data.structure.timestep_info)
        time_vec = np.linspace(0, (n_tsteps-1)*dt, n_tsteps)
        for_vel = np.zeros((n_tsteps, 3))
        for it in range(n_tsteps):
                for_vel[it, 0:3] = case_data.structure.timestep_info[it].for_vel[3:6]*180/np.pi

        for idim in range(3):
                ax[idim].plot(time_vec, for_vel[:, idim])
                ax[idim].set_ylabel(ylabels[idim])

        ax[2].set_xlabel('time [s]')
        plt.subplots_adjust(hspace=0)
        fig.suptitle(case_name + 'Angular Velocities');
        # ax.set_xlabel('X [m]')
        # ax.set_ylabel('Z [m]');
        # Set the limits for the y-axis
        val = 2 
        ax[0].set_ylim(np.min(for_vel[:, 0])- val, np.max(for_vel[:,0]) + val)
        ax[1].set_ylim(np.min(for_vel[:, 1])- val, np.max(for_vel[:,1]) + val)
        ax[2].set_ylim(np.min(for_vel[:, 2])- val, np.max(for_vel[:,2]) + val)
        plt.savefig(os.path.join(plot_file_loc,  'Angular_Vel__' + case_name + '.png')) 

#-----------------------------------------------------------------------------------
#--------------------------U_EXTERNAL-----------------------------------------------
#-----------------------------------------------------------------------------------
        u_ext_wing = []
        u_ext_other_wing = []
        u_ext_tail = []
        for it in range(n_tsteps):
                # u_temp = []
                # for ii in range(len(case_data.aero.timestep_info[it].u_ext)):
                u_ext_wing.append(case_data.aero.timestep_info[it].u_ext[0])
                u_ext_other_wing.append(case_data.aero.timestep_info[it].u_ext[1])
                u_ext_tail.append(case_data.aero.timestep_info[it].u_ext[4])
                # u_temp.append(case_data.aero.timestep_info[it].u_ext[ii])
                # if it == len(case_data.aero.timestep_info[it].u_ext) - 1:
                # u_temp = np.array(u_temp) 
                if it == n_tsteps - 1:
                        u_ext_other_wing = np.array(u_ext_other_wing)
                        u_ext_tail = np.array(u_ext_tail)
        
        coord = {}                
        for ii in range(len(case_data.structure.elements)):
                coord['{}'.format(ii)] = case_data.structure.elements[0].coordinates_ini 

        # Convert the dictionary into a NumPy array
        coord_array = np.array(list(coord.values()))


        # u_ext = np.array(u_ext)
        # df = pd.DataFrame({'xz': xz, 'time': time_vec, 'for_vel': for_vel_trans, 'for_ang_vel': for_vel, 'u_ext': u_ext})
        np.savez(os.path.join(plot_file_loc, case_name+ '_output'), xyz=xyz, time = time_vec, for_vel = for_vel_trans, for_vel_ang = for_vel, 
                 u_ext_other_wing = u_ext_other_wing, u_ext_tail = u_ext_tail, u_ext_wing = u_ext_wing, coord_array = coord_array)

if plot_sim_flag:
    plot_sim(turb_sim)

