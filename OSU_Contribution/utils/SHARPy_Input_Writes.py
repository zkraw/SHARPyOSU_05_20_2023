
import h5py as h5
import numpy as np
import os

def generate_dyn_file(input_stuff, num_node_input):
    dt = input_stuff[0]
    n_tstep= input_stuff[1]
    route= input_stuff[2]
    case_name= input_stuff[3]
    num_node = num_node_input
    amplitude = input_stuff[4]
    period = input_stuff[5]
    free_flight = input_stuff[6]

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

def generate_dyn_file_2(input_stuff, num_node_input):
    dt = input_stuff[0]
    n_tstep= input_stuff[1]
    route= input_stuff[2]
    case_name= input_stuff[3]
    num_node = num_node_input
    amplitude = input_stuff[4]
    period = input_stuff[5]
    free_flight = input_stuff[6]

    dynamic_forces_time = None
    with_dynamic_forces = False
    with_forced_vel = False
    if not free_flight:
        with_forced_vel = True

    if with_dynamic_forces:
        f1 = 100
        dynamic_forces = np.zeros((num_node, 6))
        app_node = [int(num_node_main - 1), int(num_node_main)]
        dynamic_forces[app_node, 2] = f1
        force_time = np.zeros((n_tstep,))
        limit = round(0.05 / dt)
        force_time[50:61] = 1

        dynamic_forces_time = np.zeros((n_tstep, num_node, 6))
        for it in range(n_tstep):
            dynamic_forces_time[it, :, :] = force_time[it] * dynamic_forces

    forced_for_vel = None
    if with_forced_vel:
        forced_for_vel = np.zeros((n_tstep, 6))
        forced_for_acc = np.zeros((n_tstep, 6))
        for it in range(n_tstep):
            # if dt*it < period:
            # forced_for_vel[it, 2] = 2*np.pi/period*amplitude*np.sin(2*np.pi*dt*it/period)
            # forced_for_acc[it, 2] = (2*np.pi/period)**2*amplitude*np.cos(2*np.pi*dt*it/period)

            forced_for_vel[it, 3] = 2 * np.pi / period * amplitude * np.sin(2 * np.pi * dt * it / period)
            forced_for_acc[it, 3] = (2 * np.pi / period) ** 2 * amplitude * np.cos(2 * np.pi * dt * it / period)

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

def clean_test_files(route, case_name):
    fem_file_name = os.path.join(route, '{}.fem.h5'.format(case_name)) 
    if os.path.isfile(fem_file_name):
        os.remove(fem_file_name)

    aero_file_name = os.path.join(route, '{}.aero.h5'.format(case_name)) 
    if os.path.isfile(aero_file_name):
        os.remove(aero_file_name)
        
def write_aero(aero, route, case_name):
    
    with h5.File(os.path.join(route, '{}.aero.h5'.format(case_name)), 'a') as h5file:
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

def write_fem(fem, route, case_name):
    
    with h5.File(os.path.join(route, '{}.fem.h5'.format(case_name)), 'a') as h5file:
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