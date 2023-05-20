import numpy as np
import os
import h5py as h5
from sharpy.utils.solver_interface import solver, BaseSolver, initialise_solver
import sharpy.utils.settings as settings
import sharpy.linear.src.libss as libss
import scipy.linalg as sclalg
import sharpy.utils.h5utils as h5utils
from sharpy.utils.datastructures import LinearTimeStepInfo
import sharpy.utils.cout_utils as cout
import time
import warnings

@solver
class LinDynamicSim(BaseSolver):
    """Time-domain solution of Linear Time Invariant Systems

    Uses the derived linear time invariant systems and solves it in time domain.

    Requires a ``case_name.lininput.h5`` file in the case root folder that contains:

        * ``x0`` (optional): Initial state vector
        * ``input_vec``: Input vector ``(n_tsteps, n_inputs)``.

    Note:
        This solver is seldom used in SHARPy (its focus is on nonlinear time domain aeroelasticity) hence you may
        find this solver lacking in features. If you use it, you may need to make modifications. We would greatly
        appreciate that you contribute these modifications by means of a pull request!

    """
    solver_id = 'LinDynamicSim'
    solver_classification = 'Coupled'

    settings_types = dict()
    settings_default = dict()
    settings_description = dict()

    settings_types['write_dat'] = 'list(str)'
    settings_default['write_dat'] = []
    settings_description['write_dat'] = 'List of vectors to write: ``x``, ``y``, ``u`` and/or ``t``'

    settings_types['reference_velocity'] = 'float'
    settings_default['reference_velocity'] = 1.
    settings_description['reference_velocity'] = 'Velocity to scale the structural equations when using a non-dimensional system'

    settings_default['n_tsteps'] = 10
    settings_types['n_tsteps'] = 'int'
    settings_description['n_tsteps'] = 'Number of time steps to run'

    settings_types['physical_time'] = 'float'
    settings_default['physical_time'] = 2.
    settings_description['physical_time'] = 'Time to run'

    settings_default['dt'] = 0.001
    settings_types['dt'] = 'float'
    settings_description['dt'] = 'Time increment for the solution of systems without a specified dt'

    settings_types['postprocessors'] = 'list(str)'
    settings_default['postprocessors'] = list()

    settings_types['postprocessors_settings'] = 'dict'
    settings_default['postprocessors_settings'] = dict()

    settings_table = settings.SettingsTable()
    __doc__ += settings_table.generate(settings_types, settings_default, settings_description)

    def __init__(self):

        self.data = None
        self.settings = dict()
        self.postprocessors = dict()
        self.with_postprocessors = False

        self.input_data_dict = dict()
        self.input_file_name = ""

        self.folder = None

    def initialise(self, data, custom_settings=None):

        self.data = data
        if custom_settings:
            self.settings = custom_settings
        else:
            self.settings = data.settings[self.solver_id]
        settings.to_custom_types(self.settings, self.settings_types, self.settings_default)

        # Read initial state and input data and store in dictionary
        self.read_files()

        # Output folder
        self.folder = data.output_folder + '/lindynamicsim/'
        if not os.path.exists(self.folder):
            os.makedirs(self.folder)

        # initialise postprocessors
        self.postprocessors = dict()
        if len(self.settings['postprocessors']) > 0:
            self.with_postprocessors = True
        for postproc in self.settings['postprocessors']:
            self.postprocessors[postproc] = initialise_solver(postproc)
            self.postprocessors[postproc].initialise(
                self.data, self.settings['postprocessors_settings'][postproc], caller=self)

    def run(self):

        ss = self.data.linear.ss

        n_steps = self.settings['n_tsteps']
        x0 = self.input_data_dict.get('x0', np.zeros(ss.states))
        u = self.input_data_dict['u']

        if len(x0) != ss.states:
            warnings.warn('Number of states in the initial state vector not equal to the number of states')
            x0 = np.zeros(ss.states)

        if u.shape[1] != ss.inputs:
            warnings.warn('Dimensions of the input vector not equal to the number of inputs')
            cout.cout_wrap('Number of inputs: %g' % ss.inputs, 3)
            cout.cout_wrap('Number of timesteps: %g' % n_steps, 3)
            cout.cout_wrap('Number of UVLM inputs: %g' % self.data.linear.linear_system.uvlm.ss.inputs, 3)
            cout.cout_wrap('Number of beam inputs: %g' % self.data.linear.linear_system.beam.ss.inputs, 3)
            breakpoint()

        try:
            dt = ss.dt
        except AttributeError:
            dt = self.settings['dt']

        # Total time to run
        T = n_steps*dt

        u_ref = self.settings['reference_velocity']
        # If the system is scaled:
        if u_ref != 1.:
            scaling_factors = self.data.linear.linear_system.uvlm.sys.ScalingFacts
            dt_dimensional = scaling_factors['length'] / u_ref
            T_dimensional = n_steps * dt_dimensional
            T = T_dimensional / scaling_factors['time']
            ss = self.data.linear.linear_system.update(self.settings['reference_velocity'])
        t_dom = np.linspace(0, T, n_steps)

        # Use the scipy linear solver
        sys = libss.ss_to_scipy(ss)
        cout.cout_wrap('Solving linear system using scipy...')
        t0 = time.time()
        # breakpoint()
        out = sys.output(u, t=t_dom, x0=x0)
        ts = time.time() - t0
        cout.cout_wrap('\tSolved in %.2fs' % ts, 1)

        t_out = out[0]
        x_out = out[2]
        y_out = out[1]

        if self.settings['write_dat']:
            cout.cout_wrap('Writing linear simulation output .dat files to %s' % self.folder)
            if 'y' in self.settings['write_dat']:
                np.savetxt(self.folder + '/y_out.dat', y_out)
                cout.cout_wrap('Output vector written', 2)
            if 'x' in self.settings['write_dat']:
                np.savetxt(self.folder + '/x_out.dat', x_out)
                cout.cout_wrap('State vector written', 2)
            if 'u' in self.settings['write_dat']:
                np.savetxt(self.folder + '/u_out.dat', u)
                cout.cout_wrap('Input vector written', 2)
            if 't' in self.settings['write_dat']:
                np.savetxt(self.folder + '/t_out.dat', t_out)
                cout.cout_wrap('Time domain written', 2)
            cout.cout_wrap('Success', 1)

        # Pack state variables into linear timestep info
        cout.cout_wrap('Plotting results...')
        for n in range(len(t_out)-1):
            tstep = LinearTimeStepInfo()
            tstep.x = x_out[n, :]
            tstep.y = y_out[n, :]
            tstep.t = t_out[n]
            tstep.u = u[n, :]
            self.data.linear.timestep_info.append(tstep)
            # TODO: option to save to h5

            # Pack variables into respective aero or structural time step infos (with the + f0 from lin)
            # Need to obtain information from the variables in a similar fashion as done with the database
            # for the beam case

            aero_tstep, struct_tstep = state_to_timestep(self.data, tstep.x, tstep.u, tstep.y)

            self.data.aero.timestep_info.append(aero_tstep)
            self.data.structure.timestep_info.append(struct_tstep)

            # run postprocessors
            if self.with_postprocessors:
                for postproc in self.postprocessors:
                    self.data = self.postprocessors[postproc].run(online=True)

        return self.data

    def read_files(self):

        self.input_file_name = self.data.settings['SHARPy']['route'] + '/' + self.data.settings['SHARPy']['case'] + '.lininput.h5'

        # Check that the file exists
        try:
            h5utils.check_file_exists(self.input_file_name)
            # Read and store
            with h5.File(self.input_file_name, 'r') as input_file_handle:
                self.input_data_dict = h5utils.load_h5_in_dict(input_file_handle)
        except FileNotFoundError:
            pass


def state_to_timestep(data, x, u=None, y=None):
    """
    Warnings:
        Under development

    Writes a state-space vector to SHARPy timesteps

    Args:
        data:
        x:
        u:
        y:

    Returns:

    """

    if data.settings['LinearAssembler']['linear_system_settings']['beam_settings']['modal_projection'] and \
            data.settings['LinearAssembler']['linear_system_settings']['beam_settings']['inout_coords'] == 'modes':
        modal = True
    else:
        modal = False
    # modal = True
    if data.linear.linear_system.uvlm.gust_assembler:
        start_x_aero = data.linear.linear_system.uvlm.gust_assembler.ss_gust.states
    else:
        start_x_aero = 0
    x_aero = x[start_x_aero:data.linear.linear_system.uvlm.ss.states]
    x_struct = x[-data.linear.linear_system.beam.ss.states:]
    # u_aero = TODO: external velocities
    phi = data.linear.linear_system.beam.sys.U
    Kas = data.linear.linear_system.couplings['Kas']

    # Beam output
    y_beam = x_struct

    u_q = np.zeros(data.linear.linear_system.uvlm.ss.inputs)
    if u is not None:
        u_q += u[:data.linear.linear_system.uvlm.ss.inputs]
        u_q[:y_beam.shape[0]] += y_beam
    else:
        u_q[:y_beam.shape[0]] += y_beam

    if modal:
        # add eye matrix for extra inputs
        n_modes = phi.shape[1]
        n_inputs_aero_only = len(u_q) - 2*n_modes  # Inputs to the UVLM other than structural inputs
        u_aero = Kas.dot(sclalg.block_diag(phi, phi, np.eye(n_inputs_aero_only)).dot(u_q))
    else:
        # if u_q.shape[0] !=
        # u_aero_zero = data.linear.tsaero0
        u_aero = Kas.dot(u_q)

    # Unpack input
    zeta, zeta_dot, u_ext = data.linear.linear_system.uvlm.unpack_input_vector(u_aero)

    # Also add the beam forces. I have a feeling there is a minus there as well....
    # Aero
    forces, gamma, gamma_dot, gamma_star = data.linear.linear_system.uvlm.unpack_ss_vector(
        data,
        x_n=x_aero,
        aero_tstep=data.linear.tsaero0,
        track_body=True)

    current_aero_tstep = data.aero.timestep_info[-1].copy()
    current_aero_tstep.forces = [forces[i_surf] + data.linear.tsaero0.forces[i_surf] for i_surf in
                                 range(len(gamma))]
    current_aero_tstep.gamma = [gamma[i_surf] + data.linear.tsaero0.gamma[i_surf] for i_surf in
                                range(len(gamma))]
    current_aero_tstep.gamma_dot = [gamma_dot[i_surf] + data.linear.tsaero0.gamma_dot[i_surf] for i_surf in
                                    range(len(gamma))]
    current_aero_tstep.gamma_star = [gamma_star[i_surf] + data.linear.tsaero0.gamma_star[i_surf] for i_surf in
                                     range(len(gamma))]
    current_aero_tstep.zeta = zeta
    current_aero_tstep.zeta_dot = zeta_dot
    current_aero_tstep.u_ext = u_ext

    # self.data.aero.timestep_info.append(current_aero_tstep)

    aero_forces = data.linear.linear_system.uvlm.C_to_vertex_forces.dot(x_aero)
    beam_forces = data.linear.linear_system.couplings['Ksa'].dot(aero_forces)

    if u is not None:
        u_struct = u[-data.linear.linear_system.beam.ss.inputs:]
    # y_struct = y[:self.data.linear.lsys[sys_id].lsys['LinearBeam'].ss.outputs]

    # Reconstruct the state if modal
    if modal:
        phi = data.linear.linear_system.beam.sys.U
        x_s = sclalg.block_diag(phi, phi).dot(x_struct)
    else:
        x_s = x_struct
    y_s = beam_forces #+ phi.dot(u_struct)
    # y_s = self.data.linear.lsys['LinearBeam'].sys.U.T.dot(y_struct)

    current_struct_step = data.linear.linear_system.beam.unpack_ss_vector(x_s, y_s, data.linear.tsstruct0)
    # data.structure.timestep_info.append(current_struct_step)

    return current_aero_tstep, current_struct_step
