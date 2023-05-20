import numpy as np

def t_tail(case_name):
    
    set = {}
    
    # FLIGHT CONDITIONS for AERO SOLVER
    set['u_inf'] = 15 #12 #14 for XHALE, 10 for T_Tail
    set['rho'] = 1.225 # 1.225 

    # default elevator is -0.020708598 radians 3
    set['alpha'] = 2.4724*np.pi/180 # 2.2953
    set['roll'] = 0
    set['beta'] = 0

    set['cs_deflection'] = -1.2208*np.pi/180 # 2.1124
    set['initial_thrust'] = 75 

    set['gravity'] = 'on'
    set['g_loading'] = 1 # acceleration of gravity (9.81) * factor 

    # gust settings
    set['gust_intensity'] =  0.18*set['u_inf']
    set['gust_length'] = set['u_inf'] / .307
    set['gust_offset'] = 0.1*set['u_inf']

    #DYNAMIC FLIGHT SETTINGS
    set['free_flight'] = True
    set['amplitude'] = 0
    set['period'] = 0
    if not set['free_flight']:
        case_name += '_prescribed'
        amplitude = 0*np.pi/180
        period = 3
        case_name += '_amp_' + str(amplitude).replace('.', '') + '_period_' + str(period)

    # Time for dynamic solver    
    physical_time = .1
    set['dt'] = .05
    # dt = 1/m/u_inf*tstep_factor
    set['n_tstep'] = round(physical_time/set['dt'])

    set['m'] = 8
    set['tstep_factor'] = 1

    set['num_cores']= 4
    # n_step = 100

    #numerics
    set['n_step'] = 5
    set['fsi_tolerance'] = 1e-3
    set['structural_relaxation_factor'] = 0.3
    set['relaxation_factor'] = 0.5
    set['tolerance'] = 1e-2

    set['num_cores'] = 4

    # XFIOL INPUTS

    a = 343         #speed of sound (m/s)
    u = 1.802e-05   #dynamic viscosity (kg/(m*s))

    set['xfoil'] = {
        'Mach': set['u_inf']/a,
        'n_iter': 100,
        'alfa_1': -16,
        'alfa_2': 16,
        'alfa_step': .5, 
        'u_inf': set['u_inf'],
        'rho': set['rho'],
        'u': u,
    }

    spar_cap_dict = {}
    spar_cap_dict['E_cap'] = 1.9e11
    spar_cap_dict['E_web'] = 1.3e10
    spar_cap_dict['G'] = 4.8e10

    spar_cap_dict['h_spar'] = 0.129
    spar_cap_dict['b_spar'] = 0.2

    spar_cap_dict['density'] = 1540 #163.934

    spar_cap_dict['mass_other'] = .5 # ratio of the mass that comes from sources other than the spar 
    # margin of safety calculations
    spar_cap_dict['F_tu'] = 1.35e9
    spar_cap_dict['F_cu'] = 4.85e8
    spar_cap_dict['F_tu_web'] = 7e7
    spar_cap_dict['F_su'] = 2.5e8
    
    spar_cap_dict['thickness'] = .001 # thickness of carbon not including core 
    spar_cap_dict['core'] = .005 # thickness of carbon not including core

    spar_cap_dict['core_density'] = 29 #kg/m^3
    return set, case_name, spar_cap_dict

def scaled_SI2(case_name, span):

    set = {}
    
    # FLIGHT CONDITIONS for AERO SOLVER
    set['u_inf'] = 23.7 #14 for XHALE, 10 for T_Tail
    set['rho'] = 0.1736 # 1.225 

    # default elevator is -0.020708598 radians 3
    set['alpha'] = 4.1803*np.pi/180 # 2.2953
    set['roll'] = 0
    set['beta'] = 0

    set['cs_deflection'] = 4.5411*np.pi/180 # 2.1124
    set['initial_thrust'] = 64.8168 # 37.07

    set['gravity'] = 'on'
    set['g_loading'] = 1 # acceleration of gravity (9.81) * factor 

    # gust settings
    set['gust_intensity'] =  0.18*set['u_inf']
    set['gust_length'] = set['u_inf'] / .307
    set['gust_offset'] = 0.1*set['u_inf']

    #DYNAMIC FLIGHT SETTINGS
    set['free_flight'] = True
    set['amplitude'] = 0
    set['period'] = 0
    if not set['free_flight']:
        case_name += '_prescribed'
        amplitude = 0*np.pi/180
        period = 3
        case_name += '_amp_' + str(amplitude).replace('.', '') + '_period_' + str(period)

    # Time for dynamic solver    
    physical_time = 5
    set['dt'] = .05
    # dt = 1/m/u_inf*tstep_factor
    set['n_tstep'] = round(physical_time/set['dt'])

    set['m'] = 8
    set['tstep_factor'] = 1

    set['num_cores']= 4
    # n_step = 100

    #numerics
    set['n_step'] = 5
    set['fsi_tolerance'] = 1e-3
    set['structural_relaxation_factor'] = 0.3
    set['relaxation_factor'] = 0.5
    set['tolerance'] = 1e-2

    set['num_cores'] = 4

    # XFIOL INPUTS

    a = 343         #speed of sound (m/s)
    u = 1.802e-05   #dynamic viscosity (kg/(m*s))

    set['xfoil'] = {
        'Mach': set['u_inf']/a,
        'n_iter': 100,
        'alfa_1': -16,
        'alfa_2': 16,
        'alfa_step': .5, 
        'u_inf': set['u_inf'],
        'rho': set['rho'],
        'u': u,
    }

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

    return set, case_name, spar_cap_dict