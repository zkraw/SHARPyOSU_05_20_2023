import os
import numpy as np
import sharpy.utils.algebra as algebra
import h5py as h5

def nonlinear_static(route, case_name, logs, trim_set):
    file = open(route + '/' + case_name + '.sharpy')
    lines = file.readlines()
    nu_lines = []
    for i in lines:
        current = i.replace(' ', '').split('=')
        if current[0] == 'flow':
            current[1] = 'BeamLoader, AerogridLoader, NonLinearStatic, BeamLoads, AerogridPlot, BeamPlot' + '\n'
            nu_lines.append(' = '.join(current))
            continue

        if current[0] == 'gravity':
            current[1] = str(1*9.81) + '\n'
            nu_lines.append(' = '.join(current))
            continue

        nu_lines.append(i)
    
    # making new SHARPy run file
    os.remove(route + '/' + case_name + '.sharpy')
    
    with open(route + '/' + case_name + '.sharpy', 'w') as f:
        for line in nu_lines:
            f.write(line)
    return 1

def one_g_trim(route, case_name, logs, trim_set):
    file = open(route + '/' + case_name + '.sharpy')
    lines = file.readlines()
    nu_lines = []
    for i in lines:
        current = i.replace(' ', '').split('=')
        if current[0] == 'flow':
            
            # current[1] = 'BeamLoader, AerogridLoader, StaticCoupled, BeamLoads, AerogridPlot' + '\n'
            current[1] = 'BeamLoader, AerogridLoader, StaticTrim, BeamLoads, AerogridPlot, BeamPlot' + '\n'
            nu_lines.append(' = '.join(current))
            continue

        if current[0] == 'gravity':
            current[1] = str(1*9.81) + '\n'
            nu_lines.append(' = '.join(current))
            continue

        if current[0] == 'initial_alpha': 
            if logs['alpha_log'] != []:
                current[1] = str(logs['alpha_log'][-1]*np.pi/180) + '\n'
                nu_lines.append(' = '.join(current))
            else:
                current[1] = str(trim_set['alpha']) + '\n'
                nu_lines.append(' = '.join(current))    
            continue

        if current[0] == 'initial_deflection':
            if logs['elevator_log'] != []:
                current[1] = str(logs['elevator_log'][-1]*np.pi/180) + '\n'
                nu_lines.append(' = '.join(current))
            else:
                current[1] = str(trim_set['cs_deflection']) + '\n'
                nu_lines.append(' = '.join(current))   
            continue

        if current[0] == 'initial_thrust':
            if logs['thrust_log'] != []:
                current[1] = str(logs['thrust_log'][-1]) + '\n'
                nu_lines.append(' = '.join(current))
            else:
                current[1] = str(trim_set['initial_thrust']) + '\n'
                nu_lines.append(' = '.join(current))   
            continue
        
        nu_lines.append(i) 

    os.remove(route + '/' + case_name + '.sharpy')
    
    with open(route + '/' + case_name + '.sharpy', 'w') as f:
        for line in nu_lines:
            f.write(line)
            
    return 1 # specifies g_loading so that mass can be calculated accurately

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
            
def two_g_trim(route, case_name, logs, trim_set):
    # lines of code are needed for twist minimzation as 2g and 1g are ran together
    # this should conserve cross compatability without having to specify additiona parameter
    try:
        a_log = logs['alpha_log_2g']
        e_log = logs['elevator_log_2g']
        t_log = logs['thrust_log_2g']
    except:
        a_log = logs['alpha_log']
        e_log = logs['elevator_log']
        t_log = logs['thrust_log']
        
    file = open(route + '/' + case_name + '.sharpy')
    lines = file.readlines()
    nu_lines = []
    for i in lines:
        current = i.replace(' ', '').split('=')
        if current[0] == 'flow':
            
            # current[1] = 'BeamLoader, AerogridLoader, StaticCoupled, BeamLoads, AerogridPlot' + '\n'
            current[1] = 'BeamLoader, AerogridLoader, StaticTrim, BeamLoads, AerogridPlot, BeamPlot' + '\n'
            nu_lines.append(' = '.join(current))
            continue

        if current[0] == 'gravity':
            current[1] = str(2*9.81) + '\n'
            nu_lines.append(' = '.join(current))
            continue

        if current[0] == 'initial_alpha': 
            if a_log != []:
                current[1] = str(a_log[-1]*np.pi/180) + '\n'
                nu_lines.append(' = '.join(current))
            else:
                current[1] = str(trim_set['alpha']) + '\n'
                nu_lines.append(' = '.join(current))    
            continue

        if current[0] == 'initial_deflection':
            if e_log != []:
                current[1] = str(e_log[-1]*np.pi/180) + '\n'
                nu_lines.append(' = '.join(current))
            else:
                current[1] = str(trim_set['cs_deflection']) + '\n'
                nu_lines.append(' = '.join(current))   
            continue

        if current[0] == 'initial_thrust':
            if t_log != []:
                current[1] = str(t_log[-1]) + '\n'
                nu_lines.append(' = '.join(current))
            else:
                current[1] = str(trim_set['initial_thrust']) + '\n'
                nu_lines.append(' = '.join(current))   
            continue
        
        nu_lines.append(i) 

    os.remove(route + '/' + case_name + '.sharpy')
    
    with open(route + '/' + case_name + '.sharpy', 'w') as f:
        for line in nu_lines:
            f.write(line)
            
    return 2 # specifies g_loading so that mass can be calculated accurately

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
            
def turbulent_gust(route, case_name, logs, trim_set, fem):
    file = open(route + '/' + case_name + '.sharpy')
    lines = file.readlines()
    nu_lines = []
    write_lines = []
    for i in lines:
        write_lines.append(i)
        current = i.replace(' ', '').split('=')
        if current[0] == 'flow':
            current[1] = 'BeamLoader, AerogridLoader, StaticCoupled, BeamLoads, AerogridPlot, BeamPlot, DynamicCoupled' + '\n'
            nu_lines.append(' = '.join(current))
            del write_lines[-1]
            write_lines.append(nu_lines[-1])
            continue

        if current[0] == 'gravity':
            current[1] = str(1*9.81) + '\n'
            nu_lines.append(' = '.join(current))
            del write_lines[-1]
            write_lines.append(nu_lines[-1])
            continue

        if current[0] == 'orientation':
            if logs['alpha_log'] != []:
                current[1] = str(algebra.euler2quat(np.array([trim_set['roll'], float(logs['alpha_log'][-1]*np.pi/180), trim_set['beta']]))) + '\n'
                nu_lines.append(' = '.join(current))
                del write_lines[-1]
                write_lines.append(nu_lines[-1])
            else:
                current[1] = str(algebra.euler2quat(np.array([trim_set['roll'], float(trim_set['alpha']*np.pi/180), trim_set['beta']]))) + '\n'
                nu_lines.append(' = '.join(current))
                del write_lines[-1]
                write_lines.append(nu_lines[-1])
            continue
        
        nu_lines.append(i) 

    os.remove(route + '/' + case_name + '.sharpy')
    
    with open(route + '/' + case_name + '.sharpy', 'w') as f:
        for line in write_lines:
            f.write(line)
            
    #------REWRITING ELEVATOR CONTROL INPUT TO TRIM VALUE----------------------------
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

    # creates new aero file with updated elevator angle
    len(logs['elevator_log'])
    if logs['elevator_log'] != []:
        aero['control_surface_deflection'][trim_set['tail_cs_index']] = float(logs['elevator_log'][-1]*np.pi/180)
    else:
        aero['control_surface_deflection'][trim_set['tail_cs_index']] = float(trim_set['cs_deflection']*np.pi/180)

     # creates new fem file with updated thrust values 
    for iii in trim_set['nodes_app']:
        sign = np.sign([fem['app_forces'][iii][1]])
        if logs['thrust_log'] != []:
            fem['app_forces'][iii][1] = float(sign[0]*(logs['thrust_log'][-1])) 
        else:
            fem['app_forces'][iii][1] = float(sign[0]*(trim_set['initial_thrust'])) 
    
    return aero, fem

def bad_g_trim(route, case_name, logs, trim_set):
    file = open(route + '/' + case_name + '.sharpy')
    lines = file.readlines()
    nu_lines = []
    for i in lines:
        current = i.replace(' ', '').split('=')
        if current[0] == 'flow':
            
            # current[1] = 'BeamLoader, AerogridLoader, StaticCoupled, BeamLoads, AerogridPlot' + '\n'
            current[1] = 'BeamLoader, AerogridLoader, StaticTrim, BeamLoads, AerogridPlot, BeamPlot' + '\n'
            nu_lines.append(' = '.join(current))
            continue

        if current[0] == 'gravity':
            current[1] = str(2*9.81) + '\n'
            nu_lines.append(' = '.join(current))
            continue

        if current[0] == 'initial_alpha': 
            if logs['alpha_log'] != []:
                current[1] = str(logs['alpha_log'][-1]*np.pi/180) + '\n'
                nu_lines.append(' = '.join(current))
            else:
                current[1] = str(trim_set['alpha']) + '\n'
                nu_lines.append(' = '.join(current))    
            continue

        if current[0] == 'initial_deflection':
            if logs['elevator_log'] != []:
                current[1] = str(logs['elevator_log'][-1]*np.pi/180) + '\n'
                nu_lines.append(' = '.join(current))
            else:
                current[1] = str(trim_set['cs_deflection']) + '\n'
                nu_lines.append(' = '.join(current))   
            continue

        if current[0] == 'initial_thrust':
            if logs['thrust_log'] != []:
                current[1] = str(logs['thrust_log'][-1]) + '\n'
                nu_lines.append(' = '.join(current))
            else:
                current[1] = str(trim_set['initial_thrust']) + '\n'
                nu_lines.append(' = '.join(current))   
            continue
        
        nu_lines.append(i) 

    os.remove(route + '/' + case_name + '.sharpy')
    
    with open(route + '/' + case_name + '.sharpy', 'w') as f:
        for line in nu_lines:
            f.write(line)
            
    return 2 # specifies g_loading so that mass can be calculated accurately

def any_g_trim(route, case_name, logs, trim_set):
    file = open(route + '/' + case_name + '.sharpy')
    lines = file.readlines()
    write_lines = []
    nu_lines = []
    for i in lines:
        write_lines.append(i)
        current = i.replace(' ', '').split('=')
        if current[0] == 'flow':
            
            # current[1] = 'BeamLoader, AerogridLoader, StaticCoupled, BeamLoads, AerogridPlot' + '\n'
            current[1] = 'BeamLoader, AerogridLoader, StaticTrim, BeamLoads, AerogridPlot, BeamPlot' + '\n'
            nu_lines.append(' = '.join(current))
            del write_lines[-1]
            write_lines.append(nu_lines[-1])
            continue

        if current[0] == 'gravity':
            current[1] = str(trim_set['g_loading']*9.81) + '\n'
            nu_lines.append(' = '.join(current))
            del write_lines[-1]
            write_lines.append(nu_lines[-1])
            continue

        if current[0] == 'initial_alpha': 
            if logs['alpha_log'] != []:
                current[1] = str(logs['alpha_log'][-1]*np.pi/180) + '\n'
                nu_lines.append(' = '.join(current))
                del write_lines[-1]
                write_lines.append(nu_lines[-1])
            else:
                current[1] = str(trim_set['alpha']) + '\n'
                nu_lines.append(' = '.join(current))    
                del write_lines[-1]
                write_lines.append(nu_lines[-1])
            continue

        if current[0] == 'initial_deflection':
            if logs['elevator_log'] != []:
                current[1] = str(logs['elevator_log'][-1]*np.pi/180) + '\n'
                nu_lines.append(' = '.join(current))
                del write_lines[-1]
                write_lines.append(nu_lines[-1])
            else:
                current[1] = str(trim_set['cs_deflection']) + '\n'
                nu_lines.append(' = '.join(current))   
                del write_lines[-1]
                write_lines.append(nu_lines[-1])
            continue

        if current[0] == 'initial_thrust':
            if logs['thrust_log'] != []:
                current[1] = str(logs['thrust_log'][-1]) + '\n'
                nu_lines.append(' = '.join(current))
                del write_lines[-1]
                write_lines.append(nu_lines[-1])
            else:
                current[1] = str(trim_set['initial_thrust']) + '\n'
                nu_lines.append(' = '.join(current))   
                del write_lines[-1]
                write_lines.append(nu_lines[-1])
            continue
        
        nu_lines.append(i) 

    os.remove(route + '/' + case_name + '.sharpy')
    
    with open(route + '/' + case_name + '.sharpy', 'w') as f:
        for line in write_lines:
            f.write(line)
            
    return trim_set['g_loading'] # specifies g_loading so that mass can be calculated accurately