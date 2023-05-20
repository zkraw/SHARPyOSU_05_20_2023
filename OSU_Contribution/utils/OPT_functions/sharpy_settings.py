import os
import numpy as np
import sharpy.utils.algebra as algebra
import h5py as h5

def one_g_trim(route, case_name, logs, trim_set):
    file = open(route + '/' + case_name + '.sharpy')
    lines = file.readlines()
    nu_lines = []
    for i in lines:
        current = i.replace(' ', '').split('=')
        if current[0] == 'flow':
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
    file = open(route + '/' + case_name + '.sharpy')
    lines = file.readlines()
    nu_lines = []
    for i in lines:
        current = i.replace(' ', '').split('=')
        if current[0] == 'flow':
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

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
            
def turbulent_gust(route, case_name, logs, trim_set, fem):
    file = open(route + '/' + case_name + '.sharpy')
    lines = file.readlines()
    nu_lines = []
    for i in lines:
        current = i.replace(' ', '').split('=')
        if current[0] == 'flow':
            current[1] = 'BeamLoader, AerogridLoader, StaticCoupled, BeamLoads, AerogridPlot, BeamPlot, DynamicCoupled' + '\n'
            nu_lines.append(' = '.join(current))
            continue

        if current[0] == 'gravity':
            current[1] = str(1*9.81) + '\n'
            nu_lines.append(' = '.join(current))
            continue

        if current[0] == 'orientation':
            current[1] = str(algebra.euler2quat(np.array([trim_set['roll'], float(logs['alpha_log'][-1]*np.pi/180), trim_set['beta']]))) + '\n'
            nu_lines.append(' = '.join(current))
            continue
        
        nu_lines.append(i) 

    os.remove(route + '/' + case_name + '.sharpy')
    
    with open(route + '/' + case_name + '.sharpy', 'w') as f:
        for line in nu_lines:
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
    aero['control_surface_deflection'][trim_set['tail_cs_index']] = float(logs['elevator_log'][-1]*np.pi/180)

     # creates new fem file with updated thrust values 
    for iii in trim_set['nodes_app']:
        sign = np.sign([fem['app_forces'][iii][1]]) 
        fem['app_forces'][iii][1] = float(sign[0]*(logs['thrust_log'][-1])) 
    
    return aero, fem

 
        
