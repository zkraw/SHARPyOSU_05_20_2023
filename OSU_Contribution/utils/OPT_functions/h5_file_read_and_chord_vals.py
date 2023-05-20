import numpy as np
import h5py as h5
import os

def h5_file_read(route, case):
    # getting aero data
    aero = {}
    file = os.path.join(route, case)
    with h5.File(file + '.aero.h5', 'r') as hdf:
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

    # getting fem data
    fem = {}
    with h5.File(file + '.fem.h5', 'r') as hdf:
        ls = list(hdf.keys())
        for i in ls:
            fem[i] = np.array(hdf.get(i)) 
    
    return aero, fem

def find_chord_vals(aero, beam_opt_elems, beam_opt_sym):
    chord_vals = []
    for i in range(len(beam_opt_elems)):
            if i == 0:
                chord_vals.append(aero['chord'][i][0])
                chord_vals.append(aero['chord'][i][2])
                chord_vals.append(aero['chord'][i][1])
            else:
                chord_vals.append(aero['chord'][i][2])
                chord_vals.append(aero['chord'][i][1])
    
    if beam_opt_sym == True:
        for i in range(len(beam_opt_elems)):
            if i == 0:
                chord_vals.append(aero['chord'][beam_opt_elems[-1] + i + 1][0])
                chord_vals.append(aero['chord'][beam_opt_elems[-1] + i + 1][2])
                chord_vals.append(aero['chord'][beam_opt_elems[-1] + i + 1][1])
            else:
                chord_vals.append(aero['chord'][beam_opt_elems[-1] + i + 1][2])
                chord_vals.append(aero['chord'][beam_opt_elems[-1] + i + 1][1])
    
    return chord_vals

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
def scale_chord_values(aero, beam_opt_elems, beam_opt_sym, scale):
    # finding and updating the chord values for different spans
    chord_vals = []
    for i in range(len(beam_opt_elems)):
            if i == 0:
                chord_vals.append((1/scale)*aero['chord'][i][0])
                aero['chord'][i][0] = chord_vals[-1]
                chord_vals.append((1/scale)*aero['chord'][i][2])
                aero['chord'][i][2] = chord_vals[-1]
                chord_vals.append((1/scale)*aero['chord'][i][1])
                aero['chord'][i][1] = chord_vals[-1]
            else:
                chord_vals.append((1/scale)*aero['chord'][i][2])
                aero['chord'][i][2] = chord_vals[-1]
                chord_vals.append((1/scale)*aero['chord'][i][1])
                aero['chord'][i][1] = chord_vals[-1]
    
    if beam_opt_sym == True:
        for i in range(len(beam_opt_elems)):
            if i == 0:
                chord_vals.append((1/scale)*aero['chord'][beam_opt_elems[-1] + i + 1][0])
                aero['chord'][beam_opt_elems[-1]+ i + 1][0] = chord_vals[-1]
                chord_vals.append((1/scale)*aero['chord'][beam_opt_elems[-1] + i + 1][2])
                aero['chord'][beam_opt_elems[-1] + i + 1][2] = chord_vals[-1]
                chord_vals.append((1/scale)*aero['chord'][beam_opt_elems[-1] + i + 1][1])
                aero['chord'][beam_opt_elems[-1]+ i + 1][1] = chord_vals[-1]
            else:
                chord_vals.append((1/scale)*aero['chord'][beam_opt_elems[-1] + i + 1][2])
                aero['chord'][beam_opt_elems[-1] + i + 1][2] = chord_vals[-1]
                chord_vals.append((1/scale)*aero['chord'][beam_opt_elems[-1] + i + 1][1])
                aero['chord'][beam_opt_elems[-1] + i + 1][1] = chord_vals[-1]
    
    return chord_vals