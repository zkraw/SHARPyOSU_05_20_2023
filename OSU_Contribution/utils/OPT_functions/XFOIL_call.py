import subprocess
import numpy as np

def XFOIL_call_constant_chord(XFOIL_PATH, xfoil_set):
    xfoil = subprocess.Popen(XFOIL_PATH, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True)
    
    actions = ["NACA 0015 \n ",
            #["LOAD {0}.dat \n ".format(airfoil_name),
                                    # "{0} \n ".format(airfoil_name),
                                    "PANE \n ",
                                    "OPER \n ",
                                    "Visc {0} \n ".format(xfoil_set['Re']),
                                    "Mach {0} \n ".format(xfoil_set['Mach']),
                                    "PACC \n ",
                                    "polar_file.txt \n\n ",
                                    "ITER {0} \n ".format(xfoil_set['n_iter']),
                                    "aseq {0} {1} {2} \n ". format(xfoil_set['alfa_1'], xfoil_set['alfa_2'], xfoil_set['alfa_step']),
                                    "quit"]
    command = ''.join(actions)
        
    xfoil.communicate(input=command)
    