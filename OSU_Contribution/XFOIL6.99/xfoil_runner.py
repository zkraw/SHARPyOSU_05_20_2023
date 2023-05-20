from cmath import pi
import os
from re import M
import subprocess
import numpy as np

# %% Inputs
c_root = 2      #Chord Length at root (m)
c_tip = 1      #Chord Length at tip (m)
rho = 1.225     #density (kg/m^3)
u = 1.802e-05   #dynamic viscosity (kg/(m*s))
Vel = 10         #Velocity (m/s)
a = 343         #speed of sound (m/s)
Mach = Vel / a  #Mach number
m = 5           #number of strips in wing
b = 6           #span of wing
width = b/m     #width of strip
chord = .5
Drag = np.zeros(m)


airfoil_name = "naca0015"
Re = rho * Vel * chord / u
n_iter = 100

alfa_1 = -16
alfa_2 = 16
alfa_step = .5
# %% XFOIL input file writer 
if os.path.exists("polar_file.txt"):
    os.remove("polar_file.txt")

input_file = open("input_file.in", 'w')
input_file.write("LOAD {0}.dat\n".format(airfoil_name))
input_file.write(airfoil_name + '\n')
input_file.write("PANE\n")
input_file.write("OPER\n")
input_file.write("Visc {0}\n".format(Re))
input_file.write("Mach {0}\n".format(Mach))
input_file.write("PACC\n")
input_file.write("polar_file.txt\n\n")
input_file.write("ITER {0}\n".format(n_iter))
input_file.write("aseq {0} {1} {2}\n".format(alfa_1, alfa_2, alfa_step))
input_file.write("\n\n")
input_file.write("quit\n")
input_file.close()

subprocess.call("xfoil.exe < input_file.in", shell=True)

polar_data = np.loadtxt("polar_file.txt", skiprows=12)

print(polar_data)

# Cd0 = polar_data[2]
# print(Cd0)

# dS = chord[i] * width  #

# Drag[i] = Cd0 * rho * Vel**2 * dS / 2

# print(Drag)
# Profile_drag = sum(Drag)
# print(Profile_drag)