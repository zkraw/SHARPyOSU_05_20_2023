
import sympy as sp
import numpy as np

gj = 150000
eiy = 30000
eiz = 6000000

# percentage of eiy
pereiy = .9

# percentage of gj
pergj = .9

# percentage of the mass
permass = .5

m_bar = .75
density = 1200 # aluminum 6061-T6 density

# maximum area using geometry optimizer script 
b = 0.09999999999999998 
h = 0.1171 

# ---------------------END OF INPUTS------------------------------

gj_everything = (1-pergj) * gj
print('GJ_Everything = ' + str(gj_everything))

eiy_everything = (1-pereiy) * eiy
print('EIy_Everything = ' + str(eiy_everything))

m_bar_spar = m_bar * permass
print('Mass_Everything = ' + str(m_bar * (1-permass))) 
area = m_bar_spar/density

# solving for the thickness given the area
from sympy.abc import t 
Eqn = sp.Eq(b*h - (b-2*t)*(h-2*t), area) 
ans = sp.solve(Eqn)
tt = ans[0]
print('Thickness = ' + str(tt))    

# calculating the appropriate G
j =  (4 * (b*h) ** 2) / ((2*b + 2*h)/ (tt)) # GJ from the spar only
g = (gj-gj_everything) / j 
print('G = ' + str(g))    

# calculaating the appropriate E
iy = (b*h**3)/ 12 - ((b-2*tt)*(h-2*tt)**3)/12
e = (eiy - eiy_everything) / iy 
print('E = ' + str(e))   

iz = (h*b**3)/ 12 - ((h-2*tt)*(b-2*tt)**3)/12
check = e*iz
# solving for the geometric properties by assuming a material 
eiy_test = np.linspace(30000, 15000, 15)
t_test = []
for i in range(len(eiy_test)):
    spar = eiy_test[i] - eiy_everything
    iy = spar / e
    from sympy.abc import t 
    Eqn = sp.Eq((b*h**3)/ 12 - ((b-2*t)*(h-2*t)**3)/12, iy) 

    ans = sp.solve(Eqn)
    t_test.append(ans[0])

# Aluminum 2014
e = 74.5e9
g = 27.6e9

ea = 1.5e7
gay = 1.5e5
gaz = 1.5e5
gj = 1.5e4
eiy = 3e4
eiz = 6e7

# assumed hollow rectangle 
global a, ay, az, jj, iyy, izz

a = ea / e
ay = gay/ g
az = gaz / g
jj = gj / g
iyy = eiy / e
izz = eiz / e

from sympy import symbols, Eq, solve

b, h, t = symbols('b, h, t')

eq_1 = Eq((b*h - (b-2*t)*(h-2*t)), a) # area 
# eq_2 = Eq(, ay)
# eq_3 = Eq(, az)
eq_4 = Eq((4 * (b*h) ** 2) / ((2*b + 2*h)/ (t)), j)
eq_5 = Eq((b*h**3)/ 12 - ((b-2*t)*(h-2*t)**3)/12, iy)
eq_6 = Eq((h*b**3)/ 12 - ((h-2*t)*(b-2*t)**3)/12, iz)

# def equations(vars):
#     b, h, t = vars
#     eq_4 = (4 * (b*h) ** 2) / ((2*b + 2*h)/ (t)) - jj
#     eq_5 = (b*h**3)/ 12 - ((b-2*t)*(h-2*t)**3)/12 - iyy
#     eq_6 = (h*b**3)/ 12 - ((h-2*t)*(b-2*t)**3)/12 - iz
#     return eq_4, eq_5, eq_6

# from scipy.optimize import fsolve

# b, h, t = fsolve(equations, (.3 , .1, .001))

# sol = solve([eq_4, eq_5, eq_6], (b, h, t)) #force=True, manual=True, set=True)