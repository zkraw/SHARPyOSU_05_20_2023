
import sympy as sp
import numpy as np

eiy = 532000

# percentage of eiy
pereiy = .9

# percentage of eiz
pereiz = .5

# percetage of gj
pergj = .9

# percentage of the mass
permass = .55

# maximum area using geometry optimizer script 
b =  0.11999999999999977
h =  0.13235760000000032

# Aluminum 6061 T4
# density = 2720 
# e = 68.9e9
# g = 26.2e9

# 70% intermediate carbon density and 30% neat resin 
density = 1790 * .7 + 1200 * .3
e = 155e9
pois = .309
g = e / (2*(1+ pois))

# ---------------------END OF INPUTS------------------------------

eiy_everything = 6e3
gj_everything = 3e3
mass_everything = .375

print('EIy_Everything = ' + str(eiy_everything))

# solving for the thickness given the area
from sympy.abc import t 
Eqn = sp.Eq((b*h**3)/ 12 - ((b-2*t)*(h-2*t)**3)/12, (eiy - eiy_everything)/e) 
ans = sp.solve(Eqn)
tt = ans[0]
print('Thickness = ' + str(tt))    

area = (b*h - (b-2*tt)*(h-2*tt)) 
m_bar_spar = density*area

print("m_bar =" + str(m_bar_spar + mass_everything))
print('m_bar_ev = ' + str(mass_everything))

# calculating the appropriate G
gj_spar =  (4 * (b*h) ** 2) / ((2*b + 2*h)/ (g*tt)) # GJ from the spar only
print('GJ total =' + str(gj_spar + gj_everything))
print('GJ_everything =' + str(gj_everything))

# calculating EIz
Eiz = e*((h*b**3)/ 12 - ((h-2*tt)*(b-2*tt)**3)/12)
print('EIz spar total =' + str(Eiz/pereiz))
print('EIz_everything =' + str(Eiz*(1-pereiz)))

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


# ea = 1.5e7
# gay = 1.5e5
# gaz = 1.5e5
# gj = 1.5e4
# eiy = 3e4
# eiz = 6e7

# # assumed hollow rectangle 
# global a, ay, az, jj, iyy, izz

# a = ea / e
# ay = gay/ g
# az = gaz / g
# jj = gj / g
# iyy = eiy / e
# izz = eiz / e

# from sympy import symbols, Eq, solve

# b, h, t = symbols('b, h, t')

# eq_1 = Eq((b*h - (b-2*t)*(h-2*t)), a) # area 
# # eq_2 = Eq(, ay)
# # eq_3 = Eq(, az)
# eq_4 = Eq((4 * (b*h) ** 2) / ((2*b + 2*h)/ (t)), j)
# eq_5 = Eq((b*h**3)/ 12 - ((b-2*t)*(h-2*t)**3)/12, iy)
# eq_6 = Eq((h*b**3)/ 12 - ((h-2*t)*(b-2*t)**3)/12, iz)

# def equations(vars):
#     b, h, t = vars
#     eq_4 = (4 * (b*h) ** 2) / ((2*b + 2*h)/ (t)) - jj
#     eq_5 = (b*h**3)/ 12 - ((b-2*t)*(h-2*t)**3)/12 - iyy
#     eq_6 = (h*b**3)/ 12 - ((h-2*t)*(b-2*t)**3)/12 - iz
#     return eq_4, eq_5, eq_6

# from scipy.optimize import fsolve

# b, h, t = fsolve(equations, (.3 , .1, .001))

# sol = solve([eq_4, eq_5, eq_6], (b, h, t)) #force=True, manual=True, set=True)