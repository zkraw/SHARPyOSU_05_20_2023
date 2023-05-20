import sympy as sp

b = .03
h = .0185
d = .01056
tt = .00012
e = 20.7e9
v =  .19
g = e / (2*(1+v))
density = 1335

gj_og = 593
eiy_og = 112

test = (b*h**3)/ 12 - ((b-2*tt)*(h-2*tt)**3)/12
max = (b*h**3)/ 1

area = b*h - (b-2*tt)*(h-2*tt)
mass = density*area*6

m_bar = area * 1335 # area times density gives mass per span

EIyy = test * e
eiy_everything = eiy_og - EIyy

EIzz = (((h*b**3)/ 12 - ((h-2*tt)*(b-2*tt)**3)/12) + (area + d**2))* e
I_yy = test

target = .95* 593
gj =  (4 * (b*h) ** 2) / ((2*b + 2*h)/ (g*tt)) # GJ from the spar only
gj_everything = gj_og - gj

from sympy.abc import t
Eqn = sp.Eq((b*h**3)/ 12 - ((b-2*t)*(h-2*t)**3)/12, I_yy)
ans = sp.solve(Eqn)    
print(ans)

from sympy.abc import v
Eqn2 = sp.Eq(e / (2*(1 + v)), g)
ans = sp.solve(Eqn2)
print(ans)