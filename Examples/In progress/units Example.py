import pint
u = pint.UnitRegistry()
Q = u.Quantity

c = 299792458 * u.meter / u.sec

# the following is the same
print(c.to('mph'))
print(c.to('mile/hour'))

rho = 62.4 * u.lb / u.foot**3
print(rho.to('kg/m**3'))


Rg = Q(0.08206, 'L*atm/(mole*K)')
print(Rg.to('J/mol/K'))


import numpy as np

R = 25 * u.feet
L = 5 * u.meter
V = np.pi * R**2 * L
print(V.to('meter**3'))




