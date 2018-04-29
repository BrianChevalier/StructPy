import numpy as np
import cross_sections as xs

xs1 = xs.IBeam(1, 1, 0.1, 0.1)

L = 10
p = 1
E = 29000

def constant(x, **kwargs):
	return 1

def linearup(s, **kwargs):
	return x

load = constant

def simpsons(f, a, b, n): #function, start, stop, intervals
	if n % 2 == 0:
		h = (b-a)/n
		k = 0.0
		x = a + h
		
		for i in range(1, int(n/2) + 1):
			k += 4*f(x)
			x += 2*h
		
		x = a + 2*h
		for i in range(1, n//2):
			k += 2*f(x)
			x += 2*h
		return (h/3)*(f(a) + f(b) + k)
	else:
		print('n must be even')

I0 = lambda x: p * load(x, L=L)
I1 = lambda x: p * (L-x) * load(x, L=L)
I2 = lambda x: p * (L-x)**2 * load(x, L=L)
I3 = lambda x: p * (L-x)**3 * load(x, L=L)

Int0 = simpsons(I0, 0, L, 100)
Int1 = simpsons(I0, 0, L, 100)
Int2 = simpsons(I0, 0, L, 100)
Int3 = -simpsons(I0, 0, L, 100)

z = np.array([Int0, Int1, Int2, Int3])

a = L
b = L/(E * xs1.Ix)
c = L**2/(2 * E * xs1.Ix)
d = L**3/(6 * E * xs1.Ix)

B = np.matrix([[1, 0, 0, 0, -1, 0, 0, 0], 
							 [a, 1, 0, 0, 0, -1, 0, 0],
							 [c, b, 1, 0, 0, 0, -1, 0],
							 [-d, -c, -a, 1, 0, 0, 0, -1]])
							 
fixed = [1, 1, 0, 0]
free = [0, 0, 1, 1]

BC = np.array(fixed + free)

C = B[:, BC==1]

s = np.linalg.solve(C, z)
