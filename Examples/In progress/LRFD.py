import matplotlib.pyplot as plt
import numpy as np


def factoredLoads(DL=0, LL=0, S=0, R=0, W=0, E=0, Lr=0):
	#LRFD
	LC1 = 1.4*DL
	LC2 = 1.2*DL + 1.6*LL +0.5*max(Lr,S,R)
	LC3 = 1.2*DL + 1.6*max(Lr,S,R) + max(0.5*LL,0.8*W)
	LC4 = 1.2*DL + 1.6*W + 0.5*LL + 0.5*max(Lr,S,R)
	LC5 = 1.2*DL + 1.0*E + 0.5*LL + 0.2*S
	LC6 = 0.9*DL + 1.6*W
	LC7 = 0.9*DL + 1.0*E
	
	LC = [LC1, LC2, LC3, LC4, LC5, LC6, LC7]
	
	for i, item in enumerate(LC, start=1):
		print('LC%i: %.2f' % (i, item))
	print('\n')
	
	return max(LC)

Q = factoredLoads(DL=16, LL=32)
w = factoredLoads(DL=1, LL=3)
P = factoredLoads(DL=100, LL=200)

L1 = 12
L2 = 20
Ay = 140

L = [0, L1, L1+L2]
M = [lambda x: Ay * x - 0.5 * w * x**2, 
			lambda x: Ay * x - Q*(x - L1) - 0.5 * w * x ** 2]

M1 = []
x1 = []
for i, func in enumerate(M):
	x = np.linspace(L[i], L[i+1])
	x1.extend(x)
	M1.extend([func(x) for x in x])


plt.grid(True)
plt.plot(x1, M1, lw=3)
plt.show()

		
def Cb(M, M1, L0, L1):
	# function requires list of moment expressions, list of all moments, 
	# and the starting and ending length to return the Cb for each span.
	Mmax = max(M1)
	Length = (L1-L0)
	Ma = M(L0 + Length/4)
	Mb = M(L0 + Length/2)
	Mc = M(L0 + 3*Length/4)
	
	print('Mmax: %.3f' % Mmax)
	print('Ma: %.3f' % Ma)
	print('Mb: %.3f' % Mb)
	print('Mc: %.3f' % Mc)
	
	top = 12.5 * Mmax
	bot = 2.5 * Mmax + 3*Ma + 4*Mb + 3*Mc
	
	return top / bot
	
for i, func in enumerate(M):
	print('Cb: %.3f' % Cb(func, M1, L[i], L[i+1]))
	
	
	
