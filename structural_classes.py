import cross_sections as xs
import materials as ma
import matplotlib.pyplot as plt
import numpy as np
import math


class Node():
	"""Define node class"""
	def __init__(self, x, y, z=0, n=0, cost=0, xfix=1, yfix=1, fixity='free'):
		self.x = x
		self.y = y
		self.z = z
		self.cost = cost
		self.n = n #nodal number
		
		#to be removed.
		self.xfix = xfix
		self.yfix = yfix
		
		# Assign boundary conditions
		# N: normal force
		# V: shear force
		# M: moment
		# u: axial displacement
		# w: transverse displacement
		# theta: rotation (rate of change of w)
		
		if fixity == 'free':
			self.N = 0
			self.V = 0
			self.M = 0
			self.u = 1
			self.theta = 1
			self.w = 1
		elif fixity == 'fixed':
			self.N = 1
			self.V = 1
			self.M = 1
			self.u = 0
			self.theta = 0
			self.w = 0
		elif fixity == 'pin':
			self.N = 1
			self.V = 1
			self.M = 0
			self.u = 0
			self.theta = 1
			self.w = 0
		elif fixity == 'wallslider':
			self.N = 1
			self.V = 0
			self.M = 0
			self.u = 0
			self.theta = 1
			self.w = 1
		elif fixity == 'roller':
			self.N = 0
			self.V = 1
			self.M = 0
			self.u = 1
			self.theta = 1
			self.w = 0
		elif fixity == 'cont':
			self.N = 1
			self.V = 1
			self.M = 1
			self.u = 1
			self.theta = 1
			self.w = 0
		
		#This is the displacement of the node
		#these have to be solved for
		# TO DELETE
		self.xdef = 0
		self.ydef = 0
		self.zdef = 0
		
	def __str__(self):
		string = '(%.1f, %.1f)'
		variables = (self.x, self.y)
		return string % variables
				
						
class Member(Node):
	"""Define Member class"""
	def __init__(self, SN, EN, material=ma.Steel(), cross=xs.generalSection()):		
		self.cross = cross
		self.material = material
		
		# assign start and end node properties
		self.SN = SN
		self.EN = EN
		
		#the axial force 
		self.axial = 0
		
	@property
	def vector(self):
		e1 = self.EN.x - self.SN.x
		e2 = self.EN.y - self.SN.y
		e3 = self.EN.z - self.SN.z
		return np.array([e1, e2, e3])
		
	@property
	def length(self):
		return np.linalg.norm(self.vector)
		
	@property
	def unVec(self):
		return self.vector/self.length
		
	@property
	def axialstiff(self):
		"""This defines the member axial stiffness for the direct stiffness method"""
		return (self.cross.A * self.material.E)/self.length
		
	@property
	def bendingstiff(self):
		return (self.material.E*self.cross.I)/self.length**3
		
	@property
	def torsionstiff(self):
		return (self.G*self.J)/self.length
		
	@property
	def k(self):
		"""Global member stiffness matrix for truss"""
		# direct stiffness method
		l = self.unVec[0]
		m = self.unVec[1]
		n = self.unVec[2]
		
		T = np.matrix([[l, m, 0, 0], [0, 0, l, m]])
		k = np.matrix([[1, -1], [-1, 1]])
		k = self.axialstiff * (T.T * k * T)
		return k
		
	@property
	def kframe(self):
		"""Local member stiffness matrix for frame"""
		L = self.length
		E = self.material.E
		A = self.cross.A
		I = self.cross.Ix
		
		a = (A*E)/L
		b = (E*I)/L
		c = (E*I)/L**2
		d = (E*I)/L**3
		
		v1 = [a,  0,     0,    -a, 0,     0   ]
		v2 = [0,  12*d,  6*c,  0,  -12*d, 6*c ]
		v3 = [0,  6*c,   4*b,  0,  -6*c,  2*b ]
		v4 = [-a, 0,     0,    a,  0,     0   ]
		v5 = [0,  -12*d, -6*c, 0,  12*d,  -6*c]
		v6 = [0,  6*c,   2*b,  0,  -6*c,   4*b]
		
		return np.matrix([v1, v2, v3, v4, v5, v6])
		
	@property
	def T(self):
		"""Global transformation matrix for 2-D frame"""
		l = self.unVec[0]
		m = self.unVec[1]
		n = self.unVec[2]
		
		v1 = [l,  m, 0, 0,  0, 0]
		v2 = [-m, l, 0, 0,  0, 0]
		v3 = [0,  0, 1, 0,  0, 0]
		v4 = [0,  0, 0, l,  m, 0]
		v5 = [0,  0, 0, -m, l, 0]
		v6 = [0,  0, 0, 0,  0, 1]
		
		return np.matrix([v1, v2, v3, v4, v5, v6])
		
	@property
	def trussdof(self):
		SNdof = [self.SN.n*2, self.SN.n*2+1]
		ENdof = [self.EN.n*2, self.EN.n*2+1]
		# This just tells the axial displacement numbering
		# i.e. member 1 has d1, d2, d3, d4 for the start and end nodes
		return SNdof + ENdof
		
	@property
	def frameDoF(self):
		# order u, theta, w
		n = self.SN.n #this is the node number
		sn = [3*n+1, 3*n+2, 3*n+3]
		n = self.EN.n
		en = [3*n+1, 3*n+2, 3*n+3]
		return sn + en
				
		
class Structure(Member, Node):
	"""Define structure class"""
	def __init__(self, cross, material):
		self.members = []
		self.nodes = []
		self.nNodes = 0
		self.nMembers = 0
		self.nDoF = 2
		self.defaultcross = cross
		self.defaultmaterial = material		
				
	def addNode(self, x, y, z=0, cost=0, xfix=1, yfix=1, fixity='free'):
		"""Add node to the structure"""
		n = self.nNodes #node number
		self.nodes.append(Node(x, y, z, n, cost, xfix, yfix, fixity))
		self.nNodes += 1
		
		if z !=0:
			self.nDoF = 3
		
	def addMember(self, SN, EN, material=None, cross=None):
		"""Add member to the structure"""
		SN = self.nodes[SN]
		EN = self.nodes[EN]
		
		if material == None:
			material=self.defaultmaterial
		if cross == None:
			cross = self.defaultcross
		
		self.members.append(Member(SN, EN, material, cross))
		self.nMembers += 1
		
	def plot(self, show=True):
		"""Plot the undeformed structure"""
		plt.figure(1); plt.clf(); plt.grid(True)
		
		for i, node in enumerate(self.nodes):
			plt.scatter([node.x], [node.y], color='#000000', s=100)
			R = 0.2 #length of support
			
			if node.xfix == 0:
				x1 = node.x - R
				x2 = node.x
				y1 = node.y
				y2 = node.y
				plt.plot([x1, x2], [y1, y2], color='#57d261', lw=10, zorder=-1)
				
			if node.yfix == 0:
				x1 = node.x
				x2 = node.x
				y1 = node.y - R
				y2 = node.y
				plt.plot([x1, x2], [y1, y2], color='#57d261', lw=10, zorder=-1)
				
		for i, member in enumerate(self.members):
			x1 = member.SN.x
			y1 = member.SN.y
			x2 = member.EN.x
			y2 = member.EN.y
			plt.plot([x1, x2], [y1, y2], color='#923ab8', lw=3, zorder=-1)
			
		plt.axis('equal'); 
		if show == True:
			plt.show()
		
	def plotDeformation(self, mag=10000, nfig=1):
		
		plt.figure(nfig); self.plot(show=False)
		
		for i, node in enumerate(self.nodes):
			plt.scatter([node.x + mag*node.xdef], [node.y + mag*node.ydef], color='#d80000', s=100)
			
		for i, member in enumerate(self.members):
			x1 = member.SN.x + mag*member.SN.xdef
			y1 = member.SN.y + mag*member.SN.ydef
			x2 = member.EN.x + mag*member.EN.xdef
			y2 = member.EN.y + mag*member.EN.ydef
			plt.plot([x1, x2], [y1, y2], '--', color='#b83939', lw=2, zorder=-1)
		
		plt.show()
		
	def printNodes(self):
		print('\nStructure Nodes:')
		for i, node in enumerate(self.nodes):
			string = 'Node %i: (%.1f, %.1f)'
			variables = (i, node.x, node.y)
			print(string % variables)
			
	def printMembers(self):
		print('\nStructure Members:')
		for i, member in enumerate(self.members):
			string = 'Member %i: (%i --> %i), L = %.1f, f = %.2f'
			variables = (i, member.SN.n, member.EN.n, member.length, member.axial)
			print(string % variables)

