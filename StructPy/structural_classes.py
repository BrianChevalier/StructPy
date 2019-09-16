from StructPy import cross_sections as xs
from StructPy import materials as ma
import matplotlib.pyplot as plt
import numpy as np

try:
	import yaml
except ImportError:
	raise ImportError('The library pyyaml is required to read a yaml file.')

def flatten(items):
	return sum(items, [])	

class Node(object):

	def __init__(self, x, y, n=None, cost=0, fixity='free'):
		self.x = x
		self.y = y
		self.cost = cost
		self.n = n
		self.fixity = fixity
	
	@property
	def BC(self):
		try:
			return self.__class__.fixities[self.fixity]
		except KeyError:
			validKeys = ', '.join([str(key) for key in self.__class__.fixities.keys()])
			raise ValueError(f'Not a valid nodal support type. Valid types: {validKeys}')
	
	def __str__(self):
		return f"{self.__class__.__name__}({self.x:1.1f}, {self.y:1.1f})"

class Member(object):
	"""Define Member base class"""
	
	def __init__(self, SN, EN, material=ma.Steel(), cross=xs.generalSection(), expectedaxial=None):
		self.cross = cross
		self.material = material

		# assign start and end node properties
		self.SN = SN
		self.EN = EN

		# the axial force
		#self.axial = 0
		self.expectedaxial = expectedaxial

	@property
	def vector(self):
		e1 = self.EN.x - self.SN.x
		e2 = self.EN.y - self.SN.y
		return np.array([e1, e2])

	@property
	def length(self):
		return np.linalg.norm(self.vector)

	@property
	def unVec(self):
		return self.vector/self.length

	@property
	def kglobal(self):
		"""
		Define the global stiffness matrix for a frame element
		"""
		return self.T.T @ self.k @ self.T
	
	@property
	def DoF(self):
		"""The global degree of freedom numbering for the start and end nodes"""
		nDoF = self.__class__.nDoFPerNode
		return [nDoF*self.SN.n + i for i in range(nDoF)] + [nDoF*self.EN.n + i for i in range(nDoF)]
	
	def __repr__(self):
		return f'{self.__class__.__name__}({self.SN}, {self.EN})'


class Structure(object):
	"""
	Abstract base class for Truss and Frame classes.
	"""
	def __init__(self, cross=None, material=None):
		self.members = []
		self.nodes = []
		self.nNodes = 0
		self.nMembers = 0

		if cross is None or material is None:
			raise ValueError('Please define default cross section or material type.')
		else:
			self.defaultcross = cross
			self.defaultmaterial = material

	@classmethod
	def from_yaml_file(cls, filePath):

		## Read in yaml file as a dictionary
		with open(filePath, 'r') as stream:
			try:
				data = yaml.safe_load(stream)
			except yaml.YAMLError as exc:
				print(exc)

		## Initialize structure
		xsection = xs.generalSection(**data['XSection'])
		material = ma.Custom(**data['Material'])
		s1 = cls(cross=xsection, material=material)

		## Add nodes
		nodeMap = {} # map nodes to indicies
		for i, node in enumerate(data['Nodes']):
			for key, value in node.items():
				try:
					s1.addNode(value['x'], value['y'], fixity=value['fixity'])
				except KeyError:
					s1.addNode(value['x'], value['y'])

				nodeMap.update({key: i})

		## Add members
		for member in data['Members']:
			for key, value in member.items():
				SN, EN = key.split(',')
				s1.addMember(nodeMap[SN], nodeMap[EN], expectedaxial=value['axial'])

		return s1

	def addNode(self, x, y, cost=0, fixity='free'):
		"""
		Add node to the structure
		"""
		n = self.nNodes  # node number
		node = self.__class__.NodeType(x, y, n=n, cost=cost, fixity=fixity)
		self.nodes.append(node)
		self.nNodes += 1
	
	def addMember(self, SN, EN, material=None, cross=None, expectedaxial=None):
		"""Add member to the structure"""
		SN = self.nodes[SN]
		EN = self.nodes[EN]

		if material is None:
			material=self.defaultmaterial
		if cross is None:
			cross = self.defaultcross
		
		member = self.__class__.MemberType(SN, EN, material, cross, expectedaxial=expectedaxial)
		
		self.members.append(member)
		self.nMembers += 1
	
	@property
	def BC(self):
		"""Define global boundary condition array"""
		return np.array( flatten( [list(node.BC) for node in self.nodes] ))
	
	@property
	def freeDoF(self):
		return self.BC == 1
	
	@property
	def reducedK(self):
		return self.K[np.ix_(self.freeDoF, self.freeDoF)]
	
	@property
	def nDoF(self):
		return self.__class__.nDoFPerNode * self.nNodes
	
	@property
	def K(self):
		"""Build global structure stiffness matrix"""
		
		K = np.zeros([self.nDoF, self.nDoF])
		
		for member in self.members:
			K[np.ix_(member.DoF, member.DoF)] += member.kglobal
		
		return K
	
	def isStable(self):
		eigs, vecs = np.linalg.eig(self.reducedK)
		if np.isclose(eigs, 0).any() == True:
			logging.warning(eigs)
			raise ValueError('Structure is unstable.')
	
	def solve(self, loading):
		"""Execute direct stiffness solving"""
		reducedF = loading[self.freeDoF]
		reducedD = np.linalg.solve(self.reducedK, reducedF)
		
		globalD = self.BC
		globalD[self.freeDoF] = reducedD
		
		return globalD
		
	def directStiffness(self, loading):
		"""This executes the direct stiffness method"""
		
		self.isStable()
		globalD = self.solve(loading)
		nDoFPerNode = self.__class__.nDoFPerNode
		for i, node in enumerate(self.nodes):
			node.deformation = [globalD[nDoFPerNode*node.n+i] for i in range(nDoFPerNode)]
			
		return globalD
		
	def plot(self, show=True, labels=False):
		"""
		Plot the undeformed structure
		"""
		plt.figure(1)
		plt.clf()
		plt.grid(True)

		length = 0
		for member in self.members:
			if member.length > length:
				length = member.length

		for i, node in enumerate(self.nodes):
			plt.scatter([node.x], [node.y], color='#000000', s=100)
			R = 0.2*length  # length of support
			
			if node.BC[0] == 0:
				x1 = node.x - R
				x2 = node.x
				y1 = node.y
				y2 = node.y
				plt.plot([x1, x2], [y1, y2], color='#57d261', lw=10, zorder=-1)

			if node.BC[1] == 0:
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

		if labels == True:
			for i, node in enumerate(self.nodes):
				plt.Circle((node.x, node.y), 100, color='#ffffff')
				plt.annotate(str(i), (node.x, node.y))

		plt.axis('equal');
		if show == True:
			plt.show()

	def plotDeformation(self, scale=100, nfig=1):

		plt.figure(nfig)
		self.plot(show=False)

		for i, node in enumerate(self.nodes):
			plt.scatter([node.x + scale*node.xdef], [node.y + scale*node.ydef], color='#d80000', s=100)

		for i, member in enumerate(self.members):
			x1 = member.SN.x + scale*member.SN.xdef
			y1 = member.SN.y + scale*member.SN.ydef
			x2 = member.EN.x + scale*member.EN.xdef
			y2 = member.EN.y + scale*member.EN.ydef
			plt.plot([x1, x2], [y1, y2], '--', color='#b83939', lw=2, zorder=-1)
		plt.title(f"Truss Deformation Plot (scale: {scale}X)")
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
