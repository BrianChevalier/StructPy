import StructPy.structural_classes as sc
import StructPy.cross_sections as xs
import StructPy.materials as ma
import numpy as np

def flatten(items):
	return sum(items, [])

class FrameMember(sc.Member):
	
	@property
	def k(self):
		"""Local member stiffness matrix"""
		L = self.length
		E = self.material.E
		A = self.cross.A
		I = self.cross.Ix

		a = (A*E)/L
		b = (E*I)/L
		c = (E*I)/L**2
		d = (E*I)/L**3
		
		return np.array([
			[ a,  0,     0,    -a, 0,     0    ],
			[ 0,  12*d,  6*c,  0,  -12*d, 6*c  ],
			[ 0,  6*c,   4*b,  0,  -6*c,  2*b  ],
			[ -a, 0,     0,    a,  0,     0    ],
			[ 0,  -12*d, -6*c, 0,  12*d,  -6*c ],
			[ 0,  6*c,   2*b,  0,  -6*c,   4*b ]
		])

	@property
	def T(self):
		"""Local to global transformation matrix"""
		l = self.unVec[0]
		m = self.unVec[1]

		return np.array([
			[l,  m, 0, 0,  0, 0],
			[-m, l, 0, 0,  0, 0],
			[0,  0, 1, 0,  0, 0],
			[0,  0, 0, l,  m, 0],
			[0,  0, 0, -m, l, 0],
			[0,  0, 0, 0,  0, 1]
		])

	@property
	def DoF(self):
		SN = self.SN.n  # this is the node number
		EN = self.EN.n
		
		return [
			3*SN, 3*SN + 1, 3*SN + 2,
			3*EN, 3*EN + 1, 3*EN + 2
		]
	

class Frame(sc.Structure):
	"""Frame class"""
	nDoFPerNode = 3
	
	def addMember(self, SN, EN, material=None, cross=None, expectedaxial=None):
		"""
		Add member to the structure
		"""
		SN = self.nodes[SN]
		EN = self.nodes[EN]

		if material is None:
			material=self.defaultmaterial
		if cross is None:
			cross = self.defaultcross
		
		member = FrameMember(SN, EN, material, cross, expectedaxial=expectedaxial)
			
		self.members.append(member)
		self.nMembers += 1

		return member
		
	@property
	def BC(self):
		"""Define Boundary Condition array"""
		return np.array( flatten( [ [node.xfix, node.yfix, node.theta] for node in self.nodes] ) )
	
	def directStiffness(self, loading):
		
		reducedF = loading[self.unrestrainedDoF]
		reducedD = np.linalg.solve(self.reducedK, reducedF)
		
		globalD = self.BC
		globalD[self.unrestrainedDoF] = reducedD  #substitute actual deformations
		
		for i, node in enumerate(self.nodes):
			node.xdef = globalD[3*node.n]
			node.ydef = globalD[3*node.n + 1]
			node.thetadef = globalD[3*node.n + 2]
			
		return globalD
