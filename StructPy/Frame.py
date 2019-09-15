import StructPy.structural_classes as sc
import StructPy.cross_sections as xs
import StructPy.materials as ma
import numpy as np

class FrameMember(sc.Member):
	
	@property
	def k(self):
		"""
		Local member stiffness matrix for frame
		"""
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
		"""
		Global transformation matrix for 2-D frame
		"""
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
	def kglobal(self):
		"""
		Define the global stiffness matrix for a frame element
		"""
		return self.T.T @ self.k @ self.T

	@property
	def DoF(self):
		sn = self.SN.n  # this is the node number
		en = self.EN.n
		
		return [
			3*sn, 3*sn+1, 3*sn+2,
			3*en, 3*en+1, 3*en+2
		]
	

class Frame(sc.Structure):
	
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
		BC = []
		for node in self.nodes:
			BC.extend([node.xfix, node.yfix, node.theta])		
		return np.array(BC)
	
	def directStiffness(self, loading):
		index = self.BC == 1
		reducedF = loading[index]
		
		D = np.linalg.solve(self.reducedK, reducedF)
		d = self.BC.astype('float64') #make sure its not an int.
		d[index] = D
				
		for i, node in enumerate(self.nodes):
			node.xdef = d[3*node.n]
			node.ydef = d[3*node.n + 1]
			node.thetadef = d[3*node.n + 2]
			
		return d
