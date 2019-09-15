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

		v1 = [a,  0,	 0,	-a, 0,	 0   ]
		v2 = [0,  12*d,  6*c,  0,  -12*d, 6*c ]
		v3 = [0,  6*c,   4*b,  0,  -6*c,  2*b ]
		v4 = [-a, 0,	 0,	a,  0,	 0   ]
		v5 = [0,  -12*d, -6*c, 0,  12*d,  -6*c]
		v6 = [0,  6*c,   2*b,  0,  -6*c,   4*b]

		return np.array([v1, v2, v3, v4, v5, v6])

	@property
	def T(self):
		"""
		Global transformation matrix for 2-D frame
		"""
		l = self.unVec[0]
		m = self.unVec[1]

		v1 = [l,  m, 0, 0,  0, 0]
		v2 = [-m, l, 0, 0,  0, 0]
		v3 = [0,  0, 1, 0,  0, 0]
		v4 = [0,  0, 0, l,  m, 0]
		v5 = [0,  0, 0, -m, l, 0]
		v6 = [0,  0, 0, 0,  0, 1]

		return np.array([v1, v2, v3, v4, v5, v6])

	@property
	def kglobal(self):
		"""
		Define the global stiffness matrix for a frame element
		"""
		return self.T @ self.k @ self.T.T

	@property
	def frameDoF(self):
		n = self.SN.n  # this is the node number
		sn = [3*n, 3*n+1, 3*n+2]
		n = self.EN.n
		en = [3*n, 3*n+1, 3*n+2]
		return sn + en
	

class Frame(sc.Structure):
	
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
	def K(self):
		"""Build global structure stiffness matrix"""
		K = np.zeros((3*self.nNodes, 3*self.nNodes))
		for member in self.members:
			
			SNdof = [member.SN.n*3, member.SN.n*3+1, member.SN.n*3+2]
			ENdof = [member.EN.n*3, member.EN.n*3+1, member.EN.n*3+2]
			
			ds = SNdof + ENdof
			
			for index1, i in enumerate(ds):
				for index2, j in enumerate(ds):
					#index is used for local numbering
					#i,j are used for global numbering
					#this is a local to global transformation
					K[i,j] += member.kglobal[index1, index2]
		
		return np.array(K)
		
	@property
	def BC(self):
		"""Define Boundary Condition array"""
		BC = []
		for node in self.nodes:
			BC.extend([node.xfix, node.yfix, node.theta])		
		return np.array(BC)
	
	@property
	def reducedK(self):
		index = self.BC == 1
		return self.K[index, :][:, index]
	
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
		
# 		for i, member in enumerate(self.members):
# 			l = member.unVec[0]
# 			m = member.unVec[1]
# 			ind = member.frameDoF
# 			member.axial = member.axialstiff * np.array([l, m, -l, -m]) @ np.array(d[ind]).T
			
		return d
