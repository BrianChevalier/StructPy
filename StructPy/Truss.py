from StructPy import structural_classes as sc
from StructPy import cross_sections as xs
from StructPy import materials as ma
import numpy as np
import logging

def flatten(items):
	return sum(items, [])

class TrussMember(sc.Member):
	"""
	Define Truss Member class. Only allows axial loading.
	"""
	
	@property
	def k(self):
		"""Local member stiffness matrix"""
		# direct stiffness method

		L = self.length
		E = self.material.E
		A = self.cross.A
		a = (A*E)/L

		return a * np.array([
			[1, -1],
			[-1, 1]
		])
	
	@property
	def T(self):
		"""Local to global transformation matrix"""
		l = self.unVec[0]
		m = self.unVec[1]
		return np.array([
			[l, m, 0, 0],
			[0, 0, l, m]
		])

	@property
	def DoF(self):
		SN = self.SN.n  # this is the node number
		EN = self.EN.n
		
		return [
			SN*2, SN*2 + 1,
			EN*2, EN*2 + 1
		]


class Truss(sc.Structure):	
	"""This class builds on the structure class adding truss solving 
	methods"""
	
	#Number of degrees of freedom for this structure type in the global coordinate system
	nDoFPerNode = 2
	
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
		
		member = TrussMember(SN, EN, material, cross, expectedaxial=expectedaxial)
			
		self.members.append(member)
		self.nMembers += 1

		return member

	@property
	def BC(self):
		"""Define Boundary Condition array"""
		return np.array( flatten( [ [node.xfix, node.yfix] for node in self.nodes] ) )
	
	def directStiffness(self, loading):
		"""This executes the direct stiffness method 
		for the structure given some loading"""
		
		eigs, vecs = np.linalg.eig(self.reducedK)
		if np.isclose(eigs, 0).any() == True:
			raise ValueError('Structure is unstable.')
			logging.warning(eigs)
		
		reducedF = loading[self.unrestrainedDoF]
		reducedD = np.linalg.solve(self.reducedK, reducedF)
		
		globalD = self.BC
		globalD[self.unrestrainedDoF] = reducedD
				
		for i, node in enumerate(self.nodes):
			node.xdef = globalD[2*node.n]
			node.ydef = globalD[2*node.n + 1]
		
		for i, member in enumerate(self.members):
			
			l = member.unVec[0]
			m = member.unVec[1]
			
			ind = member.DoF
			A = member.cross.A
			E = member.material.E
			L = member.length
			
			member.axial = (A*E)/L * np.array([l, m, -l, -m]) @ np.array(globalD[ind]).T
			
		return globalD
		