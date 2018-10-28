from StructPy import structural_classes as sc
from StructPy import cross_sections as xs
from StructPy import materials as ma
import numpy as np

import decimal

class Truss(sc.Structure):	
	"""This class builds on the structure class adding truss solving 
	methods"""
					
	@property
	def K(self):
		"""Build global structure stiffness matrix for a truss"""
		
		K = np.zeros((2 * self.nNodes, 2 * self.nNodes))
		
		for member in self.members:
			# this generates the correct numbering for nodal dof
			SNdof = [member.SN.n*2, member.SN.n*2+1]
			ENdof = [member.EN.n*2, member.EN.n*2+1]
			
			ds = SNdof + ENdof
			
			for index1, i in enumerate(ds):
				for index2, j in enumerate(ds):
					# index is used for local numbering
					# i,j are used for global numbering
					# this is a local to global transformation
					K[i,j] += member.k[index1, index2]
		
		return np.asmatrix(K)

	@property
	def BC(self):
		"""Define Boundary Condition array"""
		BC = []
		for node in self.nodes:
			x = node.xfix
			y = node.yfix
			BC.append(x)
			BC.append(y)
				
		return np.array(BC)
		
	def directStiffness(self, loading):
		"""This executes the direct stiffness method 
		for the structure given some loading"""
		
		if type(loading) == type([]):
			raise TypeError("input must be np.array")
		
		# rows with displacement
		index = self.BC == 1
		reducedK = self.K[index,:][:,index]
		reducedF = loading[index]
		
		D = np.linalg.solve(reducedK, reducedF)
		d = self.BC.astype('float64') #make sure its not an int.
		d[index] = D.flat
				
		for i, node in enumerate(self.nodes):
			node.xdef = d[2*node.n]
			node.ydef = d[2*node.n + 1]
		
		for i, member in enumerate(self.members):
			
			l = member.unVec[0]
			m = member.unVec[1]
			n = member.unVec[2]
			
			ind = member.trussdof
			A = member.cross.A
			E = member.material.E
			L = member.length
			
			member.axial = (A*E)/L * np.matrix([l, m, -l, -m]) * np.asmatrix(d[ind]).T
			
		return d
		
	@property
	def selfWeightAtNodes(self):
		"""
		This outputs a numpy array of the equivalent nodal loading caused by the
		self weight of the truss in the units of the self weight property
		"""
		nodes = np.zeros(self.nNodes)
		for member in self.members:
			nodes[member.SN.n] += 0.5 * member.length * member.cross.W
			nodes[member.EN.n] += 0.5 * member.length * member.cross.W
			
		return nodes
		
		
