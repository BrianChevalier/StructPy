from StructPy import structural_classes as sc
from StructPy import cross_sections as xs
from StructPy import materials as ma
import numpy as np
import logging


class TrussMember(sc.Member):
    """
    Define Truss Member class. Only allows axial loading.
    """
	
    @property
    def kglobal(self):
        """Global member stiffness matrix for truss"""
        # direct stiffness method
        l = self.unVec[0]
        m = self.unVec[1]

        L = self.length
        E = self.material.E
        A = self.cross.A
        a = (A*E)/L

        T = np.array([[l, m, 0, 0], [0, 0, l, m]])
        k = np.array([[1, -1], [-1, 1]])
        k = a * (T.T @ k @ T)
        return k

    @property
    def DoF(self):
        SNdof = [self.SN.n*2, self.SN.n*2+1]
        ENdof = [self.EN.n*2, self.EN.n*2+1]
        return SNdof + ENdof


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
		reducedF = loading[index]
		
		eigs, vecs = np.linalg.eig(self.reducedK)
		if np.isclose(eigs, 0).any() == True:
			raise ValueError('Structure is unstable.')
			logging.warning(eigs)
		
		D = np.linalg.solve(self.reducedK, reducedF)
		d = self.BC.astype('float64') #make sure its not an int.
		logging.info(d)
		d[index] = D
				
		for i, node in enumerate(self.nodes):
			node.xdef = d[2*node.n]
			node.ydef = d[2*node.n + 1]
		
		for i, member in enumerate(self.members):
			
			l = member.unVec[0]
			m = member.unVec[1]
			
			ind = member.DoF
			A = member.cross.A
			E = member.material.E
			L = member.length
			
			member.axial = (A*E)/L * np.array([l, m, -l, -m]) @ np.array(d[ind]).T
			
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
		
		
