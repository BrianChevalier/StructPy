from StructPy import structural_classes as sc
from StructPy import cross_sections as xs
from StructPy import materials as ma
import numpy as np
import logging


class TrussNode(sc.Node):

	DoFNames = ['x', 'y']
	fixities = {   # x, y
		'free'   :  (1.0, 1.0),
		'pin'    :  (0.0, 0.0),
		'roller' :  (1.0, 0.0),
		'yroller':  (0.0, 1.0),
	}
	
	@property
	def deformation_dict(self):
		return {
			'x':  self.deformation[0],
			'y':  self.deformation[1]
		}
	
class TrussMember(sc.Member):
	"""
	Define Truss Member class. Only allows axial loading.
	"""
	
	nDoFPerNode = 2
	
	@property
	def k(self):
		"""Local member stiffness matrix"""

		a = (self.cross.A * self.material.E)/self.length
		
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
	def axial(self):
		l = self.unVec[0]
		m = self.unVec[1]
		
		A = self.cross.A
		E = self.material.E
		L = self.length
		deformation = self.SN.deformation + self.EN.deformation
		return (A*E)/L * np.array([l, m, -l, -m]) @ np.array(deformation).T

class Truss(sc.Structure, sc.Planar):	
	"""This class builds on the structure class adding truss methods"""
	
	#Number of degrees of freedom for this structure type in the global coordinate system
	whatDoF = ['x', 'y']
	nDoFPerNode = len(whatDoF)
	MemberType = TrussMember
	NodeType = TrussNode
	