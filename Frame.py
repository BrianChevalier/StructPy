import structural_classes as sc
import cross_sections as xs
import materials as ma
import numpy as np


class Frame(sc.Structure):
	
	@property
	def K(self):
		"""Build global structure stiffness matrix"""
		K = np.zeros((3*self.nNodes, 3*self.nNodes))
		for member in self.members:
			x1 = member.SN
			x2 = member.EN
			ds = x1.frameDoF + x2.frameDoF
			
			for index1, i in enumerate(ds):
				for index2, j in enumerate(ds):
					#index is used for local numbering
					#i,j are used for global numbering
					#this is a local to global transformation
					K[i,j] += member.kframe[index1, index2]
		
		return np.asmatrix(K)
