import StructPy.structural_classes as sc
import StructPy.cross_sections as xs
import StructPy.materials as ma
import numpy as np

class Frame(sc.Structure):
	
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
					K[i,j] += member.kframe[index1, index2]
		
		return np.asmatrix(K)
		
	@property
	def BC(self):
		"""Define Boundary Condition array"""
		BC = []
		for node in self.nodes:
			BC.extend([node.u, node.w, node.theta])		
		return np.array(BC)
		
	def directStiffness(self, loading):
		index = self.BC == 1
		reducedK = self.K[index,:][:,index]
		reducedF = loading[index]
		
		D = np.linalg.solve(reducedK, reducedF)
		d = self.BC.astype('float64') #make sure its not an int.
		d[index] = D
				
		for i, node in enumerate(self.nodes):
			node.xdef = d[3*node.n]
			node.ydef = d[3*node.n + 1]
			node.thetadef = d[3*node.n + 2]
		
		for i, member in enumerate(self.members):
			l = member.unVec[0]
			m = member.unVec[1]
			n = member.unVec[2]
			ind = member.trussdof
			member.axial = member.axialstiff * np.matrix([l, m, -l, -m]) * np.asmatrix(d[ind]).T
			
		return d
