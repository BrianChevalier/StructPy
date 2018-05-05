import structural_classes as sc
import cross_sections as xs
import materials as ma

class Truss(sc.Structure):
	
	@property
	def K(self):
		"""Build global structure stiffness matrix for a truss"""
		
		K = np.zeros((self.gDoF, self.gDoF))
		
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
		#rows with displacement
		index = self.BC == 1
		reducedK = self.K[index,:][:,index]
		reducedF = loading[index]
		
		D = np.linalg.solve(reducedK, reducedF)
		Dflat = np.asarray(D).flatten()
				
		j = 0
		d2 = []
		for i, item in enumerate(self.d):
			if item == 1:
				d2.append(Dflat[j])
				j += 1
			else:
				d2.append(0)
		d2 = np.array(d2)
		
		for i, node in enumerate(self.nodes):
			d = node.d
			node.xdef = d2[d[0]]
			node.ydef = d2[d[1]]
		
		for i, member in enumerate(self.members):
			l = member.unVec[0]
			m = member.unVec[1]
			n = member.unVec[2]
			ind = member.disp
			member.axial = member.axialstiff * np.matrix([l, m, -l, -m]) * np.asmatrix(d2[ind]).T
			
		return d2
	

# Define material
A992 = ma.A992()

# Define cross section
xs1 = xs.generalSection(A=0.1)

# define blank structure
# the cross section and material are defaults for the members
# we will add to this later
s1 = Truss(cross=xs1, material=A992)

# Add nodes to the structure
s1.addNode(0, 0, xfix=0, yfix=0)
s1.addNode(3, 4)
s1.addNode(6, 0, xfix=0, yfix=0)

# Add members to the structure
s1.addMember(0, 1)
s1.addMember(1, 2)
s1.addMember(0, 2)

s1.plot()

Forces = np.matrix('0; 8660; 5000; 0; 0')
displacement = s1.directStiffness(Forces)
s1.plotDeformation(nfig=1, mag=1)

Forces = np.matrix('0; 0; 0; -1000; 0; 0')
s1.directStiffness(Forces)
s1.plotDeformation(nfig=2, mag=1)

s1.printMembers()
