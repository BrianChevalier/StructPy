import StructPy.structural_classes as sc
import StructPy.cross_sections as xs
import StructPy.materials as ma
import numpy as np


class FrameNode(sc.Node):
	
	fixities = {   # x,   y,   θ
		'free'   :  (1.0, 1.0, 1.0),
		'fixed'  :  (0.0, 0.0, 0.0),
		'pin'    :  (0.0, 0.0, 1.0),
		'roller' :  (1.0, 0.0, 1.0),
		'yroller':  (0.0, 1.0, 1.0),
		'slide'  :  (0.0, 1.0, 0.0)
	}
	
	def __init__(self, x, y, n=None, cost=0, fixity='free'):
		self.x = x
		self.y = y
		self.cost = cost
		self.n = n
		self.fixity = fixity
		
		try:
			self.xfix, self.yfix, self.theta = self.__class__.fixities[fixity]
		except KeyError:
			self.invalidFixityError()

class FrameMember(sc.Member):
	
	nDoFPerNode = 3
	
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
	

class Frame(sc.Structure):
	"""Frame class"""
	
	whatDoF = ['x', 'y', 'θz']
	nDoFPerNode = len(whatDoF)
	MemberType = FrameMember
	NodeType = FrameNode
	
	def directStiffness(self, loading):
		
		self.isStable()
		globalD = self.solve(loading)
		
		for i, node in enumerate(self.nodes):
			node.xdef = globalD[3*node.n]
			node.ydef = globalD[3*node.n + 1]
			node.thetadef = globalD[3*node.n + 2]
			
		return globalD
