class NodalLoad():
	
	def __init__(self, n, x=0, y=0, isGlobal=True):
		self.n = n
		self.x = x
		self.y = y
	
class ElementLoad():
	"""Distributed loading"""
	
	def __init__(self, m, type, w):
		self.m = m
		self.w = w
		self.type = type
		
	@property
	def equivalentNodalLoad(self):
		
		L = self.element.length
		w = self.w
		
		if self.type == 'constant':
			return w * L / 2
		elif self.type == 'point':
			pass


class Loading(object):
	
	def __init__(self):
		self.nodalLoads = []
		self.elementLoading = []
	
	def addNodalLoad(self, n, x=0, y=0, isGlobal=True):
		self.nodalLoads.append([
			NodalLoad(n, x=x, y=y, isGlobal=isGlobal)
		])
	
	def addElementLoad(self, m, x=0, y=0, isGlobal=True):
		self.elementLoading.append([
			ElementLoad(m, x=x, y=y, isGlobal=isGlobal)
		])
