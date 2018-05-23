

class NodalLoad():
	
	def __init__(self, node, xforce, yforce):
		self.node = node
		self.xforce = xforce
		self.yforce = yforce
	
	
	
class ElementLoad():
	"""Distributed loading"""
	
	def __init__(self, element, type, w):
		self.element = element
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
			
