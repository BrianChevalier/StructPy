
class Steel():
	
	def __init__(self, E=29000, Fy=50, Fu=85, cost=0):
		self.E = E
		self.Fy = Fy
		self.Fu = Fu
		self.cost = cost
		
		
class A992(Steel):
	
	def __init__(self):
		self.E = 29000 #ksi
		self.Fy = 50 #ksi
		self.Fu = 65 #ksi
		self.cost = 0
		self.density = 0.2836 #lb/in^3
		

class A36(Steel):
	
	def __init__(self):
		self.E = 29000
		self.Fy = 36 #ksi
		self.Fu = 58 #ksi
		self.poisson = 0.26
		
		
class Custom():
	"""
	>>> ma1 = Custom(E=29000, fy=60)
	"""
	def __init__(self, E=None, fy=None):
		"""
		>>> ma1 = Custom(E=29000, fy=60)
		>>> ma1.E
		29000
		>>> ma1.fy
		60
		"""
		self.E = E
		self.fy = fy
