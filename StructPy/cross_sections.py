import numpy as np
import math
import matplotlib.pyplot as plt
import openpyxl as xl
import StructPy.Resources.pickle_sections as pxs


class Section():
	"""This will define the cross section class that all other cross sections
	will adopt from. All other cross sections that adopt from the `Section`
	class will inherit the classes methods and properties."""

	def __init__(self):
		pass
		
	@property
	def rx(self):
		"""Describes the distribution of cross sectional area about it's centroid.
		The radius if the area was condensed into a circle without changing the
		moment of inertia."""
		return math.sqrt(self.Ix / self.A)
		
	@property
	def ry(self):
		return math.sqrt(self.Iy / self.A)
		
	@property
	def J(self):
		"""Polar moment of inertia"""
		return self.Ix + self.Iy
		
	def __str__(self):
		st0 = '%s\n' % self.name
		st1 = 'A: %.2f\n' % self.A
		st2 = 'Ix: %.2f\n' % self.Ix
		return st0 + st1 + st2
		
	def plot(self):
		"""Plot the figure"""
		plt.figure()
		plt.grid(True)
		plt.axis('equal')
		plt.title(self.name)
		plt.plot(self.xpts, self.ypts, color='#a2a2a2', lw=5, zorder=-1)
		plt.scatter(self.xpts, self.ypts, color='#000000', s=25)
		plt.show()
		
		
class AISC():
	
	def __init__(self, AISCName):		
		#the database is preserved in the pickleditem.txt file for speed
		try:
			data, labels = pxs.unPickleObject()
			self.data = data
			self.labels = labels
		except:
			pxs.main()
			data, labels = pxs.unPickleObject('pickleditem.txt')
			self.data = data
			self.labels = labels
		
		#index is the row of the shape.
		try:
			index = [i for i, s in enumerate(labels) if AISCName in s]
			self.row_number = index[0]
		except IndexError:
			raise ValueError("Shape ID not a valid AISC Shape")
		
		#assign values from the database as properties
		def item(col):
			return self.data[self.row_number][col]
		self.W = item(4)
		self.A = item(5)
		self.Ix = item(38)
		self.Iy = item(42)
		self.rx = item(41)
		self.ry = item(45)
		self.Zx = item(39)
		self.Zy = item(43)
		self.Sx = item(40)
		self.Sy = item(44)
		self.rts = item(74)
		self.ho = item(75)
		self.J = item(49)
		
						
	def printProperties(self):
		"""Print all AISC shape properties in an easy to read table."""
		row_num = self.row_number
		data = self.data
		
		print('Here are the properties')
		print('Section %s in row %i' % (data[row_num][2], row_num))
		
		for col in range(4, 84):
			colname = data[0][col]
			val = data[row_num][col]
			if val == '–':
				pass
			else:
				print('{0:<3} | {1:<8} | {2:>8}'.format(col, colname, str(val)))		
		
		
class customSection(Section):
	
	def __init__(self, xpts, ypts):
		self.xpts = xpts
		self.ypts = ypts
		
		

class generalSection(Section):
	"""This is a general cross section. Define custom properties
	you want to use without defining points"""

	def __init__(self, Ix=1, Q=1, A=1):
		self.Ix = Ix
		self.Q = Q
		self.A = A


class IBeam(Section):
	def __init__(self, b, d, tw, tf):
		# parameters
		self.name = 'I-beam'
		self.b = b
		self.d = d
		self.tw = tw
		self.tf = tf

	# Derived properties of the input parameters
	# If the parameters are updated then so are properties
	@property
	def A(self):
		A = 2 * (self.b * self.tf) + self.d * self.tw
		return A

	@property
	def Ix(self):
		var1 = (1 / 12 * self.b * (self.d + 2 * self.tf)**3
					-1 / 12 * (self.b - self.tw) * self.d**3)
		return var1
		
	@property
	def Iy(self):
		var1 = 2/12 * self.tf**3 * self.b
		var2 = 1/12 * self.tw**3 * self.d
		return var1 + var2

	@property
	def Sx(self):
		"""Section Modulus"""
		var1 = self.b * self.h**2
		var2 = self.h**3 / (self.h + 2 * self.tf)
		var3 = self.b - self.tw
		return 1 / 6 * (var1 - var2 * var3)

	@property
	def Qmax(self):
		return self.Q(0)

	def Q(self, zeta):
		# The expression for Q evaluated at zeta
		var1 = 1 / 8 * self.d**2 * self.tw - 1 / 2 * zeta**2 * self.tw
		var2 = 1 / 2 * self.b * (
			1 / 2 * self.d + self.tf)**2 - 1 / 8 * self.b * self.d**2
		return var1 + var2
		
	@property
	def xpts(self):
		b = self.b
		tw = self.tw
		return [0, b, b, b/2+tw/2, b/2+tw/2, b, b, 0, 0, b/2-tw/2, b/2-tw/2, 0, 0]
		
	@property
	def ypts(self):
		d = self.d
		tf = self.tf
		return [0, 0, tf, tf, tf+d, tf+d, 2*tf+d, 2*tf+d, tf+d, tf+d, tf, tf, 0]

	# define what happens when you call print on the object
	# def __str__(self):
	#	string = 'I-beam Parameters:\nb=%.3f\nd=%.3f\ntf=%.3f\ntw=%.3f'
	#	variables = (self.b, self.d, self.tf, self.tw)
	#	return  string % variables


class Rectangle(Section):
	def __init__(self, b, h):
		self.name = 'Rectangle'
		self.b = b
		self.h = h

	# Derived properties
	@property
	def A(self):
		A = self.b * self.h
		return A

	@property
	def Ix(self):
		Ix = 1 / 12 * self.b * self.h**3
		return Ix

	@property
	def Iy(self):
		return 1 / 12 * self.b**3 * self.h

	@property
	def Sx(self):
		return (self.b * self.h**2) / 6

	def Q(self, zeta):
		Q = 1 / 4 * self.h * (1 / 4 * self.h**2 - zeta**2)
		return Q
		
	@property
	def xpts(self):
		b = self.b
		return [0, b, b, 0, 0]
		
	@property
	def ypts(self):
		h = self.h
		return [0, 0, h, h, 0]


class Circle(Section):
	def __init__(self, r):
		self.name = 'Circle'
		self.r = r

	@property
	def d(self):
		return 2 * self.r

	@property
	def A(self):
		return math.pi * self.r**2

	@property
	def Ix(self):
		return math.pi * self.r**4 / 4

	@property
	def Sx(self):
		return (math.pi * self.d**3) / 32
		
	@property
	def xpts(self):
		r = self.r
		return [r*math.cos(theta) for theta in np.linspace(0, 2*math.pi, 100)]
		
	@property
	def ypts(self):
		r = self.r
		return [r*math.sin(theta) for theta in np.linspace(0, 2*math.pi, 100)]


class HSS(Section):
	def __init__(self, H, h, B, b):
		self.name = 'HSS'
		self.H = H
		self.h = h
		self.B = B
		self.b = b

	@property
	def ybar(self):
		return self.H / 2

	@property
	def xbar(self):
		return self.B / 2

	@property
	def A(self):
		return self.B * self.H - self.b - self.h

	@property
	def Ix(self):
		var1 = self.B * self.H**3 / 12
		var2 = self.b * self.h**3 / 12
		return var1 - var2

	@property
	def Iy(self):
		var1 = self.B**3 * self.H / 12
		var2 = self.b**3 * self.h / 12
		return var1 - var2

	@property
	def Sx(self):
		return self.Ix / self.ybar

	@property
	def Sy(self):
		return self.Iy / self.xbar


class hollowCircle(Section):
	def __init__(self, ro, ri):
		self.name = 'Hollow Circle'
		self.ro = ro
		self.ri = ri

	@property
	def do(self):
		return 2 * ro

	@property
	def di(self):
		return 2 * ri

	@property
	def A(self):
		return math.pi * (self.ro**2 - self.ri**2)

	@property
	def Ix(self):
		return math.pi / 4 * (self.ro**4 - self.ri**4)

	@property
	def Sx(self):
		var1 = math.pi * (self.do**4 - self.di**4)
		var2 = 32 * self.do
		return var1 / var2


class Triangle(Section):
	def __init__(self, b, h):
		self.b = b
		self.h = h

		@property
		def A(self):
			return 1 / 2 * self.b * self.h


if __name__ == "__main__":
	pass

