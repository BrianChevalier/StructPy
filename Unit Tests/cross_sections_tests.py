import unittest
import StructPy.cross_sections as xs


# class name doesn't matter
# must inheret from unittest.TestCase
class TestCross(unittest.TestCase):
	
	# method name must begin with 'test_'
	def test_Rectangle(self):
		xs1 = xs.Rectangle(2, 2)
		self.assertEqual(xs1.A, 2*2)
		
		# test what happens when you update
		# derived properties should dynamically be recalculated
		xs1.h = 3
		xs1.b = 3
		
		self.assertEqual(xs1.A, 3*3)
		
		
	def test_IBeam(self):
		
		b = 1
		d = 1
		tw = 0.1
		tf = 0.1
		
		Ix = (1 / 12 * b * (d + 2 * tf)**3 -1 / 12 * (b - tw) * d**3)
		A = 2 * (b * tf) + d * tw
		
		xs1 = xs.IBeam(b, d, tw, tf)
		
		self.assertEqual(xs1.Ix, Ix)
		self.assertEqual(xs1.A, A)
		
	def test_AISC(self):
		xs1 = xs.AISC('W21X44')
		A = 13
		
		self.assertEqual(xs1.A, A)
		
		xs1 = xs.AISC('C15X33.9')
		A = 8.81
		
		xs1.A = 3
		
		## TO DO: Add printing unit test.		
		
if __name__ == '__main__':
	unittest.main()
