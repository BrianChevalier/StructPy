import unittest
import StructPy.cross_sections as xs
import StructPy.structural_classes as sc
import StructPy.materials as ma
import StructPy.Truss as Truss
import numpy as np
import math

# class name doesn't matter
# must inheret from unittest.TestCase
class Struct_Classes(unittest.TestCase):
	
	def test_Truss(self):
			
			xs1 = xs.generalSection(A=1.2)
			ma1 = ma.Custom(E=30*10**6) #psi
			
			# Example 6.2.4
			T1 = Truss.Truss(cross=xs1,material=ma1)
			
			# Make nodes
			n1 = T1.addNode(0,-180)
			n2 = T1.addNode(0,0)
			n3 = T1.addNode(120,0,fixity='pin')
			n4 = T1.addNode(120,-180,fixity='wallslider')
			
			m1 = T1.addMember(0,3)
			m2 = T1.addMember(1,2)
			m3 = T1.addMember(1,0)
			m4 = T1.addMember(2,3)
			m5 = T1.addMember(0,2)
			m6 = T1.addMember(1,3)
			
			length = math.sqrt(120**2+180**2)
			
			self.assertEqual(m1.length, 120)
			self.assertEqual(m2.length, 120)
			self.assertEqual(m3.length, 180)
			self.assertEqual(m4.length, 180)
			
			self.assertEqual(m5.length, length)
			self.assertEqual(m6.length, length)
			
			self.assertEqual(m1.unVec[0],1)
			self.assertEqual(m1.unVec[1],0)
			
			loading = np.array([0,0,0,-4000,0,0,0,0])
			T1.directStiffness(loading)
			
			self.assertAlmostEqual(m1.axial[0,0],1333 + 1/3)
			self.assertAlmostEqual(m2.axial[0,0], - 1333 - 1/3)
			self.assertAlmostEqual(m3.axial[0,0], 2000)
			self.assertAlmostEqual(m4.axial[0,0], -2000)
			self.assertAlmostEqual(m5.axial[0,0], -2403.7008503093257)
			self.assertAlmostEqual(m6.axial[0,0], 2403.7008503093257)
			
if __name__ == '__main__':
	unittest.main()
