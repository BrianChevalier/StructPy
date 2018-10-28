import unittest
import StructPy.cross_sections as xs
import StructPy.structural_classes as sc
import StructPy.materials as ma
import StructPy.Truss as Truss
import numpy as np

# class name doesn't matter
# must inheret from unittest.TestCase
class Struct_Classes(unittest.TestCase):
	
	def setUp(self):
		#Truss
		xs1 = xs.Rectangle(2, 2)
		A992 = ma.A992()
		s1 = Truss.Truss(cross=xs1, material=A992)
		
		s1.addNode(0, 0, fixity='pin')
		s1.addNode(1, 1)
		s1.addNode(2, 0, fixity='roller')
		
		s1.addMember(0, 1)
		s1.addMember(1, 2)
		s1.addMember(2, 0)
		
		self.s1 = s1
		self.xs1 = xs1
		
	def tearDown(self):
		pass
	
	def test_Truss(self):
		self.assertEqual(self.s1.nMembers, 3)
		self.assertEqual(self.s1.nNodes, 3)
		self.assertEqual(self.xs1.A, self.s1.members[0].cross.A)
		self.assertEqual(self.s1.members[2].length, 2.0)	
		
	def test_solveTruss(self):
		Forces = np.matrix('0; 0; 100; 100; 0; 0')
		self.s1.directStiffness(Forces)
		#this should be a zero force member
		self.assertAlmostEqual(self.s1.members[1].axial[0], 0)
		
	def test_Frame(self):
		
		# global stiffness matrix must be symmetric
		xs1 = xs.generalSection(A=30, Ix=700)
		s1 = sc.Structure(xs1, ma.Steel(E=29000))
		s1.addNode(0, 0)
		s1.addNode(0, 144)
		s1.addNode(144, 144)
		s1.addMember(0, 1)
		s1.addMember(1, 2)
		m1 = s1.members[0]
		stiff1 = m1.T.T * m1.kframe * m1.T
		v = np.all(stiff1 == stiff1.T)
		self.assertEqual(v, True)
		
	def test_knownframe(self):
		xs1 = xs.generalSection(A=1, Ix=1)
		s1 = sc.Structure(xs1, ma.Steel(E=144))
		s1.addNode(0, 0, fixity='fixed')
		s1.addNode(24, 0, fixity='roller')
		s1.addNode(42, 0, fixity='roller')
		s1.addMember(0, 1)
		s1.addMember(2,1)
		
		
if __name__ == '__main__':
	unittest.main()
