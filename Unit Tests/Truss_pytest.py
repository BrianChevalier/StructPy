"""
This is the StructPy Truss testing suite. This uses pytest and logs information to `Unit Tests/Truss_pytest.log`. All functions beginning with `test_` are functions run by pytest, other functions aid in testing.

There are several example trusses in .yaml files. These are file formatted to store structure information. They are used for the purpose of easily testing many known solutions.
"""

from pytest import approx
import yaml
import logging

logging.basicConfig(filename='Unit Tests/Truss_pytest.log',
					level=logging.INFO,
				    format='%(message)s')
logging.info(f'Running: {__name__}')

import StructPy.cross_sections as xs
import StructPy.structural_classes as sc
import StructPy.materials as ma
import StructPy.Truss as Truss
import numpy as np

def make_structure():
	xs1 = xs.Rectangle(2, 2)
	A992 = ma.A992()
	s1 = Truss.Truss(cross=xs1, material=A992)
	
	s1.addNode(0, 0, fixity='pin')
	s1.addNode(1, 1)
	s1.addNode(2, 0, fixity='roller')
	s1.addMember(0, 1)
	s1.addMember(1, 2)
	s1.addMember(2, 0)
	return (s1, xs1)

def test_struct():
	s1, xs1 = make_structure()
	assert s1.nMembers == 3
	assert s1.nNodes == 3
	assert xs1.A == s1.members[0].cross.A
	assert s1.members[2].length == 2.0

def test_solveTruss():
	s1, xs1 = make_structure()
	Forces = np.array([0, 0, 100, 100, 0, 0])
	s1.directStiffness(Forces)
	
	#this should be a zero force member
	assert approx(s1.members[1].axial, 0.000001) == 0


def loading_from_yaml(filePath):
	with open(filePath, 'r') as stream:
		try:
			data = yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			print(exc)
	
	loading = np.zeros(2*len(data['Nodes']))
	nodeMap = {} # map nodes to indicies
	
	for index, node in enumerate(data['Nodes']):
		for name, value in node.items():
			nodeMap.update({name: index})
	
	for load in data['Loads']:
		for nodelabel, value in load.items():
			loading[2*nodeMap[nodelabel]] = value['x']
			loading[2*nodeMap[nodelabel]+1] = value['y']
	
	return loading
	

def Truss_Test_From_File(fileName):
	""" Test Example Yaml Files """
	
	s2 = Truss.Truss.from_yaml_file(fileName)
	loading = loading_from_yaml(fileName)
	s2.directStiffness(loading)
        
	for i, member in enumerate(s2.members):
		assert approx(member.axial, 0.001) == member.expectedaxial


def test_2_7_1():
	Truss_Test_From_File('Unit Tests/Ex_2.7.1.yaml')

def test_2_7_2():
	Truss_Test_From_File('Unit Tests/Ex_2.7.2.yaml')

def test_6_2_4():
	Truss_Test_From_File('Unit Tests/Ex_6.2.4.yaml')