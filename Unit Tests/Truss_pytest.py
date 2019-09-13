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
	assert approx(s1.members[1].axial[0], 0.000001) == 0


def from_yaml_file(filePath):
	
	## Read in yaml file as a dictionary
	with open(filePath, 'r') as stream:
		try:
			data = yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			print(exc)
	
	## Initialize structure
	xsection = xs.generalSection(**data['XSection'])
	material = ma.Custom(**data['Material'])
	s1 = Truss.Truss(cross=xsection, material=material)
	
	## Add nodes
	nodeMap = {} # map nodes to indicies
	for i, node in enumerate(data['Nodes']):
		for key, value in node.items():
			try:
				s1.addNode(value['x'], value['y'], fixity=value['fixity'])
			except KeyError:
				s1.addNode(value['x'], value['y'])
			
			nodeMap.update({key: i})
			
	## Add members
	for member in data['Members']:
		for key, value in member.items():
			SN, EN = key.split(',')
			s1.addMember(nodeMap[SN], nodeMap[EN], expectedaxial=value['axial'])
	
	return s1

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
	
def test_Truss_Ex_2_7_2():
	""" Example 2.7.2 Rajan, (page: 51)"""
	testName = 'Unit Tests/Ex_2.7.2.yaml'
	
	s1 = from_yaml_file(testName)
	s1.directStiffness(loading_from_yaml(testName))
	
	## log info
	logging.info(f'\nChecking "{testName}":')
	logging.info(f"{' '*10}| Computed | Expected")
	for i, member in enumerate(s1.members):
		
		logging.info(f'Member {i:<2} | {member.axial[0, 0]:<8.3} | {member.expectedaxial}')
		assert approx(member.axial[0], 0.001) == member.expectedaxial
	
def test_Truss_Ex_2_7_1():
	""" Example 2.7.1 Rajan, (page: 47)"""
	testName = 'Unit Tests/Ex_2.7.1.yaml'
	
	s2 = from_yaml_file(testName)
	loading = loading_from_yaml(testName)
	
	s2.directStiffness(loading)

	## log info
	logging.info(f'\nChecking "{testName}":')
	
	logging.info(loading)
	
	for node in s2.nodes:
		logging.info(node.ydef)
	
	logging.info(f"\n{' '*10}| Computed | Expected | Length")
	for i, member in enumerate(s2.members):
		logging.info(f'Member {i:<2} | {member.axial[0, 0]:<8.3} | {member.expectedaxial:<8} | {member.length:<8.3}')
        
	for i, member in enumerate(s2.members):
		assert approx(member.axial[0], 0.001) == member.expectedaxial