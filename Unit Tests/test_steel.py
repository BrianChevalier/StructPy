

# TODO: Make steel members. Create unit tests.

def test_bucklingCapacity(self):
		xs1 = xs.AISC('W21X44')
		s1 = sc.Structure(xs1, ma.Steel(Fy=60))
		s1.addNode(0, 0)
		s1.addNode(0, 15)
		s1.addMember(0, 1)
		member = s1.members[0]
		# almost equal because rounding.
		self.assertAlmostEqual(member.bucklingCapacity(1), 143.906124819)
		## test another example
		xs1 = xs.AISC('W16X67')
		s1 = sc.Structure(xs1, ma.Steel(Fy=60))
		s1.addNode(0, 0)
		s1.addNode(0, 15)
		s1.addMember(0, 1)
		member = s1.members[0]
		self.assertAlmostEqual(member.bucklingCapacity(1), 661.66119741611089)
	
