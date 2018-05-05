import sys
sys.path.append('..')
import structural_classes as sc
import cross_sections as xs
import materials as ma

xs1 = xs.Circle(1)
ma1 = ma.A992()

s2 = sc.Structure(cross=xs1, material=ma1)

# Add nodes to the structure
s2.addNode(0, 0, xfix=0, yfix=0)
s2.addNode(2, 0, yfix=0)
s2.addNode(1, 2)
s2.addNode(2, 2, xfix=0)

# Add members to the structure
s2.addMember(0, 1)
s2.addMember(1, 2)
s2.addMember(2, 0)
s2.addMember(1, 3)
s2.addMember(2, 3)

s2.plot()
s2.printNodes()
s2.printMembers()

