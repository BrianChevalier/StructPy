import StructPy.cross_sections as xs
import StructPy.structural_classes as sc
import StructPy.materials as ma

xs1 = xs.Circle(1)
ma1 = ma.A992()

s2 = sc.Structure(cross=xs1, material=ma1)

# Add nodes to the structure
s2.addNode(0, 0, fixity='pin')
s2.addNode(2, 0, fixity='roller')
s2.addNode(1, 2)
s2.addNode(2, 2, fixity='wallslider')

# Add members to the structure
s2.addMember(0, 1)
s2.addMember(1, 2)
s2.addMember(2, 0)
s2.addMember(1, 3)
s2.addMember(2, 3)

s2.plot()
s2.printNodes()
s2.printMembers()

