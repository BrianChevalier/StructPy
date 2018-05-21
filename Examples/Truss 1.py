import StructPy.cross_sections as xs
import StructPy.structural_classes as sc
import StructPy.materials as ma

# Define material
A992 = ma.A992()

# Define cross section
xs1 = xs.IBeam(1, 1, 0.1, 0.1)

# define blank structure
# we will add to this later
s1 = sc.Structure(cross=xs1, material=A992)

# Add nodes to the structure
s1.addNode(0, 0, fixity='pin')
s1.addNode(2, 0, fixity='roller')
s1.addNode(1, 1)

# Add members to the structure
s1.addMember(0, 1)
s1.addMember(1, 2)
s1.addMember(2, 0)	

m1 = s1.members[0]
m2 = s1.members[1]
m3 = s1.members[2]

s1.plot()

