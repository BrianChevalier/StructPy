import StructPy.cross_sections as xs
import StructPy.structural_classes as sc
import StructPy.materials as ma
import StructPy.Truss as Truss

# Define material
A992 = ma.A992()

# Define cross section
xs1 = xs.generalSection(A=1)

# define blank structure
# the cross section and material are defaults for the members
# we will add to this later
s1 = Truss.Truss(cross=xs1, material=A992)

# Add nodes to the structure
s1.addNode(0, 0, fixity='pin')
s1.addNode(1, 1)
s1.addNode(2, 0, fixity='roller')

# Add members to the structure
s1.addMember(0, 1)
s1.addMember(1, 2)
s1.addMember(0, 2)

s1.plot()

Forces = np.matrix('0; 0; 1000; 1000; 0; 0')
disp = s1.directStiffness(Forces)
s1.plotDeformation(nfig=1, mag=1)

Forces = np.matrix('0; 0; 1000; 1000; 0; 0')
s1.directStiffness(Forces)
s1.plotDeformation(nfig=2, mag=1)

s1.printMembers()
