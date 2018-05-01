import sys
sys.path.append('..')
import structural_classes as sc
import cross_sections as xs
import materials as ma

# Define material
A992 = ma.A992()

# Define cross section
xs1 = xs.generalSection(A=0.01)

# define blank structure
# the cross section and material are defaults for the members
# we will add to this later
s1 = sc.Structure(cross=xs1, material=A992)

# Add nodes to the structure
s1.addNode(0, 0, xfix=0, yfix=0)
s1.addNode(3, 4)
s1.addNode(6, 0, xfix=0, yfix=0)

# Add members to the structure
s1.addMember(0, 1)
s1.addMember(1, 2)
s1.addMember(0, 2)

s1.plot()

Forces = np.matrix('0; 0; 8660; 5000; 0; 0')
displacement = s1.directStiffness(Forces)
s1.plotDeformation(nfig=1)

Forces = np.matrix('0; 0; 0; -100000; 0; 0')
s1.directStiffness(Forces)
s1.plotDeformation(nfig=2)

s1.printMembers()
