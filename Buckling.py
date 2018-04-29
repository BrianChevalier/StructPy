import cross_sections as xs
import structural_classes as sc
import materials as ma


xs1 = xs.generalSection(A=30, Ix=700)

s1 = sc.Structure(xs1, ma.Steel(E=29000))

s1.addNode(0, 0)
s1.addNode(0, 144)
s1.addNode(144, 144)

s1.addMember(0, 1)
s1.addMember(1, 2)

m1 = s1.members[0]
stiff1 = m1.T.T * m1.kframe * m1.T

s1.plot()
