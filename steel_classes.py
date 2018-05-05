import structural_classes as sc


class steelMember(sc.Member):
	
	
	def bucklingCapacity(self, K):
		# UNITS: inches, ksi
		elasticLimit = 4.71 * math.sqrt(self.material.E / self.material.Fy)
		slenderness = (K * 12 * self.length) / self.cross.ry
		Fe = (math.pi**2 * self.material.E) / ((K * 12* self.length) / self.cross.ry)**2
	
		if slenderness < elasticLimit:
			Fcr = (0.658**(self.material.Fy / Fe)) * self.material.Fy
		else:
			Fcr = 0.877 * Fe
		
		print('Elastic Limit %f' % elasticLimit)
		print('Fe: %f' % Fe)
		print('Fcr: %f' % Fcr)
		
		nominalCapacity = Fcr * self.cross.A
		print(nominalCapacity)
		return Fcr * self.cross.A
		
	def Lr(self):
		a = self.material.E/(0.75 * self.material.Fy)
		b = self.cross.J/(self.cross.Sx * self.cross.ho)
		return 1.95 * self.rts * a * sqrt(b + sqrt(b**2 + 6.76 * a**(-2) ))
