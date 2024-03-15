import numpy as N

class Pipe_material():
	def __init__(self, Tmin, Tmax, force_T_limits=False):
		self.Tmin = Tmin
		self.Tmax = Tmax
		self.force_T_limits = force_T_limits

	@staticmethod
	def check_valid(prop):
		def wrapper_check_valid(self, T):
			if self.force_T_limits==False:
				if hasattr(T,'__len__'):
					Tlow = (T<self.Tmin).any()
					Thigh = (T>self.Tmax).any()
				else:
					Tlow = (T<self.Tmin)
					Thigh = (T>self.Tmax)
				if Tlow or Thigh:
					print ("Temperature outside of correlations range")
					return N.nan
			return prop(self, T)
		return wrapper_check_valid

class SS316(Pipe_material):
    '''
    Stainless steel 316 properties. So far only thermal conductivity extracted from the data in: http://www.inductor-jmag.ru/files/content/a129160.pdf
    '''
    def __init__(self, Tmin=273., Tmax=1600.):
        Pipe_material.__init__(self, Tmin, Tmax)
    #@Pipe_material.check_valid
    def k(self, T):
        return -2E-06*T**2. + 0.0179*T + 8.3005
    #@Pipe_material.check_valid
    def rho_E(self, T): # 10-8 ome m
        return -2.66745345e-05*T**2. + 9.10647450e-02*T + 5.20801141e+01

class Rod(Pipe_material):
    '''
    Assumed Rod material
    '''
    def __init__(self, Tmin=273., Tmax=1600.):
        Pipe_material.__init__(self, Tmin, Tmax)

    def rho_E(self, T): # 10-8 ome m
        return 18e-8
    
    def k(self, T):
        return 45 # W/m•K
    
    def Cp(self, T):
        return 0.711e3 # J/(kgK)
    
    def rho(self, T):
        return 7800 #kg/m3
    
class Rod_T(Pipe_material):
    '''
    Assumed Rod material
    '''
    def __init__(self, Tmin=273., Tmax=1600.):
        Pipe_material.__init__(self, Tmin, Tmax)

    def rho_E(self, T): # 10-8 ome m
        return 18e-8*(1+T*0.4/1200)
    
    def k(self, T):
        return 45*(1-T*0.3/1200) # W/m•K
    
    def Cp(self, T):
        return 0.711e3*(1+T*0.3/1200) # J/(kgK)
    
    def rho(self, T):
        return 7800*(1-T*0.2/1200) #kg/m3
    
    def rho_0(self):
        return 18e-8
    
    def phi_0(self):
        return 0.4/1200
    
    