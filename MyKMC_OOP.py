import numpy as np
import scipy.constants as sp


class KMC_Model:
	"""
	This class describes a simple Kinetic-Monte-Carlo (KMC) simulation of 
	the growth dynamics during Pulsed Laser Deposition (PLD) of atomically-thin
	perovskite films along the (100) crystal direction with one of the two 
	possible surface terminations being unstable against oxidation and subsequent
	evaporation/volatilization. 
	
	Attributes:
	-----------
	a : numpy.ndarray, shape = (DimSize,DimSize)
		Holds the current height of the deposited crystalline film.
		
	h : numpy.ndarray, shape = (DimSize,DimSize)
		Holds the current hopping rates of each site.
		
	e : numpy.ndarray, shape = (DimSize,DimSize)
		Holds the current evaporation rates of each site.
		
	Time : float
		The currently evolved time.
	
	StepDensity : float
		The current step density.
	
	IrO2Number : float
		The current number of IrO2-terminated lattice sites.
		
	PulseDuration : float
		The time estimated for one laser pulse and the deposition of material.
	"""

	def __init__(self, 
				  DimSize = 10, 
				  temperature = 700, 
				  DepositionRate = 1.5, 
				  PulsesPerML = 20, 
				  Esurf = 0.5, 
				  Ebond = 0.5, 
				  Eevap = 5, 
				  theta = 0, 
				  h0 = 1e7, 
				  StepDensityType = 0):
		"""
		Parameters:
		-----------
		DimSize : int (default: 10)
			The size N of the NxN matrix describing the lattice.
			
		temperature : float (default: 700)
			The employed substrate temperature
			
		DepositionRate: float(default: 1.5)
			The frequency of the laser pulses/ of the deposition, given in Hz.
			
		PulsesPerML: int (default: 20)
			The average number of pulses needed to deposit one atomic monolayer 
			(ML). Determines the deposited amount of material per Deposition event.
			
		Esurf: float (default: 0.5)
			The surface energy barrier, given in eV.
			
		Ebond : float (default: 0.5)
			The nearest-neighbour binding energy barrier, given in eV.
			
		Eevap : float (default: 5)
			The energy barrier for IrO2 evaporation, given in eV.
			
		theta : float (default: 0)
			The angle under which the RHEED electrons impinge upon the surface, 
			relative to the (100)-direction, i.e., step orientation, given in 
			degrees.
			
		h0 : float (default: 1e-7)
			The attempt rate. Typically on the order of 1e-7.
			
		StepDensityType : int, 0 or 1 (default: 0)
			Specifies the formula for calculation of the step density.
			0: Use the formula by Achutharaman et al., 
			   Phys. Rev. B, 50, 8122 (1994) 
			1: Use the formula by Clarke et al., Surf. Sc. 255 91, (1991)
		"""
		self.temperature = temperature
		self.DepositionRate = DepositionRate
		self.PulsesPerML = PulsesPerML
		self.Esurf = Esurf
		self.Ebond = Ebond
		self.Eevap = Eevap
		self.theta = theta
		self.h0 = h0
		self.StepDensityType = StepDensityType
		
		#Matrix describing the crystal surface
		self.a = np.zeros((DimSize,DimSize)) 

		self.Time = 0
		self.StepDensity = self._CalcStepDensity()
		self.IrO2Number=self.a.size
		self.PulseDuration = 1e-7

		self._h = np.zeros((DimSize,DimSize)) #Matrix holding the hopping parameters
		self._e = np.zeros((DimSize,DimSize)) #Matrix holding the evaporation parameters
		self._gh = np.zeros(DimSize)          #Group matrix hopping
		self._ge = np.zeros(DimSize)			#Group Matrix evaporation
		self._NN = np.zeros(9)                #NN[i]: Number of sites with i nearest-neighbours
		
		# Initialize the private attributes _h, _e, _gh, _ge, _NN
		self._NN.fill(0)
		for i in range(np.size(self.a, 0)):
			for j in range(np.size(self.a, 1)):
				# compute hopping rates h[i][j] from a[i][j]
				self._h[i][j]= self.h0 * self._HoppingProb(self._CheckForNN(i,j))
				# compute evaporation rates e[i][j] from a[i][j] 
				self._e[i][j]= self.h0 * self._EvaporateProb(i,j)  
				# Initialize nearest-neighbours
				self._NN[int(2*self._CheckForNN(i,j))]+=1
			#sum up hopping and evaporation rates in each group
			self._gh[i] = np.sum(self._h[i][:])
			self._ge[i] = np.sum(self._e[i][:])  
		self._gh = np.cumsum(self._gh)
		self._ge = np.cumsum(self._ge)
		
		
	
	def Deposition(self):
		"""
		Calling this method will randomly increase a fraction (1/PulsesPerML)
		of the lattice sites by one. This models the pulsed deposition caused
		by one high-intensity laser pulse during PLD. The deposition is assumed
		to occur within a time span PulseDuration (default: 1e-7 s).
		"""
		#self.PulsesPerML: Number of depositions required for one nominal monolayer
		for i in range(0,np.size(self.a,0)): 
			for j in range(0,np.size(self.a,1)):
				x = np.random.random()
				if x<=1/self.PulsesPerML: self.a[i][j]+=1
		
		# Update private attributes		
		self.StepDensity = self._CalcStepDensity()
		self._NN.fill(0)
		for i in range(np.size(self.a, 0)):
			for j in range(np.size(self.a, 1)):
				self._h[i][j]= self.h0 * self._HoppingProb(self._CheckForNN(i,j)) # compute hopping rates h[i][j] from a[i][j]
				self._e[i][j]= self.h0 * self._EvaporateProb(i,j)  # compute evaporation rates e[i][j] from a[i][j]
				self._NN[int(2*self._CheckForNN(i,j))]+=1
			self._gh[i] = np.sum(self._h[i][:])
			self._ge[i] = np.sum(self._e[i][:])  #sum up hopping and evaporation rates in each group g[i]
		self._gh = np.cumsum(self._gh)
		self._ge = np.cumsum(self._ge)
		
		self.Time += self.PulseDuration
		

	def Event(self):
		"""
		Calling this method will increment the simulation by one event. 
		The event can either be Diffusion of one perovskite ABO3 unit cell or the 
		Evaporation of a half unit cell (BO2), resulting in a local AO-termination.
		The type of event and the lattice site is determined randomnly.
		"""
		(i, j, di, dj) = self._GetEventType()
		#print(i,j,di,dj)
		mh, me =[], []
		size = (np.shape(self.a))
		#hopping event, left/right
		if di !=0 and dj == 0: 
			mh.append([(i-di)% size[0],    j% size[1]])
			mh.append([i% size[0],     (j-1)% size[1]])
			mh.append([i% size[0],         j% size[1]])
			mh.append([i% size[0],     (j+1)% size[1]])
			mh.append([(i+di)% size[0],(j-1)% size[1]])
			mh.append([(i+di)% size[0],    j% size[1]])
			mh.append([(i+di)% size[0],(j+1)% size[1]])
			mh.append([(i+2*di)% size[0],  j% size[1]])
	
			me.append([i% size[0],          j% size[1]])
			me.append([(i+di)% size[0],     j% size[1]])        
		
		#hopping event, up/down
		elif dj!=0 and di==0: 
			mh.append([i% size[0],    (j-dj)% size[1]])
			mh.append([(i-1)% size[0],     j% size[1]])
			mh.append([i% size[0],         j% size[1]])
			mh.append([(i+1)% size[0],     j% size[1]])
			mh.append([(i-1)% size[0],(j+dj)% size[1]])
			mh.append([i% size[0],    (j+dj)% size[1]])
			mh.append([(i+1)% size[0],(j+dj)% size[1]])
			mh.append([i% size[0],  (j+2*dj)% size[1]])
			
			me.append([i% size[0],         j% size[1]])
			me.append([i% size[0],    (j+dj)% size[1]])
		
		#evaporation event
		elif di == 0 and dj == 0: 
			mh.append([i% size[0],         j% size[1]])
			mh.append([(i+1)% size[0],     j% size[1]])
			mh.append([(i-1)% size[0],     j% size[1]])
			mh.append([i% size[0],     (j+1)% size[1]])
			mh.append([i% size[0],     (j-1)% size[1]])
			
			me.append([i% size[0],         j% size[1]])
	
		mh.sort(key=lambda tup:tup[0]) # m: List of entries that are changed in h[i][j], where i = mh[][0], j=mh[][1]
		me.sort(key=lambda tup:tup[0]) # m: List of entries that are changed in e[i][j], where i = me[][0], j=me[][1]
	
	# Save NN and StepDensity before Update of a[i[j]]                
		SD_bef, SD_aft = 0, 0
		for kl in mh:
			k, l = kl[0], kl[1]
			self._NN[int(2*self._CheckForNN(k, l))] -= 1 #update nearest-neighbour-table
			SD_bef += self._CalcStepDensity_perSite(k, l)
	
	# Update a[i][j]    
		if di !=0 or dj != 0: #hopping event
			self.a[i % size[0], j % size[1]]-=1
			self.a[(i+di) % size[0], (j+dj) % size[1]]+=1
		
		elif di == 0 and dj == 0: #evaporation event
			self.a[i][j] -= 0.5
			self.IrO2Number -= 1  # Update IrO2Number
		gh_up = np.zeros(np.size(self._gh))
		ge_up = np.zeros(np.size(self._ge))       
	        
	# Update h[i][j]
		for kl in mh:
			k, l = kl[0], kl[1]
			h_bef = self._h[k][l]
			self._h[k][l] = self.h0 * self._HoppingProb(self._CheckForNN(k, l))
			gh_up[k] += (self._h[k][l] - h_bef)
			self._NN[int(2*self._CheckForNN(k, l))] += 1 #update nearest-neighbour-table
			SD_aft += self._CalcStepDensity_perSite(k, l)
	
	# Update e[i][j]    
		for kl in me:
			k, l = kl[0], kl[1]
			e_bef = self._e[k][l]
			self._e[k][l] = self.h0 * self._EvaporateProb(k,l)
			ge_up[k] += (self._e[k][l] - e_bef)
	
	# Update ge[i] and gh[i]     
		self._gh += np.cumsum(gh_up)
		self._ge += np.cumsum(ge_up)        
	
	# Update StepDensity   
		self.StepDensity += (SD_aft - SD_bef)
	
		if abs(np.sum(self._e) - self._ge[-1]) > 1e-4*np.sum(self._e) and self._ge[-1] > 1e-9:
		#if self._e[-1] > 1e-10:
			print("t = ",self.Time," s: Error in ge calculation!") #Recalculate ge from scratch!")
			print("np.sum(e)=",np.sum(self._e))
			print("ge[-1] =",self._ge[-1])
			print("diff = ",abs((np.sum(self._e) - self._ge[-1]))) 
			raise Exception("Error in Calculation")  
	# Update Time
		self.Time += self._GetTimeStep()
	
	
	
	def PrintInfo(self):
		"""
		Returns a string with the most relevant physical parameters describing
		the KMC simulation.
		"""
		info = "T = "+"{:1.0f}".format(self.temperature)+", "
		info += "DepRate="+"{:1.1f}".format(self.DepositionRate)+", "
		info += "PulsesperML="+"{:1.0f}".format(self.PulsesPerML)+", "
		info += "Esurf="+"{:1.2f}".format(self.Esurf)+", "
		info += "Ebond="+"{:1.2f}".format(self.Ebond)+", "
		info += "Eevap="+"{:1.2f}".format(self.Eevap)+", "
		info += "N="+"{:1.0f}".format(self.DimSize)+", "
		info += "h0="+"{:.1E}".format(self.h0)
		return info
		
	# Monte-Carlo-Step: Create random numbers and return event, site, and direction
	def _GetEventType(self):   # Checked! returns(weighted) random site
		r = np.random.random()*(self._gh[-1]+self._ge[-1]) # random number times total rate g[-1]=np.sum(h)
		if r <= self._gh[-1]:  # hopping!
			# Get Site
			i = np.searchsorted(self._gh, r, side='right') # binary chopping --> look for group with the sought-after site --> i
			if i>0: r -= self._gh[i-1]
			j = np.searchsorted(np.cumsum(self._h[i]), r) # binary chopping --> look for sought-after site within the group --> j
			
			#Get Direction
			direction=[(1,0),(-1,0),(0,1),(0,-1)]
			direc = direction[np.random.randint(0,4)] 
			di, dj = direc[0], direc[1]

		elif r > self._gh[-1] and r < self._gh[-1]+self._ge[-1]: # evaporation!
			r -= self._gh[-1]
			i = np.searchsorted(self._ge, r, side='right') # binary chopping --> look for group with the sought-after site --> i
			if i!=0: r -= self._ge[i-1]

			j = np.searchsorted(np.cumsum(self._e[i]), r) # binary chopping --> look for sought-after site within the group --> j
			di, dj = 0, 0
		return (i, j, di, dj)


	# Calculate the discrete time evolution increment
	def _GetTimeStep(self):
		hop = np.array([self.h0 * self._HoppingProb(n) for n in [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]])
		Rtot= sum(self._NN*hop) + self.h0 * self.IrO2Number * np.exp(-self.Eevap*sp.e/(sp.k*self.temperature))   # GEHÖRT HIER EINE 1 HIN?
		return -np.log(np.random.random())/Rtot
	
	# Calculates the total number of IrO2-terminated sites
	def _CalcIrO2Number(self):
		IrO2=0
		for i in range(np.size(self.a, 0)):
			for j in range(np.size(self.a, 1)):
				IrO2+=self._CalcIrO2Number_perSite(i,j)
		return IrO2
	
	# Checks whether the site (i,j) is IrO2-terminated
	def _CalcIrO2Number_perSite(self,i,j):
		if self.a[i][j] % 1 == 0: return 1
		else: return 0 
	
	# Calculates the total step density
	def _CalcStepDensity(self):
		SD = 0
		for i in range(np.size(self.a, 0)):
			for j in range(np.size(self.a, 1)):
				SD+=self._CalcStepDensity_perSite(i,j)
		return SD
	
	# Calculates the number of steps at site (i,j)
	def _CalcStepDensity_perSite(self,i,j):
		SD = 0 
		if self.StepDensityType == 0:
			SD += abs(self.a[i % np.size(self.a, 0), j % np.size(self.a, 1)] - self.a[i % np.size(self.a, 0), (j+1) % np.size(self.a, 1)])*np.cos(self.theta*np.pi/180)
			SD += abs(self.a[i % np.size(self.a, 0), j % np.size(self.a, 1)] - self.a[(i+1) % np.size(self.a, 0), j % np.size(self.a, 1)])*np.sin(self.theta*np.pi/180)
		elif self.StepDensityType == 1:
			SD += (1-self._Delta(self.a[i % np.size(self.a, 0), j % np.size(self.a, 1)], self.a[i % np.size(self.a, 0), (j+1) % np.size(self.a, 1)]))*np.cos(self.theta*np.pi/180)
			SD += (1-self._Delta(self.a[i % np.size(self.a, 0), j % np.size(self.a, 1)], self.a[(i+1) % np.size(self.a, 0), j % np.size(self.a, 1)]))*np.sin(self.theta*np.pi/180)						
		return SD/self.a.size
	
	# Helper function for _CalcStepDensity_perSite
	def _Delta(self,a,b):
		if a==b:
			return 1
		elif a!=b:
			return 0
	# Calculates the number of nearest neighbors at site (i,j)
	def _CheckForNN(self,i,j):
		n=0
		m=[]
		m.append(((i-1) % np.size(self.a, 0), j % np.size(self.a, 1)))
		m.append(((i+1) % np.size(self.a, 0), j % np.size(self.a, 1)))
		m.append((i % np.size(self.a, 0), (j-1) % np.size(self.a, 1)))
		m.append((i % np.size(self.a, 0), (j+1) % np.size(self.a, 1)))

		for kl in m:
			if self.a[i, j] - self.a[kl] <= 0:
				n += 1
			elif self.a[i, j] - self.a[kl] == 0.5:
				n += 0.5
		return n  
	
	# Calculates the Arrhenius-type hopping rate/propability of site (i,j)
	def _HoppingProb(self,n):
		return np.exp(-(self.Esurf+n*self.Ebond)*sp.e/(sp.k*self.temperature))   

	# Calculates the Arrhenius-type evaporation rate of site (i,j) - if it is IrO2-terminated
	def _EvaporateProb(self, i, j):
		if self.a[i][j] <= 0: return 0
		elif self.a[i][j] >0:
			if self.a[i][j]%1 !=0:
				return 0
			if self.a[i][j]%1 == 0:
				return np.exp(-self.Eevap*sp.e/(sp.k*self.temperature))

    
    
    
    

    
    
    
    
    

    
    
    
    
    
    
    
    
        
        




