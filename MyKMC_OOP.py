# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import scipy.constants as sp


class KMC_Model:
	"""
	This class describes a simple Kinetic-Monte-Carlo (KMC) simulation of 
	the growth dynamics during Pulsed Laser Deposition (PLD) of atomically-thin
	perovskite films along the (100) crystal direction with one of the two 
	possible surface terminations being unstable against oxidation and subsequent
	evaporation/volatilization. 
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
			  SaveFile=False, 
			  SaveVideo=False, 
			  StepDensityType = 0):
		
		self.temperature = temperature
		self.DepositionRate = DepositionRate
		self.PulsesPerML = PulsesPerML
		self.Esurf = Esurf
		self.Ebond = Ebond
		self.Eevap = Eevap
		self.theta = theta
		self.h0 = h0
		self.DimSize = DimSize
		self.SaveFile = SaveFile
		self.SaveVideo = SaveVideo
		self.StepDensityType = StepDensityType
		
		self.PulseDuration = 1e-7

		self.a = np.zeros((DimSize,DimSize)) #Matrix describing the crystal surface
		self.h = np.zeros((DimSize,DimSize)) #Matrix holding the hopping parameters
		self.e = np.zeros((DimSize,DimSize)) #Matrix holding the evaporation parameters
		self.gh = np.zeros(DimSize)          # group matrix
		self.ge = np.zeros(DimSize)
		self.NN = np.zeros(9)                ## NN[i]: Number of sites with i nearest-neighbours
		self.Time = 0
		self.StepDensity = self._CalcStepDensity()
		self.IrO2Number=self.a.size
		
		self.NN.fill(0)
		for i in range(np.size(self.a, 0)):
			for j in range(np.size(self.a, 1)):
				self.h[i][j]= self.h0 * self._HoppingProb(self._CheckForNN(i,j)) # compute hopping rates h[i][j] from a[i][j]
				self.e[i][j]= self.h0 * self._EvaporateProb(i,j)  # compute evaporation rates e[i][j] from a[i][j]
				self.NN[int(2*self._CheckForNN(i,j))]+=1
			self.gh[i] = np.sum(self.h[i][:])
			self.ge[i] = np.sum(self.e[i][:])  #sum up hopping and evaporation rates in each group g[i]
		self.gh = np.cumsum(self.gh)
		self.ge = np.cumsum(self.ge)
		
		
	
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
				
		self.StepDensity = self._CalcStepDensity()
		self.NN.fill(0)
		for i in range(np.size(self.a, 0)):
			for j in range(np.size(self.a, 1)):
				self.h[i][j]= self.h0 * self._HoppingProb(self._CheckForNN(i,j)) # compute hopping rates h[i][j] from a[i][j]
				self.e[i][j]= self.h0 * self._EvaporateProb(i,j)  # compute evaporation rates e[i][j] from a[i][j]
				self.NN[int(2*self._CheckForNN(i,j))]+=1
			self.gh[i] = np.sum(self.h[i][:])
			self.ge[i] = np.sum(self.e[i][:])  #sum up hopping and evaporation rates in each group g[i]
		self.gh = np.cumsum(self.gh)
		self.ge = np.cumsum(self.ge)
		
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
		if di !=0 and dj == 0: #hopping event
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
			
		elif dj!=0 and di==0: #hopping event
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
	
		elif di == 0 and dj == 0: #evaporation event
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
			self.NN[int(2*self._CheckForNN(k, l))] -= 1 #update nearest-neighbour-table
			SD_bef += self._CalcStepDensity_perSite(k, l)
	
	# Update a[i][j]    
		if di !=0 or dj != 0: #hopping event
			self.a[i % size[0], j % size[1]]-=1
			self.a[(i+di) % size[0], (j+dj) % size[1]]+=1
		
		elif di == 0 and dj == 0: #evaporation event
			self.a[i][j] -= 0.5
			self.IrO2Number -= 1  # Update IrO2Number
		gh_up = np.zeros(np.size(self.gh))
		ge_up = np.zeros(np.size(self.ge))       
	        
	# Update h[i][j]
		for kl in mh:
			k, l = kl[0], kl[1]
			h_bef = self.h[k][l]
			self.h[k][l] = self.h0 * self._HoppingProb(self._CheckForNN(k, l))
			gh_up[k] += (self.h[k][l] - h_bef)
			self.NN[int(2*self._CheckForNN(k, l))] += 1 #update nearest-neighbour-table
			SD_aft += self._CalcStepDensity_perSite(k, l)
	
	# Update e[i][j]    
		for kl in me:
			k, l = kl[0], kl[1]
			e_bef = self.e[k][l]
			self.e[k][l] = self.h0 * self._EvaporateProb(k,l)
			ge_up[k] += (self.e[k][l] - e_bef)
	
	# Update ge[i] and gh[i]     
		self.gh += np.cumsum(gh_up)
		self.ge += np.cumsum(ge_up)        
	
	# Update StepDensity   
		self.StepDensity += (SD_aft - SD_bef)
	
		if abs(np.sum(self.e) - self.ge[-1]) > 1e-4*np.sum(self.e) and self.ge[-1] > 1e-9:
		#if self.e[-1] > 1e-10:
			print("t = ",self.Time," s: Error in ge calculation!") #Recalculate ge from scratch!")
			print("np.sum(e)=",np.sum(self.e))
			print("ge[-1] =",self.ge[-1])
			print("diff = ",abs((np.sum(self.e) - self.ge[-1]))) 
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
		
	
	def _GetEventType(self):   # Checked! returns(weighted) random site
		r = np.random.random()*(self.gh[-1]+self.ge[-1]) # random number times total rate g[-1]=np.sum(h)
		if r <= self.gh[-1]:  # hopping!
			# Get Site
			i = np.searchsorted(self.gh, r, side='right') # binary chopping --> look for group with the sought-after site --> i
			if i>0: r -= self.gh[i-1]
			j = np.searchsorted(np.cumsum(self.h[i]), r) # binary chopping --> look for sought-after site within the group --> j
			
			#Get Direction
			direction=[(1,0),(-1,0),(0,1),(0,-1)]
			direc = direction[np.random.randint(0,4)] 
			di, dj = direc[0], direc[1]

		elif r > self.gh[-1] and r < self.gh[-1]+self.ge[-1]: # evaporation!
			r -= self.gh[-1]
			i = np.searchsorted(self.ge, r, side='right') # binary chopping --> look for group with the sought-after site --> i
			if i!=0: r -= self.ge[i-1]

			j = np.searchsorted(np.cumsum(self.e[i]), r) # binary chopping --> look for sought-after site within the group --> j
			di, dj = 0, 0
		return (i, j, di, dj)



	def _GetTimeStep(self):
		hop = np.array([self.h0 * self._HoppingProb(n) for n in [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]])
		Rtot= sum(self.NN*hop) + self.h0 * self.IrO2Number * np.exp(-self.Eevap*sp.e/(sp.k*self.temperature))   # GEHÖRT HIER EINE 1 HIN?
		return -np.log(np.random.random())/Rtot

	def _CalcIrO2Number(self):
		IrO2=0
		for i in range(np.size(self.a, 0)):
			for j in range(np.size(self.a, 1)):
				IrO2+=self._CalcIrO2Number_perSite(i,j)
		return IrO2

	def _CalcIrO2Number_perSite(self,i,j):
		if self.a[i][j] % 1 == 0: return 1
		else: return 0 

	def _CalcStepDensity(self):
		SD = 0
		for i in range(np.size(self.a, 0)):
			for j in range(np.size(self.a, 1)):
				SD+=self._CalcStepDensity_perSite(i,j)
		return SD

	def _CalcStepDensity_perSite(self,i,j):
		SD = 0 
		if self.StepDensityType == 0:
			SD += abs(self.a[i % np.size(self.a, 0), j % np.size(self.a, 1)] - self.a[i % np.size(self.a, 0), (j+1) % np.size(self.a, 1)])*np.cos(self.theta*np.pi/180)
			SD += abs(self.a[i % np.size(self.a, 0), j % np.size(self.a, 1)] - self.a[(i+1) % np.size(self.a, 0), j % np.size(self.a, 1)])*np.sin(self.theta*np.pi/180)
		elif self.StepDensityType == 1:
			SD += (1-self._Delta(self.a[i % np.size(self.a, 0), j % np.size(self.a, 1)], self.a[i % np.size(self.a, 0), (j+1) % np.size(self.a, 1)]))*np.cos(self.theta*np.pi/180)
			SD += (1-self._Delta(self.a[i % np.size(self.a, 0), j % np.size(self.a, 1)], self.a[(i+1) % np.size(self.a, 0), j % np.size(self.a, 1)]))*np.sin(self.theta*np.pi/180)						
		return SD/self.a.size

	def _Delta(self,a,b):
		if a==b:
			return 1
		elif a!=b:
			return 0

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

	def _HoppingProb(self,n):
		return np.exp(-(self.Esurf+n*self.Ebond)*sp.e/(sp.k*self.temperature))   

	def _EvaporateProb(self, i, j):
		if self.a[i][j] <= 0: return 0
		elif self.a[i][j] >0:
			if self.a[i][j]%1 !=0:
				return 0
			if self.a[i][j]%1 == 0:
				return np.exp(-self.Eevap*sp.e/(sp.k*self.temperature))

    
    
    
    

    
    
    
    
    

    
    
    
    
    
    
    
    
        
        




