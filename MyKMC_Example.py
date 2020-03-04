"""
This is an example how to simulate the growth dynamics using the KMC_Model class
"""
import numpy as np
import matplotlib.pyplot as plt

Example = KMC_Model(DimSize=200, 
					temperature=1000, 
					DepositionRate=1.5, 
					PulsesPerML=15,
					Esurf=0.5, 
					Ebond=0.5, 
					Eevap=5, 
					theta=0, 
					h0=1e7, 
					StepDensityType = 0)

# Number of total laer pulses in this simulation
TotalPulses = 120
# Time resolution of the CCD camera: TimeStep
# Average the step density over this time span
TimeStep = 0.05

# Lists that will contain the time, step density and IrO2 density
StepDensity = [Example.StepDensity]
Time = [Example.Time]
IrO2Density = [Example.IrO2Number/Example.a.size]

# Buffer lists for averaging over TimeStep
TimeBuffer = []
StepDensityBuffer = []
IrO2Buffer=[]

# Set counter of laser pulses to 1
PulseNumber = 1
while PulseNumber <= TotalPulses:	
	# Deposition of material
	Example.Deposition()
	
	# Write time, step density and IrO2-density in buffer lists
	TimeBuffer.append(Example.Time)
	StepDensityBuffer.append(Example.StepDensity)
	IrO2Buffer.append(Example.IrO2Number/Example.a.size)
	
	# Start time evolution between laser pulses, i.e., deposition events
	while Example.Time < PulseNumber/Example.DepositionRate: 
		# Diffusion or evaporation event
		Example.Event()
		# Write time, step density and IrO2-density in buffer lists
		TimeBuffer.append(Example.Time)
		StepDensityBuffer.append(Example.StepDensity)
		IrO2Buffer.append(Example.IrO2Number/Example.a.size)
		
		# Calculate average values of time, step density and IrO2 density if
		# the evolved time is equal to TimeStep
		if(TimeBuffer[-1]-TimeBuffer[0]>TimeStep):
			Time.append(np.mean(TimeBuffer))
			StepDensity.append(np.mean(StepDensityBuffer))
			IrO2Density.append(np.mean(IrO2Buffer))
			del TimeBuffer[:-1]
			del StepDensityBuffer[:-1]
			del IrO2Buffer[:-1]
	# Increase Pulse counter
	PulseNumber+=1

# Calculate RHEED intensity from step density
Intensity=[]
Intensity[:]=[1-x for x in StepDensity]

# Plot evolution of Rheed intensity and IrO2-density
plt.plot(Time,Intensity, Time, IrO2Density)
plt.xlabel("Time (s)")
plt.ylabel("Simulated RHEED intensity and IrO2 density")