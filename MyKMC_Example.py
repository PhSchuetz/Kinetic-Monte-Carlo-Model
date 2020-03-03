# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 11:06:10 2020

@author: Philipp
"""
import numpy as np
#import MyKMC_OOP as MyKMC

Sim = KMC_Model(temperature=1000, DepositionRate=1.5, PulsesPerML=15,Esurf=0.5, Ebond=0.5, Eevap=5, theta=0, h0=1e7, DimSize=50, SaveFile=False, SaveVideo=False, StepDensityType = 0)

TotalPulses = 15
TimeStep = 0.05

StepDensity = [Sim.StepDensity]
Time = [Sim.Time]
IrO2Density = [Sim.IrO2Number/Sim.a.size]

PulseNumber = 1
TimeBuffer = []
StepDensityBuffer = []
IrO2Buffer=[]

while PulseNumber <= TotalPulses:	
	Sim.Deposition()
	
	TimeBuffer.append(Sim.Time)
	StepDensityBuffer.append(Sim.StepDensity)
	IrO2Buffer.append(Sim.IrO2Number/Sim.a.size)

	while Sim.Time < PulseNumber/Sim.DepositionRate: 
		Sim.Event()
		TimeBuffer.append(Sim.Time)
		StepDensityBuffer.append(Sim.StepDensity)
		IrO2Buffer.append(Sim.IrO2Number/Sim.a.size)
		
		if(TimeBuffer[-1]-TimeBuffer[0]>TimeStep):
			Time.append(np.mean(TimeBuffer))
			StepDensity.append(np.mean(StepDensityBuffer))
			IrO2Density.append(np.mean(IrO2Buffer))
			del TimeBuffer[:-1]
			del StepDensityBuffer[:-1]
			del IrO2Buffer[:-1]

#	
	PulseNumber+=1
	
Intensity=[]
Intensity[:]=[1-x for x in StepDensity]