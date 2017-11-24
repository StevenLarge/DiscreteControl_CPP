#This is the general optimization scheme for Discrete optimal protocol using a random search method
#
#Steven Large
#November 23rd 2017

import os
import random
from math import *
import matplotlib.pylot as plt

def LoadCorrelationArray(Path,Filename):
	
	CorrArray = []

	CompleteName = os.path.join(Path,Filename)
	file1 = open(CompleteName, 'r')
	TempData = file1.readlines()
	file1.close()

	for index1 in range(len(TempData)):
		CorrArray.append([])

		Parsed = TempData[index1].split()

		for index2 in range(len(TempData[index1])):
			CorrArray[index1].append(eval(Parsed[index2]))

	return CorrArray

def LoadVector(Path,Filename):

	Vector = []

	CompleteName = os.path.join(Path,Filename)
	file1 = open(CompleteName,'r')
	TempData = file1.readlines()
	file1.close()

	for index in range(len(TempData)):
		Vector.append(eval(TempData[index]))

	return Vector


def CalculateCost(TimeStep,CPStep,LambdaVals):

	return Cost


def RandomPerturbation(TimeStep,CPStep,LambdaVals):

	return TimeStepNew,CPStepNew,LambdaValsNew


def SaveProtocol(TimeStep,CPStep,LambdaVals,Path,Filename):

	CompleteName = os.path.join(Path,Filename)

	file1 = open(CompleteName,'w')

	



CorrArray = LoadCorrelationArray(ReadPath,FilenameMesh)
LagTime = LoadVector(ReadPath,FilenameTime)
CPVal = LoadVector(ReadPath,FilenameCP)

dCP = CPVal[1]-CPVal[0]
dT = LagTime[1]-LagTime[0]

NumSteps = 9 	 											#10 CP vals, but 9 steps			
TotalTime = 100
Lambda = 0

TimeStep = []
CPStep = []
LambdaIncrement = int((2.0/10.0)/dCP)
LambdaVals = [0]

for index in range(NumSteps):
	TimeStep.append((100.0/8.0)/dT) 						#8 CP values that need time allocated to them
	CPStep.append((2.0/10.0)/dCP) 				
	LambdaVals.append(int(Lambda + LambdaIncrement))

Counter = 0

Cost = CalculateCost(TimeStep,CPStep,LambdaVals)
CostTracker = []
CostTracker.append(Cost)

while Counter < 10000000:

	TimeStep,CPStep,LambdaVals = RandomPerturbation(TimeStep,CPStep,LambdaVals)

	TestCost = CalculateCost(TimeStep,CPStep,LambdaIncrement)

	if TestCost < Cost:
		Cost = TestCost
		CostTracker.append(Cost)

	Counter = Counter + 1


SaveProtocol(TimeStep,CPStep,LambdaVals)

plt.plot(CostTracker)
plt.show()
plt.close()



