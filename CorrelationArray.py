#This python script generates a master Correlation array file that contains a lambda,Delta t mesh as well as the Force Variance as a function of CP
#
#Steven Large
#November 22nd

import os
from math import *


def ReadFile(Path,ReadName):

	CompleteName = os.path.join(Path,ReadName)
	
	LagTime = []
	Correlation = []

	file1 = open(CompleteName,'r')
	TempData = file1.readlines()
	file1.close()

	for index in range(len(TempData)-2):
		parsed = TempData[index+2].split()

		LagTime.append(eval(parsed[0]))
		Correlation.append(eval(parsed[1]))

	return LagTime,Correlation


def WriteCorrelationArray(Path,Filename,DataArray):

	CompleteName = os.path.join(Path,Filename)

	file1 = open(CompleteName,'w')

	for index1 in range(len(DataArray)):
		for index2 in range(len(DataArray[index1])):
			file1.write('%lf\t' % DataArray[index1][index2])
		file1.write('\n')

	file1.close()


def WriteFisherLandscape(Path,Filename,DataArray):

	CompleteName = os.path.join(Path,Filename)

	file1 = open(CompleteName,'w')

	for index in range(len(DataArray)):
		file1.write('%lf\n' % DataArray[index][0])

	file1.close()


def WriteLagTime(Path,Filename,LagTime):

	CompleteName = os.path.join(Path,Filename)

	file1 = open(CompleteName,'w')

	for index in range(len(LagTime)):
		file1.write('%lf\n' % LagTime[index])

	file1.close()



ReadPath = "EquilibriumData/"
WritePath = "CorrelationMesh/"

CPValsStr = []
CPVals = [-1.00]
CP = -0.99

while CP <= 1.0:
	CPVals.append(round(CP,2))
	CP = CP + 0.01

CPVals.append(1.0)

for index in range(len(CPVals)):
	if len(str(abs(round(CPVals[index],2)))) == 4:
		CPValsStr.append(str(round(CPVals[index],2)) + str("0000"))
	else:
		CPValsStr.append(str(round(CPVals[index],2)) + str("00000"))


FilenameBase = 'BistableCorrelation_CP_'
WriteNameMesh = 'CorrelationMesh.dat'
WriteNameLagTime = 'LagTime.dat'
WriteNameCP = 'CPVals.dat'
WriteNameFisher = 'FisherInformation.dat'

CorrelationMesh = []

for index in range(len(CPVals)):

	ReadName= FilenameBase + CPValsStr[index] + '.dat'

	LagTime,Correlation = ReadFile(ReadPath,ReadName)

	CorrelationMesh.append(Correlation)


WriteCorrelationArray(WritePath,WriteNameMesh,CorrelationMesh)
WriteFisherLandscape(WritePath,WriteNameFisher,CorrelationMesh)
WriteLagTime(WritePath,WriteNameLagTime,LagTime)
WriteLagTime(WritePath,WriteNameCP,CPVals)







