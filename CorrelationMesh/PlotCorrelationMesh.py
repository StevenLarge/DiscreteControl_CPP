#This python script plots the correlation mesh plot in Matplotlib
#
#Steven Large
#November 22nd 2017


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

def ReadArray(Filename):

	CorrelationArray = []

	file1 = open(Filename,'r')
	TotalData = file1.readlines()
	file1.close()

	for index1 in range(len(TotalData)):
		Parsed = TotalData[index1].split()
		CorrelationArray.append([])
		for index2 in range(len(Parsed)):
			CorrelationArray[index1].append(eval(Parsed[index2]))

	return CorrelationArray


def ReadVector(Filename):

	Data = []

	file1 = open(Filename,'r')
	TempData = file1.readlines()
	file1.close()

	for index in range(len(TempData)):
		Parsed = TempData[index].split()
		Data.append(eval(Parsed[0]))

	return Data


ReadNameArray = "CorrelationMesh.dat"
ReadNameTime = "LagTime.dat"
ReadNameCP = "CPVals.dat"

CorrArray = ReadArray(ReadNameArray)
LagTime = ReadVector(ReadNameTime)
CPVals = ReadVector(ReadNameCP)


LagTimeNP = np.zeros(int(len(LagTime)/2))
CPValsNP = np.zeros(len(CPVals))

print("Check 1 --- \n\n")

for index in range(int(len(LagTime)/2)):
	LagTimeNP[index] = LagTime[index]

print("Check 2 --- \n\n")

for index in range(len(CPVals)):
	CPValsNP[index] = CPVals[index]

print("Check 3 --- \n\n")

CorrMeshNP = np.zeros((len(CorrArray),int(len(CorrArray[0])/2)))

print("Check 4 --- \n\n")

for index1 in range(len(CorrArray)):
	for index2 in range(int(len(CorrArray[0])/2)):
		CorrMeshNP[index1,index2] = CorrArray[index1][index2]

print("Check 5 --- \n\n")

LagTimeNP,CPValsNP = np.meshgrid(LagTimeNP,CPValsNP)

print("Check 6 --- \n\n")

print np.shape(LagTimeNP)
print np.shape(CPValsNP)
print np.shape(CorrMeshNP)

fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(LagTimeNP,CPValsNP,CorrMeshNP, rstride=1, cstride=1, cmap='plasma')

hfont = {'fontname':'Times New Roman'}

ax.text2D(0.30,0.95, "Autocorrelation Mesh", transform=ax.transAxes, fontsize=20, **hfont)

ax.set_xlim(0,250)
ax.set_ylim(-1,1)
ax.set_zlim(0,3.5)

ax.set_xlabel('Lag Time', fontsize=16, **hfont)
ax.set_ylabel('Trap Minima', fontsize=16, **hfont)
ax.set_zlabel('Force Autocorrelation', fontsize=16, **hfont)

plt.savefig('AutoCorrMesh.pdf',format='pdf')

plt.show()
plt.close()



