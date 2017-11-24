/* This C++ function calculates the CP values and Switch times for an Optimal step placement protocol (Using Fisher Information)  */

#include <fstream>
#include <iostream>
#include <random>
#include <cmath>
#include <string>
#include <ctime>

using namespace std;


/* Global variable declarations for Langevin integrator parameters */

double TrapStrength = 1.0;
double DampingVal = 0.25;
double beta = 1.0;
double dt = 0.1;
double mass = 1.0;


/* Global variables for Bistable well parameters */

double kL = 0.2;
double kR = 0.2;
double DeltaE = 0.0;
double X_m = 10.0;
double BistableTrap = 0.3;


/* Global variables for simulation parameters */

double DeltaCP = 20.0;
double CPMin = -10.0;
double CPMax = 10.0;
double dX = 0.01; 									//Spatial discretization parameter
double dX_Rough = 0.25; 							//Spatial disretization for autocorrelation calculation
double BoltzmannBuffer = 10.0;
double Threshold = 0.0001;							//Gradient Descent cost change threshold below which optimization is assumed to have converged
double LearningRate = 1.0;

double MAXTIME = 1000000;
double MAXSTEPS = 1500;


/* Function Prototypes */

void OptimalSpaceProtocolFisher(int NumSteps, double ProtocolTime);
void FisherInformation(double * FisherArray, double * FisherCPArray);
void ContinuousFisherInformation(double * FisherArray, double * ContinuousFisherArray, double * ContinuousCP);


/* Main OptimalStep protocol Generator Code */

int main() {
	
	double ProtocolTime = 0.1;
	int NumSteps = 5;

	//while(ProtocolTime < MAXTIME) {

	//	NumSteps = 5;

	//	while(NumSteps < MAXSTEPS) {

			OptimalSpaceProtocolFisher(NumSteps,ProtocolTime);

	//		NumSteps = NumSteps*2;

	//	}

	//	ProtocolTime = ProtocolTime*10;
	//}
}


/* Routine to generate Naive protocol sequences */

void OptimalSpaceProtocolFisher(int NumSteps, double ProtocolTime) {

	string Filename = "Protocols/OptimalSpace/OptimalSpace_FisherSteps_" + std::to_string(NumSteps) + "_Time_" + std::to_string(ProtocolTime) + ".dat";

	double * CPVals;
	CPVals = new double [NumSteps];
	double * SwitchTimes;
	SwitchTimes = new double [NumSteps];

	double StepTime = ProtocolTime/double(NumSteps - 2);
	double BufferTime = 100;

	CPVals[0] = -10.0;
	SwitchTimes[0] = BufferTime;

	for(int k = 1 ; k < NumSteps-1 ; k++) {
		SwitchTimes[k] = StepTime;
	}
	SwitchTimes[NumSteps-1] = BufferTime;

	int ArrayLength = int(DeltaCP/dX_Rough);
	int ContinuousArrayLength = int(DeltaCP/dX) + 1;

	double * FisherArray;
	FisherArray = new double [ArrayLength];
	double * CPArray;
	CPArray = new double [ArrayLength];

	double * ContinuousFisherArray;
	ContinuousFisherArray = new double [ContinuousArrayLength];
	double * CumulativeFisherArray;
	CumulativeFisherArray = new double [ContinuousArrayLength];
	double * ContinuousCP;
	ContinuousCP = new double [ContinuousArrayLength];

	double DeltaFisher;
	double FisherAcc;
	double TempCP = CPVals[0];

	for(int k = 0 ; k < ArrayLength ; k++){
		CPArray[k] = TempCP;
		TempCP += dX_Rough;
	}

	FisherInformation(FisherArray, CPArray);

	ContinuousFisherInformation(FisherArray, ContinuousFisherArray, ContinuousCP);

	FisherAcc = 0;

	for(int k = 0 ; k < ContinuousArrayLength ; k++) {
		CumulativeFisherArray[k] = FisherAcc;
		FisherAcc += ContinuousFisherArray[k];
		//printf("%lf\n", ContinuousFisherArray[k]);
	}

	DeltaFisher = CumulativeFisherArray[ContinuousArrayLength-1]/double(NumSteps - 1);
	FisherAcc = DeltaFisher;

	double FisherCost;
	double TempCost;
	int IndexLabel = 0;

	printf("\nDeltaFisher --> %lf\n\n", DeltaFisher);

	for(int k = 1 ; k < NumSteps ; k++) {
		FisherCost = 99999;

		for(int i = 1 ; i < ArrayLength ; i++) {
			TempCost = (CumulativeFisherArray[i] - FisherAcc)*(CumulativeFisherArray[i] - FisherAcc);

			if(TempCost < FisherCost){
				FisherCost = TempCost;
				IndexLabel = i;
			}
		}

		CPVals[k] = ContinuousCP[IndexLabel];
		FisherAcc += DeltaFisher;

		printf("FisherAcc --> %lf\n", FisherAcc);

	}

	std::ofstream WriteFileOptimalStepFisher;

	WriteFileOptimalStepFisher.open(Filename);
	WriteFileOptimalStepFisher << "CPValue\tDelta T\n\n";
	for(int k = 0 ; k < NumSteps ; k++) {
		WriteFileOptimalStepFisher << CPVals[k] << "\t" << SwitchTimes[k] << "\n";
	}
	WriteFileOptimalStepFisher.close();

	delete CPVals;
	delete SwitchTimes;
	delete FisherArray;
	delete ContinuousFisherArray;
	delete ContinuousCP;
	delete CumulativeFisherArray;
	delete CPArray;
}


/* Associated Numerical Routines */

void FisherInformation(double * FisherArray, double * FisherCPArray) {

	string FilenameBase = "EquilibriumData/BistableCorrelation_CP_";
	std::ifstream ReadFile;

	double CPFisher = -10.0;
	double TimeDummy;
	string dummy;

	int Counter = 0;

	while(CPFisher <= CPMax) {

		string Filename_Full = FilenameBase + std::to_string(CPFisher) + ".dat";
		ReadFile.open(Filename_Full);

		ReadFile >> dummy;
		ReadFile >> dummy;
		ReadFile >> TimeDummy;

		ReadFile >> FisherArray[Counter];
		FisherCPArray[Counter] = CPFisher;

		printf("CPFisher --> %lf\n", FisherArray[Counter]);

		Counter = Counter + 1;

		CPFisher = CPFisher + 0.25;

	}
}

void ContinuousFisherInformation(double * FisherArray, double * ContinuousFisherArray, double * ContinuousCP) {

		
	double CPFisher = -10.0;
	double ContinuousCP_Acc = -10.0;
	double Correction;
	double Slope;

	int Counter = 0;
	int OuterCounter = 0;
	int InnerCounter;

	while(CPFisher < CPMax) {

		InnerCounter = 1;

		//printf("%lf\n", FisherArray[OuterCounter]);

		while(ContinuousCP_Acc < CPFisher + dX_Rough) {
			ContinuousCP[Counter] = ContinuousCP_Acc;
			Slope = (FisherArray[OuterCounter+1] - FisherArray[OuterCounter])/dX_Rough;
			Correction = Slope*dX*InnerCounter;
			ContinuousFisherArray[Counter] = FisherArray[OuterCounter] + Correction;
			ContinuousCP_Acc += dX;
			Counter += 1;
			InnerCounter += 1;
		}

		OuterCounter += 1;
		CPFisher += dX_Rough;

	}
}




