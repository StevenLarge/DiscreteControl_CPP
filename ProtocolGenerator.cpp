/* This C++ function calculates the CP values and Switch times for each of the protocol classes */

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

double MAXTIME = 100000;
double MAXSTEPS = 1000;



/* Function Prototypes */

void NaiveProtocol(int NumSteps, double ProtocolTime);
void OptimalSpaceProtocol(int NumSteps, double ProtocolTime);
void OptimalSpaceProtocolFisher(int NumSteps, double ProtocolTime);
void OptimalTimeProtocol(int NumSteps, double ProtocolTime);
void OptimalProtocol(int NumSteps, double ProtocolTime);
double RelativeEntropy(double CP1, double CP2);
void FisherInformation(double * FisherArray);
void LoadCorrelationArray(double * CorrelationArray, double * CPVals);
void OptimizeSwitchTimes(double * SwitchTimes, double * CorrelationArray, double ProtocolTime, int NumTimes);
double GradientDescentTime(double * SwitchTimes, double * CorrelationArray, double ProtocolTime);
void ContinuousFisherInformation(double * FisherArray, double * ContinuousFisherArray, double * ContinuousCP);



/* Main Driving Routine */

int main(){

	double ProtocolTime = 0.1;
	int NumSteps = 5;


	while(ProtocolTime < MAXTIME) {

		while(NumSteps < MAXSTEPS) {

			NaiveProtocol(NumSteps,ProtocolTime);
			OptimalSpaceProtocolFisher(NumSteps,ProtocolTime);
			//OptimalSpaceProtocol(NumSteps,ProtocolTime);
			OptimalTimeProtocol(NumSteps,ProtocolTime);
			//OptimalProtocol(NumSteps,ProtocolTime);

			NumSteps = NumSteps*2;

		}

		ProtocolTime = ProtocolTime*10;

	}

}



/* Protocol Generator Routines */

void NaiveProtocol(int NumSteps, double ProtocolTime) {

	string Filename = "Protocols/Naive_Steps" + std::to_string(NumSteps) + "_Time" + std::to_string(ProtocolTime) + ".dat";

	double * CPVals;
	CPVals = new double [NumSteps];
	double * SwitchTimes; 									//Note that these SwitchTimes are Delta quantities, not absolute
	SwitchTimes = new double [NumSteps];

	double StepSize = DeltaCP/double(NumSteps - 1);
	double StepTime = ProtocolTime/double(NumSteps - 2);
	double BufferTime = 100;

	double CurrentCP = StepSize;

	CPVals[0] = -10.0;
	SwitchTimes[0] = BufferTime;

	for(int k = 1 ; k < NumSteps ; k++) {
		CPVals[k] = CurrentCP;
		CurrentCP += StepSize;
	}

	for(int k = 1 ; k < NumSteps-1 ; k++) {
		SwitchTimes[k] = StepTime;
	}
	SwitchTimes[NumSteps-1] = BufferTime;

	std::ofstream WriteFileNaive;

	WriteFileNaive.open(Filename);
	WriteFileNaive << "CPValue\tDelta T\n\n";
	for(int k = 0 ; k < NumSteps ; k++) {
		WriteFileNaive << CPVals[k] << "\t" << SwitchTimes[k] << "\n";
	}
	WriteFileNaive.close();
}


void OptimalSpaceProtocol(int NumSteps, double ProtocolTime) {

	string Filename = "Protocols/OptimalSpace_Steps" + std::to_string(NumSteps) + "_Time" + std::to_string(ProtocolTime) + ".dat";

	double * CPVals;
	CPVals = new double [NumSteps];
	double * SwitchTimes; 									//Note that these SwitchTimes are Delta quantities, not absolute
	SwitchTimes = new double [NumSteps - 1];
}


void OptimalSpaceProtocolFisher(int NumSteps, double ProtocolTime) {

	string Filename = "Protocols/OptimalSpace_FisherSteps" + std::to_string(NumSteps) + "_Time" + std::to_string(ProtocolTime) + ".dat";

	double * CPVals;
	CPVals = new double [NumSteps];
	double * SwitchTimes;
	SwitchTimes = new double [NumSteps - 1];

	double StepTime = ProtocolTime/double(NumSteps - 2);
	double BufferTime = 100;

	CPVals[0] = -10.0;
	SwitchTimes[0] = BufferTime;

	for(int k = 1 ; k < NumSteps-1 ; k++) {
		SwitchTimes[k] = StepTime;
	}
	SwitchTimes[NumSteps-1] = BufferTime;

	double ArrayLength = int(DeltaCP/dX_Rough);
	double ContinuousArrayLength = int(DeltaCP/dX);

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
	double TempCP = -10.0;

	for(int k = 0 ; k < ArrayLength ; k++){
		CPArray[k] = TempCP;
		TempCP += dX_Rough;
	}

	FisherInformation(FisherArray, CPArray);

	ContinuousFisherInformation(FisherArray, ContinuousFisherArray, ContinuousCP);

	for(int k = 0 ; k < ContinuousArrayLength ; k++) {
		for(int i = 0 ; i <= k ; i++){
			CumulativeFisherArray[i] += ContinuousFisherArray[i];
		}
	}

	DeltaFisher = CumulativeFisherArray[ContinuousArrayLength-1]/double(NumSteps - 1);
	FisherAcc = DeltaFisher;

	double FisherCost;
	double IndexLabel = 0;

	for(int k = 1 ; k < NumSteps ; k++) {
		FisherCost = 99999;

		for(int i = 1 ; i < ArrayLength ; i++) {
			TempCost = CumulativeFisherArray[i] - FisherAcc;

			if(TempCost < FisherCost){
				FisherCost = TempCost;
				IndexLabel = i;
			}
		}

		CPVal[k] = ContinuousCP[IndexLabel];
		FisherAcc += DeltaFisher;

	}

	std::ofstream WriteFileOptimalStepFisher;

	WriteFileOptimalStepFisher.open(Filename);
	WriteFileOptimalStepFisher << "CPValue\tDelta T\n\n";
	for(int k = 0 ; k < NumSteps ; k++) {
		WriteFileOptimalStepFisher << CPVals[k] << "\t" << SwitchTimes[k] << "\n";
	}
	WriteFileOptimalStepFisher.close();
}


void OptimalTimeProtocol(int NumSteps, double ProtocolTime) {

	string Filename = "Protocols/OptimalTime_Steps" + std::to_string(NumSteps) + "_Time" + std::to_string(ProtocolTime) + ".dat";

	double * CPVals;
	CPVals = new double [NumSteps];
	double * SwitchTimes; 									//Note that these SwitchTimes are Delta quantities, not absolute
	SwitchTimes = new double [NumSteps - 1];

	double StepSize = DeltaCP/double(NumSteps - 1);

	double CurrentCP = StepSize;

	CPVals[0] = -10.0;

	for(int k = 1 ; k < NumSteps ; k++) {
		CPVals[k] = CurrentCP;
		CurrentCP += StepSize;
	}

	int LagSteps = int(MAXLAG/dt);
	int ArrayLength = int(DeltaCP/dX);

	double * CorrelationArray;
	CorrelationArray = new double [LagSteps][ArrayLength];

	for(int k = 0 ; k < LagSteps ; k++) {
		for(int i = 0 ; i < ArrayLength ; i++) {
			CorrelationArray[k][i] = 0.0;
		}
	}

	LoadCorrelationArray(CorrelationArray,CPVals);

	OptimizeSwitchTimes(SwitchTimes, CorrelationArray, ProtocolTime, NumSteps);

	std::ofstream WriteFileOptimalStep;

	WriteFileOptimalStep.open(Filename);
	WriteFileOptimalStep << "CPValue\tDelta T\n\n";
	for(int k = 0 ; k < NumSteps ; k++) {
		WriteFileOptimalStep << CPVals[k] << "\t" << SwitchTimes[k] << "\n";
	}
	WriteFileOptimalStep.close();
}


void OptimalProtocol(int NumSteps, double ProtocolTime) {

	string Filename = "Protocols/Optimal_Steps" + std::to_string(NumSteps) + "_Time" + std::to_string(ProtocolTime) + ".dat";

	double * CPVals;
	CPVals = new double [NumSteps];
	double * SwitchTimes; 									//Note that these SwitchTimes are Delta quantities, not absolute
	SwitchTimes = new double [NumSteps-1];
}



/* Auxillary Numerical routines */

double RelativeEntropy(double CPVal1, double CPVal2) {

	double ArrayLength = (CPRange + 2*BolzmannBuffer)/dX;

	double Partition = 0.0;
	double Energy1;
	double Energy2;

	double * Boltzmann;
	Boltzmann = new double [ArrayLength];
	double * DeltaEnergyArray;
	DeltaEnergyArray = new double [ArrayLength];

	double * RelativeEntropy;
	RelativeEntropy = new double [ArrayLength];

	double CP = CPMin - BoltzmannBuffer;

	for(int k = 0 ; k < ArrayLength ; k++) {
		Energy1 = TotalEnergy(CPVal1);
		Energy2 = TotalEnergy(CPVal2);
		Boltzmann[k] = exp(-beta*Energy1);
		DeltaEnergyArray[k] = Energy2 - Energy1;
		Partition += Boltzmann[k];
		CP += dX;
	}

	for(int k = 0 ; k < ArrayLength ; k++) {
		Boltzmann[k] = Boltzmann[k]/Partition;
	}

	for(int k = 0 ; k < ArrayLength ; k++) {
		RelativeEntropy += DeltaEnergyArray[k]*Boltzmann[k]*dX;
	}

	return RelativeEntropy;
}

void FisherInformation(double * FisherArray, double * FisherCPArray) {

	string FilenameBase = "BistableCorrelation_CP_";
	std::ifstream ReadFile;

	double CPFisher = -10.0;
	double TimeDummy;
	string dummy;

	int Counter = 0;

	while(CPFisher <= CPMax) {

		string Filename_Full = FilenameBase + std::to_string(CPFisher) + ".dat";
		ReadFile.open(Filename_Fill);

		ReadFile >> dummy;
		ReadFile >> dummy;
		ReadFile >> TimeDummy;

		ReadFile >> FisherArray[Counter];
		FisherCPArray[Counter] = CPFisher;

		Counter = Counter + 1;

		CPFisher = CPFisher + 0.25;

	}
}

void ContinuousFisherInformation(double * FisherArray, double * ContinuousFisherArray, double * Continuous CP) {

		
	double CPFisher = -10.0;
	double ContinuousCP_Acc = -10.0;
	double Correction;

	int Counter = 0;
	int OuterCounter = 0;

	while(CPFisher < CPMax) {

		while(ContinuousCP_Acc < CPFisher + dX_Rough) {
			ContinuousCP[Counter] = ContinuousCP_Acc;
			Slope = (FisherArray[OuterCounter+1] - FisherArray[OuterCounter])/dX_Rough;
			Correction = Slope*dX;
			ContinuousFisherArray[Counter] = FisherArray[OuterCounter] + Correction;
			ContinuousCP_Acc += dX;
			Counter += 1;
		}

	OuterCounter += 1;
	CPFisher += dX_Rough

	}
}


void LoadCorrelationArray(double * CorrelationArray, double * CPVals) {

	string FilenameBase = "BistableCorrelation_CP_";
	std::ofstream ReadFile;

	double CPFisher = -10.0;

	while(CPFisher <= CPMax) {

		string Filename_Full = FilenameBase + std::to_string(CPFisher) + ".dat";

		/* THE REST OF THIS ROUTINE NEEDS TO READ IN EACH OF THE CORRELATION FUNCTIONS TO EACH ARRAY ROW (ONLY FOR THE CPVALS CHOSEN) */

	}
}



/* Optimization Algorithms */

void OptimizeSwitchTimes(double * SwitchTimes, double * CorrelationArray, double ProtocolTime, int NumTimes) {

	//Choose initial time allocation to be uniform

	double BufferTime = 100.0;
	SwitchTimes[0] = BufferTime;
	SwitchTimes[NumTimes - 1] = BufferTime;

	double UniformTime = ProtocolTime/double(NumTimes - 2);

	for(int k = 1 ; k < NumTimes - 1 ; k++) {
		SwitchTimes[k] = UniformTime;
	}

	double Cost;
	double OldCost; 				
	int Condition = 0;
	OldCost = GradientDescentTime(SwitchTimes,CorrelationArray, ProtocolTime);

	while(Condition==0){
		Cost = GradientDescentTime(SwitchTimes);
		if (Cost - OldCost < Threshold){
			break;
		}
	}
}

double GradientDescentTime(double * SwitchTimes, double * CorrelationArray, double ProtocolTime) {

//This routine will perform a single step of gradient descent on the correlation times and preserve the fixed-time constraint

}

//double GradientDescentFull(double * SwitchTimes, double * CorrelationArray, double * CPVals)








