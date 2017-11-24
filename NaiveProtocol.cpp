/* This C++ code generates a series of Naive Protocol sequences for Discrete Optimal Control Simualtions */

#include <fstream>
#include <iostream>
#include <random>
#include <cmath>
#include <string>
#include <ctime>

using namespace std;



/* Global variables for simulation parameters */

double DeltaCP = 20.0;
double CPMin = -10.0;
double CPMax = 10.0;

double MAXTIME = 100000;
double MAXSTEPS = 1500;


/* Function Prototypes */

void NaiveProtocol(int NumSteps, double ProtocolTime);


/* Main Naive protocol Generator Code */

int main() {
	
	double ProtocolTime = 0.1;
	int NumSteps = 5;

	int Iterations = 20;

	while(ProtocolTime <= MAXTIME) {
	//for(int k = 0 ; k < Iterations ; k++) {

		NumSteps = 5;

		while(NumSteps <= MAXSTEPS) {
		//for(int j = 0 ; j < Iterations ; j++) {

			NaiveProtocol(NumSteps,ProtocolTime);

			NumSteps = NumSteps*2.0;

		}

		ProtocolTime = ProtocolTime*10.0;
	}

}


/* Routine to generate Naive protocol sequences */


void NaiveProtocol(int NumSteps, double ProtocolTime) {

	string Filename = "Protocols/Naive/Naive_Steps_" + std::to_string(NumSteps) + "_Time_" + std::to_string(ProtocolTime) + ".dat";

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









