/* Equilibrium Sampling Driver Implementation written in C++ with a Correlation function calculator built in */

#include <fstream>
#include <iostream>
#include <random>
#include <cmath>
#include <string>
#include <ctime>

using namespace std;



/* Global declaration of random device */

std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<> d(0,1);
double GaussRandom;



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
double BistableTrap = 1.0;



/* Global variables for simulation parameters */

double EqTime = 500;
double MAXLAG = 500;
double TotalStats = 10000000;
//double CP = -10.0;
//double CPMax = 10.0;
double CP = -1.0;
double CPMax = 1.0;
//double dX = 0.25;
double dX = 0.01;
int HUN_MILLION = 100000000;
int TEN_MILLION = 10000000;
int ONE_MILLION = 1000000;



/* Subroutine prototypes */

void Langevin(double * position, double * velocity, double * time);
void LangevinBistable(double * position, double * velocity, double * time);
double ForceParticleTrap(double position);
double ForceParticleBistable(double position);
double BistableWell(double position);
double ForceParticleBistableTrap(double position);

void TrajectoryTracker(double SampleLength);
void TrajectoryCorrelation(double * ForceArray, double CP, int SampleLength2);



/* Main Routine */

int main(){

	//TrajectoryTracker(1000000);
	
	while(CP <= CPMax){

		cout << "CPVal --> " << CP << "\r" << flush;

		string WriteNameCorr = "EquilibriumData/BistableCorrelation_CP_" + std::to_string(CP) + ".dat";

		int ArrayLength = int(MAXLAG/dt) + 1;

		double * Corr;
		Corr = new double [ArrayLength];
		double * LagTime;
		LagTime = new double [ArrayLength];
		int Accumulator = 0;


		double * ForceWindow;
		ForceWindow = new double [ArrayLength];
		double ForceAccumulator = 0;
		double ForceMean;

		double time = 0;
		double position = CP;
		double velocity = 0;

		double * timePointer = &time;
		double * positionPointer = &position;
		double * velocityPointer = &velocity;


		/* Run Langevin integrator to equilibrate sytem at CP */

		while(time <= EqTime) {
			LangevinBistable(positionPointer, velocityPointer, timePointer);
		}

		time = 0;


		/* Initialize the ForceWindow and LagTime arrays */

		LagTime[0] = time;
		ForceWindow[0] = ForceParticleBistableTrap(position);
		ForceAccumulator += ForceWindow[0];

		for(int k = 1 ; k < ArrayLength ; k++) {
			LangevinBistable(positionPointer, velocityPointer, timePointer);
			ForceWindow[k] = ForceParticleBistableTrap(position);
			LagTime[k] = time;
			ForceAccumulator += ForceWindow[k];
		}


		/* Calculate the first entries in the Correlation array */

		for(int k = 0 ; k < ArrayLength ; k++) {
			Corr[k] = ForceWindow[0]*ForceWindow[k];
		}
	

		/* Calculate the Correlation array with TotalStats number of entries at each position */

		for(int k = 0 ; k < TotalStats ; k++) {
			LangevinBistable(positionPointer, velocityPointer, timePointer);

			for(int i = 0 ; i < ArrayLength-1 ; i++) {
				ForceWindow[i] = ForceWindow[i+1];
			}

			ForceWindow[ArrayLength-1] = ForceParticleBistableTrap(position);
			ForceAccumulator += ForceWindow[ArrayLength-1];

			for(int j = 0 ; j < ArrayLength ; j++) {
				Corr[j] += ForceWindow[0]*ForceWindow[j];
			}
		}

		/* Average Correlation array entries and Subtract off the Force Variance */
	
		ForceMean = ForceAccumulator/(TotalStats + 1);

		for(int k = 0 ; k < ArrayLength ; k++){
			Corr[k] = (Corr[k]/double(TotalStats + 1)) - ForceMean*ForceMean;
		}

		/* Write Data to the Output File */

		std::ofstream Writefile;

		Writefile.open(WriteNameCorr);
		Writefile << "LagTime\tCorrelation\n\n";
		for(int k = 0 ; k < ArrayLength ; k++){
			Writefile << LagTime[k] << "\t" << Corr[k] << "\n";
		}
		Writefile.close();

		CP += dX;

		delete Corr;
		delete LagTime;
		delete ForceWindow;

		}

	}



/* Trajectory tracking routine, this prints a long stationary trajectory to make sure the data looks correct */

void TrajectoryTracker(double SampleLength) {

	CP = -10.0;

	int SampleLength2 = int(SampleLength/dt);

	double * positionArray;
	positionArray = new double [SampleLength2];
	double * timeArray;
	timeArray = new double [SampleLength2];
	double * ForceArray;
	ForceArray = new double [SampleLength2];

	double position = CP;
	double velocity = 0;
	double time = 0;

	double * positionPointer;
	double * velocityPointer;
	double * timePointer;

	positionPointer = &position;
	velocityPointer = &velocity;
	timePointer = &time;

	for(int k = 0 ; k < SampleLength2 ; k++) {
		LangevinBistable(positionPointer,velocityPointer,timePointer);
		positionArray[k] = position;
		timeArray[k] = time;
		ForceArray[k] = ForceParticleBistableTrap(position);
	}

	TrajectoryCorrelation(ForceArray, CP, SampleLength2);

	std:ofstream Writefile;

	Writefile.open("EquilibriumData/Trajectory_m10.dat");
	Writefile << "Time\tPosition\n\n";
	for(int k = 0 ; k < SampleLength2 ; k++) {
		Writefile << timeArray[k] << "\t" << positionArray[k] << "\n";
	}
	Writefile.close();

	Writefile.open("EquilibriumData/TrajectoryForce_m10.dat");
	Writefile << "Time\tForce\n\n";
	for(int k = 0 ; k < SampleLength2 ; k++) {
		Writefile << timeArray[k] << "\t" << ForceArray[k] << "\n";
	}
	Writefile.close();


	CP = 0.0;

	position = CP;
	velocity = 0;
	time = 0;

	for(int k = 0 ; k < SampleLength2 ; k++) {
		LangevinBistable(positionPointer,velocityPointer,timePointer);
		positionArray[k] = position;
		timeArray[k] = time;
		ForceArray[k] = ForceParticleBistableTrap(position);
	}

	TrajectoryCorrelation(ForceArray, CP, SampleLength2);

	Writefile.open("EquilibriumData/Trajectory_0.dat");
	Writefile << "Time\tPosition\n\n";
	for(int k = 0 ; k < SampleLength2 ; k++) {
		Writefile << timeArray[k] << "\t" << positionArray[k] << "\n";
	}
	Writefile.close();

	Writefile.open("EquilibriumData/TrajectoryForce_0.dat");
	Writefile << "Time\tForce\n\n";
	for(int k = 0 ; k < SampleLength2 ; k++) {
		Writefile << timeArray[k] << "\t" << ForceArray[k] << "\n";
	}
	Writefile.close();

}


void TrajectoryCorrelation(double * ForceArray, double CP, int SampleTime) {

	int CorrLagMAX = int(500/dt);

	double * LagTime;
	LagTime = new double [CorrLagMAX];

	double * Correlation;
	Correlation = new double [CorrLagMAX];

	double TimeCounter = 0.0;

	for(int k = 0 ; k < CorrLagMAX ; k++) {
		Correlation[k] = 0.0;
		LagTime[k] = TimeCounter;
		TimeCounter += dt;
	}

	for(int k = 0 ; k < SampleTime-CorrLagMAX ; k++) {
		for(int i = 0 ; i < CorrLagMAX ; i++) {
			Correlation[i] += ForceArray[k]*ForceArray[k+i];
		}
	}

	for(int k = 0 ; k < CorrLagMAX ; k++) {
		Correlation[k] = Correlation[k]/double(SampleTime-CorrLagMAX);
	}

	std:ofstream Writefile;

	Writefile.open("EquilibriumData/ForceCorrelation_" + std::to_string(int(CP)) + ".dat");
	Writefile << "LagTime" << "\t" << "Correlation" << "\n\n";
	for(int k = 0 ; k < CorrLagMAX ; k++) {
		Writefile << LagTime[k] << "\t" << Correlation[k] << "\n";
	}
	Writefile.close();

}



/* Langevin Integrators */

void Langevin(double * position, double * velocity, double * time){

	GaussRandom = d(gen);

	*velocity = sqrt(DampingVal)*(*velocity) + sqrt((1-DampingVal)/(beta*mass))*GaussRandom;
	*velocity = *velocity + 0.5*dt*ForceParticleTrap(*position)/mass;
	*position = *position + 0.5*dt*(*velocity);

	*time += dt;

	GaussRandom = d(gen);

	*position = *position + 0.5*dt*(*velocity);
	*velocity = *velocity + 0.5*dt*ForceParticleTrap(*position)/mass;
	*velocity = sqrt(DampingVal)*(*velocity) + sqrt((1-DampingVal)/(beta*mass))*GaussRandom;

}


void LangevinBistable(double * position, double * velocity, double * time){

	GaussRandom = d(gen);

	*velocity = sqrt(DampingVal)*(*velocity) + sqrt((1-DampingVal)/(beta*mass))*GaussRandom;
	*velocity = *velocity + 0.5*dt*ForceParticleBistable(*position)/mass;
	*position = *position + 0.5*dt*(*velocity);

	*time += dt;

	GaussRandom = d(gen);

	*position = *position + 0.5*dt*(*velocity);
	*velocity = *velocity + 0.5*dt*ForceParticleBistable(*position)/mass;
	*velocity = sqrt(DampingVal)*(*velocity) + sqrt((1-DampingVal)/(beta*mass))*GaussRandom;

}



/* Auxillary Numerical Routines */

double ForceParticleTrap(double position){

	double Force = -TrapStrength*(position - CP);
	return Force;
}

double ForceParticleBistable(double position){

	double EnergyLeft;
	double EnergyRight;
	double ForceBistable;

	EnergyLeft = BistableWell(position - dX);
	EnergyRight = BistableWell(position + dX);

	ForceBistable = -1.0*(EnergyRight - EnergyLeft)/(2*dX);

	return ForceBistable;

}

double ForceParticleBistableTrap(double position) {

	double Force = -BistableTrap*(position - CP);
	return Force;
}

double BistableWell(double position){

	double BistableEnergy;
	double CPEnergy;
	double TotalEnergy;

	BistableEnergy = (-1.0/beta)*log(exp(-0.5*beta*kL*((position + X_m)*(position + X_m))) + exp(-0.5*beta*kR*((position - X_m)*(position - X_m)) - beta*DeltaE));
	CPEnergy = 0.5*BistableTrap*(position - CP)*(position - CP);

	TotalEnergy = BistableEnergy + CPEnergy;

	return TotalEnergy;

}




