#include <iostream>

#include "BFieldModels.h"
#include "Loop.h"
#include "mcDTubeSupportClass.h"
#include "Shell.h"
#include "Tube.h"
#include "utils.h"

int main(){
	std::cout.precision(9);	// sets the number of significant digits of cout
	
	double carP[3] = {0.,0.,0.}; 		// point to calculate field at in cartesian coordinates
	double cylP[3];						// point to calculate field at in cylindrical coordinates
	carPToCylP(carP,cylP);				// converting the point in cartesian coor to cylindrical coor
	double BCarVec[3] = {0.,0.,0.}; 	// placeholder for the field in cartesian coordinates
	double BCylVec[3] = {0.,0.,0.}; 	// placeholder for the field in cylindrical coordinates
	
	const int N_t = 1000;
	double time;
	double time_squared;
	double mean;
	double stdev;
	
	std::cout << "\n";
	std::cout << "Calculating field at point:\n";
	printVec(carP,"P");
	std::cout << "\n";
	
	//////////////////// LOOP ////////////////////
	std::cout << "CALCULATING MODELS FOR A CURRENT LOOP\n";
	
	const Loop loop = Loop();		// a loop class that is instantiated with the relevant parameters
	const int McDOrder = 7;	// number of terms to use in the McDonald model
	const int N_BS = 1000; 	// number of segments to be used in the Biot-Savart model
	
	
	
	std::cout << "Using the (exact) Simple Analytic Model (SAM):\n";
	time = 0;
	time_squared = 0;
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		loopExactSAM(loop,cylP,BCylVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}	
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";
	
	std::cout << "Using the McDonald model:\n";
	time = 0;
	time_squared = 0;
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		mcDonald(loop,McDOrder,carP,BCarVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}	
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";



	std::cout << "Using the Biot-Savart model:\n";
	time = 0;
	time_squared = 0;
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		loopBiotSavart(loop,N_BS,carP,BCarVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}	
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";
	
	std::cout << "\n";
	
	//////////////////// SHELL ////////////////////
	std::cout << "CALCULATING MODELS FOR A SHELL\n";
	
	
	const Shell shell = Shell();
	//~ McDOrder = 7; // I would like the option to change here
	const int N_z = 30;
	const int NG_z = 3;
	
	std::cout << "Using the (exact) Conway model:\n";
	time = 0;
	time_squared = 0;
    assert(shell.getx() == 0 && shell.gety() == 0); // the solenoid has to be centered around the axis
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		Conway1D(shell,cylP,BCylVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}	
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";
	
	
	
	std::cout << "Using the McDonald model:\n";
	time = 0;
	time_squared = 0;
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		mcDonald(shell,McDOrder,cylP,BCylVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}	
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";
	
	
	
	std::cout << "Using the N-Wire model:\n";
	time = 0;
	time_squared = 0;
	for(int i=0; i<N_t; i++){
		BCylVec[0] = 0.; BCylVec[1] = 0.; BCylVec[2] = 0.;
		auto start = std::chrono::steady_clock::now();
		NWire(shell,N_z,cylP,BCylVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}	
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";
	
	
	
	std::cout << "Using the Gaussian Quadrature model:\n";
	time = 0;
	time_squared = 0;
	for(int i=0; i<N_t; i++){
		BCylVec[0] = 0.; BCylVec[1] = 0.; BCylVec[2] = 0.;
		auto start = std::chrono::steady_clock::now();
		GaussianQuadratureLoop(shell,N_z,NG_z,cylP,BCylVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}	
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";
	
	
		
	std::cout << "\n";
	
	//////////////////// FINITE SOLENOID ////////////////////
	std::cout << "CALCULATING MODELS FOR A FINITE SOLENOID\n";
	
	const Tube tube = Tube();
	const int McDOrder2 = 3; // I would like the option to change here. PUT IN A WARNING FOR HIGH ORDERS
	const int N_rho = 4;
	//~ const int N_z = 30;
	//~ const int N_BS = 10000; // I would like the option to change here
	const int NG_rho = 1;
	//~ const int NG_z = 3;
	McD_Tube_Support McDSup = McD_Tube_Support(McDOrder2);
	
	
	
	std::cout << "Using the detailed Biot-Savart model:\n";
	time = 0;
	time_squared = 0;
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		Helix(tube,N_rho,N_z,N_BS,carP,BCarVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}	
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";
	
	
	
	std::cout << "Using the McDonald model:\n";
	time = 0;
	time_squared = 0;
	
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		mcDonald(tube,cylP,BCylVec,McDSup);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}	
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";
	
	
	
	std::cout << "Using the N-Wire model:\n";
	time = 0;
    time_squared = 0;
	for(int i=0; i<N_t; i++){
		BCylVec[0] = 0.; BCylVec[1] = 0.; BCylVec[2] = 0.;
		auto start = std::chrono::steady_clock::now();
		NWire(tube,N_rho,N_z,cylP,BCylVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}	
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";
	
	
	
	std::cout << "Using the Gaussian Quadrature model with Loops:\n";
	time = 0;
	time_squared = 0;
	for(int i=0; i<N_t; i++){
		BCylVec[0] = 0.; BCylVec[1] = 0.; BCylVec[2] = 0.;
		auto start = std::chrono::steady_clock::now();
		GaussianQuadratureLoop(tube,N_rho,N_z,NG_rho,NG_z,cylP,BCylVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}	
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";
	
	
	
	std::cout << "Using the Gaussian Quadrature model with Shells:\n";
	time = 0;
	time_squared = 0;
	for(int i=0; i<N_t; i++){
		BCylVec[0] = 0.; BCylVec[1] = 0.; BCylVec[2] = 0.;
		auto start = std::chrono::steady_clock::now();
		GaussianQuadratureShell(tube,N_rho,NG_rho,cylP,BCylVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}	
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";
	
	
	
	std::cout << "\n";
}
