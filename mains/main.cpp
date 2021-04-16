#include <iostream>

#include "BFieldModels.h"
#include "Loop.h"
#include "Shell.h"
#include "Tube.h"
#include "utils.h"

int main(){
	std::cout.precision(15);	// sets the number of significant digits of cout
	
	double carP[3] = {0.,0.,0.}; 		// point to calculate field at in cartesian coordinates
	double cylP[3];						// point to calculate field at in cylindrical coordinates
	double sphP[3];						// point to calculate field at in spherical coordinates
	carPToCylP(carP,cylP);				// converting the point in cartesian coor to cylindrical coordinates
	carPToSphP(carP,sphP);				// converting the point in cartesian coor to spherical coordiantes
	double BCarVec[3] = {0.,0.,0.}; 	// placeholder for the field in cartesian coordinates
	double BCylVec[3] = {0.,0.,0.}; 	// placeholder for the field in cylindrical coordinates
	double BSphVec[3] = {0.,0.,0.}; 	// placeholder for the field in spherical coordinates
	
	const int N_t = 1;
	double time;
	double time_squared;
	double mean;
	double stdev;
	
	std::cout << "\n";
	std::cout << "Calculating field at point:\n";
	printVec(carP,"P");
	//~ printVec(cylP,"P");
	//~ printVec(sphP,"P");
	std::cout << "\n";
	
	//////////////////// LOOP ////////////////////
	std::cout << "CALCULATING MODELS FOR A CURRENT LOOP\n";
	const double R = (0.0463701-0.04125)/2.0+0.04125;
	const double I = 120*600;
	const double x = 0;
	const double y = 0;
	const double z = 0.;
	const int McDOrder = 7;	// number of terms to use in the McDonald model
	const int N_BS = 1000; 	// number of segments to be used in the Biot-Savart model
		
	std::cout << "Using the (exact) Simple Analytic Model (SAM):\n";
	time = 0;
	time_squared = 0;
	SimpleAnalyticModel SAM = SimpleAnalyticModel(R,I,x,y,z);
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		SAM.getB(cylP,BCylVec);
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
	McD_Loop mcD_Loop = McD_Loop(McDOrder,R,I,x,y,z);
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		mcD_Loop.getB(cylP,BCylVec);
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
	BiotSavart_Loop BS_L = BiotSavart_Loop(N_BS,R,I,x,y,z);
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		BS_L.getB(carP,BCarVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}
	carVecToCylVec(BCarVec,carP,BCylVec);	
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";
	
	std::cout << "\n";
	
	//////////////////// SHELL ////////////////////
	std::cout << "CALCULATING MODELS FOR A SHELL\n";
	
	//~ const double R = (0.0463701-0.04125)/2.0+0.04125;
	const int N = 120;
	const double i = 600;
	const double L = 0.0346811;
	//~ const double x = 0;
	//~ const double y = 0;
	//~ const double z = 0;
	//~ const Shell shell = Shell(R,N,i,L,x,y,z); // a shell class that is instantiated with the relevant parameters
	//~ McDOrder = 7; // I would like the option to change here
	const int N_z = 30;
	const int NG_z = 3;
	
	std::cout << "Using the (exact) Conway model:\n";
	time = 0;
	time_squared = 0;
	Conway1D conway1D = Conway1D(R,N,i,L,x,y,z);
    assert(x == 0 && y == 0); // the solenoid has to be centered around the axis
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		conway1D.getB(cylP,BCylVec);
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
	McD_Shell mcDShell = McD_Shell(McDOrder,R,N,i,L,x,y,z);
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		mcDShell.getB(cylP,BCylVec);
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
	NWire_Shell nWire = NWire_Shell(N_z,R,N,i,L,x,y,z);
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		nWire.getB(cylP,BCylVec);
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
	GaussianQuadratureLoops_Shell GQL_S = GaussianQuadratureLoops_Shell(N_z,NG_z,R,N,i,L,x,y,z);
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		GQL_S.getB(cylP,BCylVec);
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
	
	const double R1 = 0.04125;
	const double R2 = 0.04637;
	const double L2 = 0.03468;
	//~ const int N = 4*30;
	//~ const double L = 0.03468;
	//~ const double x = 0;
	//~ const double y = 0;
	//~ const double z = 0;
	
	const int McDOrder2 = 3; // I would like the option to change here. PUT IN A WARNING FOR HIGH ORDERS
	const int N_rho = 4;
	//~ const int N_z = 30;
	//~ const int N_BS = 10000; // I would like the option to change here
	const int NG_rho = 1;
	//~ const int NG_z = 3;	
	const double lambda = 0.866;
	
	
	std::cout << "Using the detailed Biot-Savart model:\n";
	time = 0;
	time_squared = 0;
	Helix helix = Helix(N_z,N_rho,N_BS,R1,R2,N,i,L2,x,y,z);
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		helix.getB(carP,BCarVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}
	carVecToCylVec(BCarVec,carP,BCylVec);		
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";
	
	std::cout << "Using the TAVP model:\n";
	time = 0;
	time_squared = 0;
	TAVP tavp = TAVP(lambda,R,I,x,y,z);
	for(int i=0; i<N_t; i++){
		//~ std::cout << "i = " << i << "\n";
		auto start = std::chrono::steady_clock::now();
		tavp.getB(sphP,BSphVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}	
	sphVecToCarVec(BSphVec,sphP,BCarVec);
	carVecToCylVec(BCarVec,carP,BCylVec);	
	printVec(BCylVec,"B");
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
	std::cout << "\n";
	
	std::cout << "Using the McDonald model:\n";
	time = 0;
	time_squared = 0;
	McD_Tube mcD_Tube = McD_Tube(McDOrder2,R1,R2,N,i,L2,x,y,z);
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		mcD_Tube.getB(cylP,BCylVec);
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
    NWire_Tube nWire_Tube = NWire_Tube(N_z,N_rho,NG_z,NG_rho,R1,R2,N,i,L2,x,y,z);
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		nWire_Tube.getB(cylP,BCylVec);
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
	GaussianQuadratureLoops_Tube GQL_T = GaussianQuadratureLoops_Tube(N_z,N_rho,NG_z,NG_rho,R1,R2,N,i,L2,x,y,z);
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		GQL_T.getB(cylP,BCylVec);
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
	GaussianQuadratureShells_Tube GQS_T = GaussianQuadratureShells_Tube(N_rho,NG_rho,R1,R2,N,i,L2,x,y,z);
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		GQS_T.getB(cylP,BCylVec);
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
