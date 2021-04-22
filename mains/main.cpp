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
	
	const double R1 = 0.04125;			// inner radius
	const double R2 = 0.0463701;		// outer radius
	const double R = (R2-R1)/2.0+R1;	// radius to represent an object without radial extension
	const double R_TAVP = R2*1.015355; 	// this radius will make B_z match B_z of the helix BS model (R_TAVP = R2_helix*1.015355)
	const double L = 0.0346811; 		// Length
	const double i = 600; 				// current per loop/current in the wire
	const int N_z = 30; 				// number of wire turns/loops per layer/shell
	const int N_rho = 4; 				// number of layers/shell
	const int N_wires = N_z*N_rho;		// total number of wires
	const double I = N_z*N_rho*i;		// "total" amount of current
	const double x = 0.;				// x coordinate of magnet centre
	const double y = 0.;				// y coordinate of magnet centre
	const double z = 0.;				// z coordinate of magnet centre	
	
	int McDOrder;	// number of terms to use in the McDonald model
	int N_BS; 		// number of segments to be used in the Biot-Savart model
	int NG_z; 		// number of loops used to represent a layer in the GQ methods
	int NG_rho;		// number of layers used to represent a magnet in the GQ methods
	double lambda;	// parameter for TAVP model
	
	const int N_t = 10;	// number of test of each method
	double time;
	double time_squared;
	double mean;
	double stdev;
	
	const std::string path = "";
	
	std::cout << "\n";
	std::cout << "Calculating field at point:\n";
	printVec(carP,"P");
	//~ printVec(cylP,"P");
	//~ printVec(sphP,"P");
	std::cout << "\n";
	
	//////////////////// LOOP ////////////////////
	std::cout << "CALCULATING MODELS FOR A CURRENT LOOP\n";
	
	McDOrder = 7;	// number of terms to use in the McDonald model
	N_BS = 1000; 	// number of segments to be used in the Biot-Savart model
		
	std::cout << "Using the (exact) Simple Analytic Model (SAM):\n";
	time = 0;
	time_squared = 0;
	//~ double[N_t][3] cylPArr;
	//~ double[N_t][3] BCylVecArr;
	SimpleAnalyticModel SAM = SimpleAnalyticModel(R,I,x,y,z);
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		SAM.getB(cylP,BCylVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
		
		//~ cylPArr[i][0] = cylP[0];
		//~ cylPArr[i][1] = cylP[1];
		//~ cylPArr[i][2] = cylP[2];
		//~ BCylVecArr[i][0] = BCylVec[0];
		//~ BCylVecArr[i][1] = BCylVec[1];
		//~ BCylVecArr[i][2] = BCylVec[2];
	}	
	printVec(BCylVec,"B");
	//~ writeBFieldToFile author_SAM = new writeToFile(path,"SAM_I" + std::to_string(I));
	//~ author-> write(cylPArr,BCylVecArr,N_t);			
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
	
	NG_z = 3;
	McDOrder = 7;	// number of terms to use in the McDonald model
	N_BS = 1000; 	// number of segments to be used in the Biot-Savart model
	
	
	
	std::cout << "Using the (exact) Conway model:\n";
	time = 0;
	time_squared = 0;
	Conway1D conway1D = Conway1D(R,N_wires,i,L,x,y,z);
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
	McD_Shell mcDShell = McD_Shell(McDOrder,R,N_wires,i,L,x,y,z);
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
	NWire_Shell nWire = NWire_Shell(N_z,R,N_wires,i,L,x,y,z);
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
	GaussianQuadratureLoops_Shell GQL_S = GaussianQuadratureLoops_Shell(N_z,NG_z,R,N_wires,i,L,x,y,z);
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

	McDOrder = 5;
	N_BS = 10000;
	NG_rho = 1;
	NG_z = 3;	
	lambda = 0.866;
	
	
	std::cout << "Using the detailed Biot-Savart model:\n";
	time = 0;
	time_squared = 0;
	Helix helix = Helix(N_z,N_rho,N_BS,R1,R2,N_wires,i,L,x,y,z);
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
	TAVP tavp = TAVP(lambda,R_TAVP,I,x,y,z);
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
	McD_Tube mcD_Tube = McD_Tube(McDOrder,R1,R2,N_wires,i,L,x,y,z);
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
    NWire_Tube nWire_Tube = NWire_Tube(N_z,N_rho,NG_z,NG_rho,R1,R2,N_wires,i,L,x,y,z);
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
	GaussianQuadratureLoops_Tube GQL_T = GaussianQuadratureLoops_Tube(N_z,N_rho,NG_z,NG_rho,R1,R2,N_wires,i,L,x,y,z);
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
	GaussianQuadratureShells_Tube GQS_T = GaussianQuadratureShells_Tube(N_rho,NG_rho,R1,R2,N_wires,i,L,x,y,z);
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
