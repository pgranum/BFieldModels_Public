#include <iostream>

#include "BFieldModels.h"
#include "Loop.h"
#include "Shell.h"
#include "Tube.h"
#include "utils.h"

int main(){
	std::cout.precision(15);	// sets the number of significant digits of cout
	
	
	
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
	const double z = 0;					// z coordinate of magnet centre	
	
	int McDOrder;	// number of terms to use in the McDonald model
	int N_BS; 		// number of segments to be used in the Biot-Savart model
	int NG_z; 		// number of loops used to represent a layer in the GQ methods
	int NG_rho;		// number of layers used to represent a magnet in the GQ methods
	double lambda;	// parameter for TAVP model
	
	// setting the code up to run multiple times to get statistics of the computation time
	const int N_t = 1;			// number of evaluations of each method within the position loop
	const int N_t2 = 1;			// number of evaluations of each method outside the position loop
	double time;
	double time_squared;
	double mean;
	double stdev;
	
	//~ const std::string path = "/home/magn5452/Data";
	const std::string path = "/home/penielse/BFieldModels_Public/BinFilesNt1Np10000";
	
	// setting the code up to loop over multiple points in space along a straight line

	const int N_p = 100;								// number of points along the line (number of segments = N_p-1 )
	
	const bool rhoOrZ = false; // true is rho false is z
	const bool onOrOff = true; // true is on false is off
	const double z_bound = 3; // The maximum value of R1 of the paths on the axis

	double z_min;
	double z_max;
	double rho_min;
	double rho_max;
	std::string suffix;

	if (rhoOrZ) {
		if (onOrOff) {
			z_min = 0;							
			z_max = 0;					
			rho_min = 0;							
			rho_max = R1;							
			suffix = "/rho/on";
		} else {
			z_min = 1./2.*R1;							
			z_max = 1./2.*R1;							
			rho_min = 0;							
			rho_max = R1;							
			suffix = "/rho/off";
		}
	} else {
		if (onOrOff) {
			z_min = -z_bound*R1;							
			z_max = z_bound*R1;							
			rho_min = 0;							
			rho_max = 0;							
			suffix = "/z/on";
		} else {
			z_min = -z_bound*R1;							
			z_max = z_bound*R1;							
			rho_min = 1./2.*R1;							
			rho_max = 1./2.*R1;							
			suffix = "/z/off";
		}
	}
	
	const double z_d = (z_max-z_min)/((double)N_p-1.0);	// length of a stright line segment
	const double rho_d = (rho_max-rho_min)/((double)N_p-1.0);	// length of a stright line segment

	double carP[3] = {R1*0.9,0,0}; 	// point to calculate field at in cartesian coordinates
	double cylP[3];								// point to calculate field at in cylindrical coordinates
	double sphP[3];								// point to calculate field at in spherical coordinates
	carPToCylP(carP,cylP);						// converting the point in cartesian coor to cylindrical coordinates
	carPToSphP(carP,sphP);						// converting the point in cartesian coor to spherical coordiantes
	double BCarVec[3] = {0.,0.,0.}; 			// placeholder for the field in cartesian coordinates
	double BCylVec[3] = {0.,0.,0.}; 			// placeholder for the field in cylindrical coordinates
	double BSphVec[3] = {0.,0.,0.}; 			// placeholder for the field in spherical coordinates		
	

	//////////////////////////// LOOP ////////////////////////
	
	McDOrder = 7;	// number of terms to use in the McDonald model
	//~ N_BS = 1000; 	// number of segments to be used in the Biot-Savart model I HAVE MOVED THIS FURTHER DOWN
	
	//~ std::cout << "Using the (exact) Simple Analytic Model (SAM):\n";
	//~ time = 0;
	//~ time_squared = 0;
	//~ SimpleAnalyticModel SAM = SimpleAnalyticModel(R,I,x,y,z);		
	//~ for(int i=0; i<N_t; i++){
		//~ auto start = std::chrono::steady_clock::now();
		//~ SAM.getB(cylP,BCylVec);
		//~ auto end = std::chrono::steady_clock::now();
		//~ double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		//~ time += t;
		//~ time_squared += t*t;
	//~ }	
	//~ printVec(BCylVec,"B");
	//~ cylVecToCarVec(BCylVec,cylP,BCarVec);
	//~ mean = time/(double)N_t;
	//~ stdev = sqrt( time_squared / (double)N_t - mean * mean );
	
	//////////////////////////// SHELL ////////////////////////

	McDOrder = 7;	// number of terms to use in the McDonald model
		
				
		
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
	cylVecToCarVec(BCylVec,cylP,BCarVec);
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );
	
	NG_z = 4;
	std::cout << "Using the Gaussian Quadrature model:\n";
	time = 0;
	time_squared = 0;
	GaussianQuadratureLoops_Shell GQL_S3 = GaussianQuadratureLoops_Shell(N_z,NG_z,R,N_wires,i,L,x,y,z);
	for(int i=0; i<N_t; i++){
		auto start = std::chrono::steady_clock::now();
		GQL_S3.getB(cylP,BCylVec);
		auto end = std::chrono::steady_clock::now();
		double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		time += t;
		time_squared += t*t;
	}	
	printVec(BCylVec,"B");
	cylVecToCarVec(BCylVec,cylP,BCarVec);
	mean = time/(double)N_t;
	stdev = sqrt( time_squared / (double)N_t - mean * mean );

	//////////////////////////// FINITE SOLENOID ////////////////////////

	McDOrder = 5;
	N_BS = 10000;
	//~ NG_rho = 1;	//MOVED FURTHER DOWN
	//~ NG_z = 3;	//MOVED FURTHER DOWN
	//~ lambda = 0.866; //MOVED FURTHER DOWN
		
		
	//~ std::cout << "Using the detailed Biot-Savart model:\n";
	//~ time = 0;
	//~ time_squared = 0;
	//~ Helix helix = Helix(N_z,N_rho,N_BS,R1,R2,N_wires,i,L,x,y,z);
	//~ for(int i=0; i<N_t; i++){
		//~ auto start = std::chrono::steady_clock::now();
		//~ helix.getB(carP,BCarVec);
		//~ auto end = std::chrono::steady_clock::now();
		//~ double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		//~ time += t;
		//~ time_squared += t*t;
	//~ }
	//~ carVecToCylVec(BCarVec,carP,BCylVec);		
	//~ printVec(BCylVec,"B");
	//~ mean = time/(double)N_t;
	//~ stdev = sqrt( time_squared / (double)N_t - mean * mean );	
	

	//~ lambda = 0.866;
	//~ std::cout << "Using the TAVP model:\n";
	//~ time = 0;
	//~ time_squared = 0;
	//~ TAVP tavp866 = TAVP(lambda,R_TAVP,I,x,y,z);
	//~ for(int i=0; i<N_t; i++){
		//~ //std::cout << "i = " << i << "\n";
		//~ auto start = std::chrono::steady_clock::now();
		//~ tavp866.getB(sphP,BSphVec);
		//~ auto end = std::chrono::steady_clock::now();
		//~ double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		//~ time += t;
		//~ time_squared += t*t;
	//~ }	
	//~ sphVecToCarVec(BSphVec,sphP,BCarVec);
	//~ carVecToCylVec(BCarVec,carP,BCylVec);
	//~ printVec(BCylVec,"B");
	//~ mean = time/(double)N_t;
	//~ stdev = sqrt( time_squared / (double)N_t - mean * mean );
}
