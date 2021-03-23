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
	
	std::cout << std::endl;
	std::cout << "Calculating field at point:" << std::endl;
	printVec(carP,"P");
	std::cout << std::endl;
	
	//////////////////// LOOP ////////////////////
	std::cout << "CALCULATING MODELS FOR A CURRENT LOOP" << std::endl;
	
	Loop loop = Loop();		// a loop class that is instantiated with the relevant parameters
	const int McDOrder = 7;	// number of terms to use in the McDonald model
	const int N_BS = 1000; 	// number of segments to be used in the Biot-Savart model
	
	std::cout << "Using the (exact) Simple Analytic Model (SAM):" << std::endl;
	loopExactSAM(loop,cylP,BCylVec);
	printVec(BCylVec,"B");
	
	std::cout << "Using the McDonald model:" << std::endl;
	mcDonald(loop,McDOrder,carP,BCarVec);
	printVec(BCarVec,"B");

	std::cout << "Using the Biot-Savart model:" << std::endl;
	loopBiotSavart(loop,N_BS,carP,BCarVec);
	printVec(BCarVec,"B");
	
	std::cout << std::endl;
	
	//////////////////// SHELL ////////////////////
	std::cout << "CALCULATING MODELS FOR A SHELL" << std::endl;
	
	
	Shell shell = Shell();
	//~ McDOrder = 7; // I would like the option to change here
	const int N_z = 30;
	const int NG_z = 3;
	
	std::cout << "Using the (exact) Conway model:" << std::endl;
	Conway1D(shell,cylP,BCylVec);
	printVec(BCylVec,"B");
	
	std::cout << "Using the McDonald model:" << std::endl;
	mcDonald(shell,McDOrder,cylP,BCylVec);
	printVec(BCylVec,"B");
	
	
	BCylVec[0] = 0.; BCylVec[1] = 0.; BCylVec[2] = 0.;
	std::cout << "Using the N-Wire model:" << std::endl;
	NWire(shell,N_z,cylP,BCylVec);
	printVec(BCylVec,"B");
	
	BCylVec[0] = 0.; BCylVec[1] = 0.; BCylVec[2] = 0.;
	std::cout << "Using the Gaussian Quadrature model:" << std::endl;
	GaussianQuadratureLoop(shell,N_z,NG_z,cylP,BCylVec);
	printVec(BCylVec,"B");
		
	std::cout << std::endl;
	
	//////////////////// FINITE SOLENOID ////////////////////
	std::cout << "CALCULATING MODELS FOR A FINITE SOLENOID" << std::endl;
	
	Tube tube = Tube();
	const int McDOrder2 = 3; // I would like the option to change here. PUT IN A WARNING FOR HIGH ORDERS
	const int N_rho = 4;
	//~ const int N_z = 30;
	//~ const int N_BS = 10000; // I would like the option to change here
	const int NG_rho = 1;
	//~ const int NG_z = 3;
	McD_Tube_Support McDSup = McD_Tube_Support(McDOrder2);
	
	std::cout << "Using the detailed Biot-Savart model:" << std::endl;
	Helix(tube,N_rho,N_z,N_BS,carP,BCarVec);
	printVec(BCylVec,"B");
	
	std::cout << "Using the McDonald model:" << std::endl;
	mcDonald(tube,McDOrder2,cylP,BCylVec,McDSup);
	printVec(BCylVec,"B");
	
	BCylVec[0] = 0.; BCylVec[1] = 0.; BCylVec[2] = 0.;
	std::cout << "Using the N-Wire model:" << std::endl;
	NWire(tube,N_rho,N_z,cylP,BCylVec);
	printVec(BCylVec,"B");
	
	BCylVec[0] = 0.; BCylVec[1] = 0.; BCylVec[2] = 0.;
	std::cout << "Using the Gaussian Quadrature model with Loops:" << std::endl;
	GaussianQuadratureLoop(tube,N_rho,N_z,NG_rho,NG_z,cylP,BCylVec);
	printVec(BCylVec,"B");
	
	BCylVec[0] = 0.; BCylVec[1] = 0.; BCylVec[2] = 0.;
	std::cout << "Using the Gaussian Quadrature model with Shells:" << std::endl;
	GaussianQuadratureShell(tube,N_rho,NG_rho,cylP,BCylVec);
	printVec(BCylVec,"B");
	
	std::cout << std::endl;
}
