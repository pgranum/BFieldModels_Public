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
	
	const bool rhoOrZ = true; // true is rho false is z
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
			z_min = 1/2*R1;							
			z_max = 1/2*R1;							
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
			rho_min = 1/2*R1;							
			rho_max = 1/2*R1;							
			suffix = "/z/off";
		}
	}
	
	const double z_d = (z_max-z_min)/((double)N_p-1.0);	// length of a stright line segment
	const double rho_d = (rho_max-rho_min)/((double)N_p-1.0);	// length of a stright line segment


	// setting up files for the code to write results to
	writeBFieldToFile author_SAM (path + "/SAM" + suffix,"SAM");
	writeToFile author_SAM_t(path,"SAM_t");
	
	std::cout << path + "/McDonald/Loop" + suffix << std::endl;

	writeBFieldToFile author_McDLoop1(path + "/McDonald/loop" + suffix,"McD1");
	writeToFile author_McDLoop1_t(path,"McD_loop1_t");
	writeBFieldToFile author_McDLoop2(path + "/McDonald/loop" + suffix,"McD2");
	writeToFile author_McDLoop2_t(path,"McD_loop2_t");
	writeBFieldToFile author_McDLoop3(path + "/McDonald/loop" + suffix,"McD3");
	writeToFile author_McDLoop3_t(path,"McD_loop3_t");
	writeBFieldToFile author_McDLoop4(path + "/McDonald/loop" + suffix,"McD4");
	writeToFile author_McDLoop4_t(path,"McD_loop4_t");
	writeBFieldToFile author_McDLoop5(path + "/McDonald/loop" + suffix,"McD5");
	writeToFile author_McDLoop5_t(path,"McD_loop5_t");
	writeBFieldToFile author_McDLoop6(path + "/McDonald/loop" + suffix,"McD6");
	writeToFile author_McDLoop6_t(path,"McD_loop6_t");
	writeBFieldToFile author_McDLoop7(path + "/McDonald/loop" + suffix,"McD7");
	writeToFile author_McDLoop7_t(path,"McD_loop7_t");
	
	writeBFieldToFile author_BSLoop10(path + "/Biot_Savart" + suffix,"Biot10");
	writeToFile author_BSLoop10_t(path,"Biot10_t");
	writeBFieldToFile author_BSLoop100(path + "/Biot_Savart" + suffix,"Biot100");
	writeToFile author_BSLoop100_t(path,"Biot100_t");
	writeBFieldToFile author_BSLoop1000(path + "/Biot_Savart" + suffix,"Biot1000");
	writeToFile author_BSLoop1000_t(path,"Biot1000_t");
	
	writeBFieldToFile author_Conway1D(path + "/Conway" + suffix,"Conway");
	writeToFile author_Conway1D_t(path,"Conway");
	
	writeBFieldToFile author_McDShell1(path + "/McDonald/shell" + suffix,"McD1");
	writeToFile author_McDShell1_t(path,"McD_shell1_t");
	writeBFieldToFile author_McDShell2(path + "/McDonald/shell" + suffix,"McD2");
	writeToFile author_McDShell2_t(path,"McD_shell2_t");
	writeBFieldToFile author_McDShell3(path + "/McDonald/shell" + suffix,"McD3");
	writeToFile author_McDShell3_t(path,"McD_shell3_t");
	writeBFieldToFile author_McDShell4(path + "/McDonald/shell" + suffix,"McD4");
	writeToFile author_McDShell4_t(path,"McD_shell4_t");
	writeBFieldToFile author_McDShell5(path + "/McDonald/shell" + suffix,"McD5");
	writeToFile author_McDShell5_t(path,"McD_shell5_t");
	writeBFieldToFile author_McDShell6(path + "/McDonald/shell" + suffix,"McD6");
	writeToFile author_McDShell6_t(path,"McD_shell6_t");
	writeBFieldToFile author_McDShell7(path + "/McDonald/shell" + suffix,"McD7");
	writeToFile author_McDShell7_t(path,"McD_shell7_t");
	
	writeBFieldToFile author_NWireShell(path + "/NWire/shell" + suffix,"NWire");
	writeToFile author_NWireShell_t(path + "/NWire/shell" + suffix,"NWire_t");
	
	writeBFieldToFile author_GQLoopsShell3(path + "/GQ/loopshell" + suffix,"GQloopshell3");
	writeToFile author_GQLoopsShell3_t(path ,"GQloopshell3_t");
	writeBFieldToFile author_GQLoopsShell4(path + "/GQ/loopshell" + suffix,"GQloopshell4");
	writeToFile author_GQLoopsShell4_t(path,"GQloopshell4_t");
	writeBFieldToFile author_GQLoopsShell5(path + "/GQ/loopshell" + suffix,"GQloopshell5");
	writeToFile author_GQLoopsShell5_t(path,"GQloopshell5_t");
	
	writeBFieldToFile author_Helix(path + "/Helix" + suffix,"Helix");
	writeToFile author_Helix_t(path,"Helix_t");

	writeBFieldToFile author_TAVP866(path + "/TAVP" + suffix,"TAVP866");
	writeToFile author_TAVP866_t(path,"TAVP866_t");
	writeBFieldToFile author_TAVP902(path + "/TAVP" + suffix,"TAVP902");
	writeToFile author_TAVP902_t(path,"TAVP902_t");
	
	writeBFieldToFile author_McDTube1(path + "/McDonald/tube" + suffix,"McD1");
	writeToFile author_McDTube1_t(path,"McD_tube1_t");
	writeBFieldToFile author_McDTube2(path + "/McDonald/tube" + suffix,"McD2");
	writeToFile author_McDTube2_t(path,"McD_tube2_t");
	writeBFieldToFile author_McDTube3(path + "/McDonald/tube" + suffix,"McD3");
	writeToFile author_McDTube3_t(path,"McD_tube3_t");
	writeBFieldToFile author_McDTube4(path + "/McDonald/tube" + suffix,"McD4");
	writeToFile author_McDTube4_t(path,"McD_tube4_t");
	writeBFieldToFile author_McDTube5(path + "/McDonald/tube" + suffix,"McD5");
	writeToFile author_McDTube5_t(path,"McD_tube5_t");
	
	writeBFieldToFile author_NWireTube(path + "/NWire/tube"+ suffix,"NWire");
	writeToFile author_NWireTube_t(path,"NWire_t");
	
	writeBFieldToFile author_GQLoopsTube13(path + "/GQ/looptube" + suffix,"GQlooptube13");
	writeToFile author_GQLoopsTube13_t(path,"GQlooptube13_t");
	writeBFieldToFile author_GQLoopsTube14(path + "/GQ/looptube" + suffix,"GQlooptube14");
	writeToFile author_GQLoopsTube14_t(path,"GQlooptube14_t");
	writeBFieldToFile author_GQLoopsTube15(path + "/GQ/looptube" + suffix,"GQlooptube15");
	writeToFile author_GQLoopsTube15_t(path,"GQlooptube15_t");
	writeBFieldToFile author_GQLoopsTube23(path + "/GQ/looptube" + suffix,"GQlooptube23");
	writeToFile author_GQLoopsTube23_t(path,"GQlooptube23_t");
	
	writeBFieldToFile author_GQShellsTube1(path + "/GQ/shelltube" + suffix,"GQshelltube1");
	writeToFile author_GQShellsTube1_t(path,"GQshelltube1_t");
	writeBFieldToFile author_GQShellsTube2(path + "/GQ/shelltube" + suffix,"GQshelltube2");
	writeToFile author_GQShellsTube2_t(path,"GQshelltube2_t");
	writeBFieldToFile author_GQShellsTube3(path + "/GQ/shelltube" + suffix,"GQshelltube3");
	writeToFile author_GQShellsTube3_t(path,"GQshelltube3_t");
	writeBFieldToFile author_GQShellsTube4(path + "/GQ/shelltube" + suffix,"GQshelltube4");
	writeToFile author_GQShellsTube4_t(path,"GQshelltube4_t");
	
	for(int n_t=0; n_t<N_t2; n_t++){
		std::cout << "N_t = " << n_t << std::endl;
	for(int n_p=0; n_p<N_p; n_p++){
		
		double carP[3] = {rho_min + n_p*rho_d,0.,z_min + n_p*z_d}; 	// point to calculate field at in cartesian coordinates
		double cylP[3];								// point to calculate field at in cylindrical coordinates
		double sphP[3];								// point to calculate field at in spherical coordinates
		carPToCylP(carP,cylP);						// converting the point in cartesian coor to cylindrical coordinates
		carPToSphP(carP,sphP);						// converting the point in cartesian coor to spherical coordiantes
		double BCarVec[3] = {0.,0.,0.}; 			// placeholder for the field in cartesian coordinates
		double BCylVec[3] = {0.,0.,0.}; 			// placeholder for the field in cylindrical coordinates
		double BSphVec[3] = {0.,0.,0.}; 			// placeholder for the field in spherical coordinates		
	
		std::cout << "\n";
		std::cout << "Calculating field at point:\n";
		printVec(carP,"P_car");
		printVec(cylP,"P_cyl");
		//~ printVec(sphP,"P");
		std::cout << "\n";
		
		//////////////////// LOOP ////////////////////
		std::cout << "CALCULATING MODELS FOR A CURRENT LOOP\n";
		
		McDOrder = 7;	// number of terms to use in the McDonald model
		//~ N_BS = 1000; 	// number of segments to be used in the Biot-Savart model I HAVE MOVED THIS FURTHER DOWN
			
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
		if(n_t == 0){author_SAM.write(cylP,BCylVec);}
		author_SAM_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		
		
		for(int n_McD=1; n_McD<McDOrder+1; n_McD++){
			std::cout << "Using the McDonald model:\n";
			time = 0;
			time_squared = 0;
			McD_Loop mcD_Loop = McD_Loop(n_McD,R,I,x,y,z);
			for(int i=0; i<N_t; i++){
				auto start = std::chrono::steady_clock::now();
				mcD_Loop.getB(cylP,BCylVec);
				auto end = std::chrono::steady_clock::now();
				double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
				time += t;
				time_squared += t*t;
			}	
			printVec(BCylVec,"B");
			if(n_McD == 1){
				if(n_t == 0){author_McDLoop1.write(cylP,BCylVec);}
				author_McDLoop1_t.write(time);
			}else if(n_McD == 2){
				if(n_t == 0){author_McDLoop2.write(cylP,BCylVec);}
				author_McDLoop2_t.write(time);			
			}else if(n_McD == 3){
				if(n_t == 0){author_McDLoop3.write(cylP,BCylVec);}
				author_McDLoop3_t.write(time);
			}else if(n_McD == 4){
				if(n_t == 0){author_McDLoop4.write(cylP,BCylVec);}
				author_McDLoop4_t.write(time);	
			}else if(n_McD == 5){
				if(n_t == 0){author_McDLoop5.write(cylP,BCylVec);}
				author_McDLoop5_t.write(time);	
			}else if(n_McD == 6){
				if(n_t == 0){author_McDLoop6.write(cylP,BCylVec);}
				author_McDLoop6_t.write(time);	
			}else if(n_McD == 7){
				if(n_t == 0){author_McDLoop7.write(cylP,BCylVec);}
				author_McDLoop7_t.write(time);				
			}else{
				std::cout << "You are asking for too high an order on the McDLoop\n";
			}
			
			mean = time/(double)N_t;
			stdev = sqrt( time_squared / (double)N_t - mean * mean );
			//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		}
	
	
	
		N_BS = 10; 	// number of segments to be used in the Biot-Savart model
		std::cout << "Using the Biot-Savart model:\n";
		time = 0;
		time_squared = 0;
		BiotSavart_Loop BS_L10 = BiotSavart_Loop(N_BS,R,I,x,y,z);
		for(int i=0; i<N_t; i++){
			auto start = std::chrono::steady_clock::now();
			BS_L10.getB(carP,BCarVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}
		carVecToCylVec(BCarVec,carP,BCylVec);	
		printVec(BCylVec,"B");
		if(n_t == 0){author_BSLoop10.write(cylP,BCylVec);}
		author_BSLoop10_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		
		
		
		N_BS = 100; 	// number of segments to be used in the Biot-Savart model
		std::cout << "Using the Biot-Savart model:\n";
		time = 0;
		time_squared = 0;
		BiotSavart_Loop BS_L100 = BiotSavart_Loop(N_BS,R,I,x,y,z);
		for(int i=0; i<N_t; i++){
			auto start = std::chrono::steady_clock::now();
			BS_L100.getB(carP,BCarVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}
		carVecToCylVec(BCarVec,carP,BCylVec);	
		printVec(BCylVec,"B");
		if(n_t == 0){author_BSLoop100.write(cylP,BCylVec);}
		author_BSLoop100_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		
		
		
		N_BS = 1000; 	// number of segments to be used in the Biot-Savart model
		std::cout << "Using the Biot-Savart model:\n";
		time = 0;
		time_squared = 0;
		BiotSavart_Loop BS_L1000 = BiotSavart_Loop(N_BS,R,I,x,y,z);
		for(int i=0; i<N_t; i++){
			auto start = std::chrono::steady_clock::now();
			BS_L1000.getB(carP,BCarVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}
		carVecToCylVec(BCarVec,carP,BCylVec);	
		printVec(BCylVec,"B");
		if(n_t == 0){author_BSLoop1000.write(cylP,BCylVec);}
		author_BSLoop1000_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		
		std::cout << "\n";
		
		//////////////////// SHELL ////////////////////
		std::cout << "CALCULATING MODELS FOR A SHELL\n";
		
		//~ NG_z = 3;	// MOVED FURTHER DOWN
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
		if(n_t == 0){author_Conway1D.write(cylP,BCylVec);}
		author_Conway1D_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
	
	
	
		for(int n_McD=1; n_McD<McDOrder+1; n_McD++){
			std::cout << "Using the McDonald model:\n";
			time = 0;
			time_squared = 0;
			McD_Shell mcDShell = McD_Shell(n_McD,R,N_wires,i,L,x,y,z);
			for(int i=0; i<N_t; i++){
				auto start = std::chrono::steady_clock::now();
				mcDShell.getB(cylP,BCylVec);
				auto end = std::chrono::steady_clock::now();
				double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
				time += t;
				time_squared += t*t;
			}	
			printVec(BCylVec,"B");
			if(n_McD == 1){
				if(n_t == 0){author_McDShell1.write(cylP,BCylVec);}
				author_McDShell1_t.write(time);
			}else if(n_McD == 2){
				if(n_t == 0){author_McDShell2.write(cylP,BCylVec);}
				author_McDShell2_t.write(time);			
			}else if(n_McD == 3){
				if(n_t == 0){author_McDShell3.write(cylP,BCylVec);}
				author_McDShell3_t.write(time);
			}else if(n_McD == 4){
				if(n_t == 0){author_McDShell4.write(cylP,BCylVec);}
				author_McDShell4_t.write(time);	
			}else if(n_McD == 5){
				if(n_t == 0){author_McDShell5.write(cylP,BCylVec);}
				author_McDShell5_t.write(time);	
			}else if(n_McD == 6){
				if(n_t == 0){author_McDShell6.write(cylP,BCylVec);}
				author_McDShell6_t.write(time);	
			}else if(n_McD == 7){
				if(n_t == 0){author_McDShell7.write(cylP,BCylVec);}
				author_McDShell7_t.write(time);				
			}else{
				std::cout << "You are asking for too high an order on the McDShell\n";
			}
			mean = time/(double)N_t;
			stdev = sqrt( time_squared / (double)N_t - mean * mean );
			//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		}
		
		
		
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
		if(n_t == 0){author_NWireShell.write(cylP,BCylVec);}
		author_NWireShell_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		
		
		NG_z = 3;
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
		if(n_t == 0){author_GQLoopsShell3.write(cylP,BCylVec);}
		author_GQLoopsShell3_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";	
		
		
		
		NG_z = 4;
		std::cout << "Using the Gaussian Quadrature model:\n";
		time = 0;
		time_squared = 0;
		GaussianQuadratureLoops_Shell GQL_S4 = GaussianQuadratureLoops_Shell(N_z,NG_z,R,N_wires,i,L,x,y,z);
		for(int i=0; i<N_t; i++){
			auto start = std::chrono::steady_clock::now();
			GQL_S4.getB(cylP,BCylVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}	
		printVec(BCylVec,"B");
		if(n_t == 0){author_GQLoopsShell4.write(cylP,BCylVec);}
		author_GQLoopsShell4_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";	
		
		
		
		NG_z = 5;
		std::cout << "Using the Gaussian Quadrature model:\n";
		time = 0;
		time_squared = 0;
		GaussianQuadratureLoops_Shell GQL_S5 = GaussianQuadratureLoops_Shell(N_z,NG_z,R,N_wires,i,L,x,y,z);
		for(int i=0; i<N_t; i++){
			auto start = std::chrono::steady_clock::now();
			GQL_S5.getB(cylP,BCylVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}	
		printVec(BCylVec,"B");
		if(n_t == 0){author_GQLoopsShell5.write(cylP,BCylVec);}
		author_GQLoopsShell5_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";	
		
			
		std::cout << "\n";	
		
		
		//////////////////// FINITE SOLENOID ////////////////////
		std::cout << "CALCULATING MODELS FOR A FINITE SOLENOID\n";
	
		McDOrder = 4;
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
		//~ if(n_t == 0){author_Helix.write(cylP,BCylVec);}
		//~ author_Helix_t.write(time);
		//~ mean = time/(double)N_t;
		//~ stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//~ std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
		//~ std::cout << "\n";
		
		
		lambda = 0.866;
		std::cout << "Using the TAVP model:\n";
		time = 0;
		time_squared = 0;
		TAVP tavp866 = TAVP(lambda,R_TAVP,I,x,y,z);
		for(int i=0; i<N_t; i++){
			//std::cout << "i = " << i << "\n";
			auto start = std::chrono::steady_clock::now();
			tavp866.getB(sphP,BSphVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}	
		sphVecToCarVec(BSphVec,sphP,BCarVec);
		carVecToCylVec(BCarVec,carP,BCylVec);	
		printVec(BCylVec,"B");
		if(n_t == 0){author_TAVP866.write(cylP,BCylVec);}
		author_TAVP866_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//~ std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
		//~ std::cout << "\n";
		
		
		
		lambda = 0.902;
		std::cout << "Using the TAVP model:\n";
		time = 0;
		time_squared = 0;
		TAVP tavp902 = TAVP(lambda,R_TAVP,I,x,y,z);
		for(int i=0; i<N_t; i++){
			//~ std::cout << "i = " << i << "\n";
			auto start = std::chrono::steady_clock::now();
			tavp902.getB(sphP,BSphVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}	
		sphVecToCarVec(BSphVec,sphP,BCarVec);
		carVecToCylVec(BCarVec,carP,BCylVec);	
		printVec(BCylVec,"B");
		if(n_t == 0){author_TAVP902.write(cylP,BCylVec);}
		author_TAVP902_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		
		
		
		for(int n_McD=1; n_McD<McDOrder+1; n_McD++){
			std::cout << "Using the McDonald model:\n";
			time = 0;
			time_squared = 0;
			McD_Tube mcD_Tube = McD_Tube(n_McD,R1,R2,N_wires,i,L,x,y,z);
			for(int i=0; i<N_t; i++){
				auto start = std::chrono::steady_clock::now();
				mcD_Tube.getB(cylP,BCylVec);
				auto end = std::chrono::steady_clock::now();
				double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
				time += t;
				time_squared += t*t;
			}	
			printVec(BCylVec,"B");
			if(n_McD == 1){
				if(n_t == 0){author_McDTube1.write(cylP,BCylVec);}
				author_McDTube1_t.write(time);
			}else if(n_McD == 2){
				if(n_t == 0){author_McDTube2.write(cylP,BCylVec);}
				author_McDTube2_t.write(time);			
			}else if(n_McD == 3){
				if(n_t == 0){author_McDTube3.write(cylP,BCylVec);}
				author_McDTube3_t.write(time);
			}else if(n_McD == 4){
				if(n_t == 0){author_McDTube4.write(cylP,BCylVec);}
				author_McDTube4_t.write(time);	
			}else if(n_McD == 5){
				if(n_t == 0){author_McDTube5.write(cylP,BCylVec);}
				author_McDTube5_t.write(time);			
			}else{
				std::cout << "You are asking for too high an order on the McDTube\n";
			}
			mean = time/(double)N_t;
			stdev = sqrt( time_squared / (double)N_t - mean * mean );
			//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		}
		
		
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
		if(n_t == 0){author_NWireTube.write(cylP,BCylVec);}
		author_NWireTube_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		
		
		NG_rho = 1;	
		NG_z = 3;	
		std::cout << "Using the Gaussian Quadrature model with Loops:\n";
		time = 0;
		time_squared = 0;
		GaussianQuadratureLoops_Tube GQL_T13 = GaussianQuadratureLoops_Tube(N_z,N_rho,NG_z,NG_rho,R1,R2,N_wires,i,L,x,y,z);
		for(int i=0; i<N_t; i++){
			auto start = std::chrono::steady_clock::now();
			GQL_T13.getB(cylP,BCylVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}	
		printVec(BCylVec,"B");
		if(n_t == 0){author_GQLoopsTube13.write(cylP,BCylVec);}
		author_GQLoopsTube13_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		
		
		
		NG_rho = 1;	
		NG_z = 4;	
		std::cout << "Using the Gaussian Quadrature model with Loops:\n";
		time = 0;
		time_squared = 0;
		GaussianQuadratureLoops_Tube GQL_T14 = GaussianQuadratureLoops_Tube(N_z,N_rho,NG_z,NG_rho,R1,R2,N_wires,i,L,x,y,z);
		for(int i=0; i<N_t; i++){
			auto start = std::chrono::steady_clock::now();
			GQL_T14.getB(cylP,BCylVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}	
		printVec(BCylVec,"B");
		if(n_t == 0){author_GQLoopsTube14.write(cylP,BCylVec);}
		author_GQLoopsTube14_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		
		
		
		NG_rho = 1;	
		NG_z = 5;	
		std::cout << "Using the Gaussian Quadrature model with Loops:\n";
		time = 0;
		time_squared = 0;
		GaussianQuadratureLoops_Tube GQL_T15 = GaussianQuadratureLoops_Tube(N_z,N_rho,NG_z,NG_rho,R1,R2,N_wires,i,L,x,y,z);
		for(int i=0; i<N_t; i++){
			auto start = std::chrono::steady_clock::now();
			GQL_T15.getB(cylP,BCylVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}	
		printVec(BCylVec,"B");
		if(n_t == 0){author_GQLoopsTube15.write(cylP,BCylVec);}
		author_GQLoopsTube15_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		
		NG_rho = 2;	
		NG_z = 3;	
		std::cout << "Using the Gaussian Quadrature model with Loops:\n";
		time = 0;
		time_squared = 0;
		GaussianQuadratureLoops_Tube GQL_T23 = GaussianQuadratureLoops_Tube(N_z,N_rho,NG_z,NG_rho,R1,R2,N_wires,i,L,x,y,z);
		for(int i=0; i<N_t; i++){
			auto start = std::chrono::steady_clock::now();
			GQL_T23.getB(cylP,BCylVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}	
		printVec(BCylVec,"B");
		if(n_t == 0){author_GQLoopsTube23.write(cylP,BCylVec);}
		author_GQLoopsTube23_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";

		
		NG_rho = 1;
		std::cout << "Using the Gaussian Quadrature model with Shells:\n";
		time = 0;
		time_squared = 0;
		GaussianQuadratureShells_Tube GQS_T1 = GaussianQuadratureShells_Tube(N_rho,NG_rho,R1,R2,N_wires,i,L,x,y,z);
		for(int i=0; i<N_t; i++){
			auto start = std::chrono::steady_clock::now();
			GQS_T1.getB(cylP,BCylVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}	
		printVec(BCylVec,"B");
		if(n_t == 0){author_GQShellsTube1.write(cylP,BCylVec);}
		author_GQShellsTube1_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		
		
		
		NG_rho = 2;
		std::cout << "Using the Gaussian Quadrature model with Shells:\n";
		time = 0;
		time_squared = 0;
		GaussianQuadratureShells_Tube GQS_T2 = GaussianQuadratureShells_Tube(N_rho,NG_rho,R1,R2,N_wires,i,L,x,y,z);
		for(int i=0; i<N_t; i++){
			auto start = std::chrono::steady_clock::now();
			GQS_T2.getB(cylP,BCylVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}	
		printVec(BCylVec,"B");
		if(n_t == 0){author_GQShellsTube2.write(cylP,BCylVec);}
		author_GQShellsTube2_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		
		
		
		NG_rho = 3;
		std::cout << "Using the Gaussian Quadrature model with Shells:\n";
		time = 0;
		time_squared = 0;
		GaussianQuadratureShells_Tube GQS_T3 = GaussianQuadratureShells_Tube(N_rho,NG_rho,R1,R2,N_wires,i,L,x,y,z);
		for(int i=0; i<N_t; i++){
			auto start = std::chrono::steady_clock::now();
			GQS_T3.getB(cylP,BCylVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}	
		printVec(BCylVec,"B");
		if(n_t == 0){author_GQShellsTube3.write(cylP,BCylVec);}
		author_GQShellsTube3_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		
		
		
		NG_rho = 4;
		std::cout << "Using the Gaussian Quadrature model with Shells:\n";
		time = 0;
		time_squared = 0;
		GaussianQuadratureShells_Tube GQS_T4 = GaussianQuadratureShells_Tube(N_rho,NG_rho,R1,R2,N_wires,i,L,x,y,z);
		for(int i=0; i<N_t; i++){
			auto start = std::chrono::steady_clock::now();
			GQS_T4.getB(cylP,BCylVec);
			auto end = std::chrono::steady_clock::now();
			double t = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
			time += t;
			time_squared += t*t;
		}	
		printVec(BCylVec,"B");
		if(n_t == 0){author_GQShellsTube4.write(cylP,BCylVec);}
		author_GQShellsTube4_t.write(time);
		mean = time/(double)N_t;
		stdev = sqrt( time_squared / (double)N_t - mean * mean );
		//std::cout << "Average calc time = " << mean << " +/- " <<  stdev <<" s\n";
			//std::cout << "\n";
		
		
		
		std::cout << "\n";
		
	}
	}
}
