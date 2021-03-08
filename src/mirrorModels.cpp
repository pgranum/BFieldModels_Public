#include "mirrorModels.h"

//NWire models

void exactMirrorModel(std::string path, double carP1[3], double carP2[3]){
	// this model represents the magnet with an NxM wire grid, and calculates the total field
	// as the sum of the individual wire contributions. The field of a wire is calculated
	// with the SAM or McD method (see what method is commented in)	
	
	double coilCentre[3]{0.,0.,0.};
	
	// Different radii
	const double innerRad = 0.04125;				// 43.0658 mm AG, 41.25 mm A2
	const double outerRad = 0.04637;				// 46.37 mm A2
	const double width_rho = outerRad - innerRad;
	//~ double centreRad = (width_rho)/2.0 + innerRad;
	
	const int N_rho = 4; 							// 8 layers AG, 4 A2
	const int N_z = 30; 							// 90 windings per layer AG, 30 A2
	const int N_wires = N_rho*N_z;
	const double width_z = 0.03468; 				// 34.6811 mm A2
	const double delta_rho = width_rho/N_rho; 	// 533.2 um AG (table value), 1156 um A2 (derived value)
	const double delta_z = width_z/N_z; 			// 385.3 um AG (table value), 1280 um A2 (derived value)
	const double I = 120*600; 					// The total current, given by the number of wires times the current per wire: 710*100 A AG, 120*600 A A2
	
	// define the line to calculate field on:	
	const int N_P = 100; //number of points on line
	double BArr[N_P][3];
	double carPArr[N_P][3];
	
	for(int k = 0; k < N_P; k++){
		BArr[k][0] = 0.0;
		BArr[k][1] = 0.0;
		BArr[k][2] = 0.0;
	}	
	
	// calculate field for a full mirror coil	
	for(int i = 0; i < N_rho; i++){
		for(int j = 0; j < N_z; j++){
			// The radius and axial position of the loop is set, so the wires are distributed evenly over the whole volume
			// The current of the loop is set, so the currents in the loop add up to the total current
			// std::cout << innerRad + delta_rho/2.0 + i*delta_rho << std::endl;
			Loop loop_i = Loop(innerRad + delta_rho/2.0 + i*delta_rho, I/N_wires, 0, 0, coilCentre[2] - 0.5*width_z + delta_z/2.0 + j*delta_z);
			
			double BArr_i[N_P][3];			
			SAMWrapper(loop_i, carP1, carP2, N_P, carPArr, BArr_i);
			//mcDonaldWrapper(loop_i, 7, carP1, carP2,N_P, carPArr, BArr_i);

			for(int k = 0; k < N_P; k++){
				BArr[k][0] += BArr_i[k][0];
				BArr[k][1] += BArr_i[k][1];
				BArr[k][2] += BArr_i[k][2];
			}		
		}
	}
	
	//~ printVec(BArr[0],"B1");
	//~ printVec(BArr[N_P-1],"B2");
	//~ printVec(carPArr[0],"r1");
	//~ printVec(carPArr[N_P-1],"r2");
	
	//writeBFieldToFile author_exact(path,"NWire_1_30");
	//author_exact.write(carPArr, BArr, N_P);
}

void exactMirrorModel_vecClass(std::string path){
	// this model represents the magnet with an NxM wire grid, and calculates the total field
	// as the sum of the individual wire contributions. The field of a wire is calculated
	// with the SAM or McD method (see what method is commented in)	
	CylindricalVector coilCentre{0., 0., 0.};

	double I = 1; // 100 A AG

	// Different radii
	double innerRad = 0.04125;				// 43.0658 mm AG, 41.25 mm A2
	double outerRad = 0.04637;				// 46.37 mm A2
	double width_rho = outerRad - innerRad;
	//~ double centreRad = (width_rho)/2.0 + innerRad;
	
	int N_rho = 4; 							// 8 layers AG, 4 A2
	int N_z = 30; 							// 90 windings per layer AG, 30 A2
	int N_wires = N_rho*N_z;
	double width_z = 0.03468; 				// 34.6811 mm A2
	double delta_rho = width_rho/N_rho; 	// 533.2 um AG (table value), 1156 um A2 (derived value)
	double delta_z = width_z/N_z; 			// 385.3 um AG (table value), 1280 um A2 (derived value)
	
	// define the line to calculate field on:
	double phi = 0.0;
	double z = 0.0;
	double rho_min = 0.0;
	double rho_max = innerRad*0.8;
	int N_P = 500; //number of points on line
	double BArr[N_P][3];
	for(int k = 0; k < N_P; k++){
		BArr[k][0] = 0.0;
		BArr[k][1] = 0.0;
		BArr[k][2] = 0.0;
	}
	double carPArr[N_P][3];
	
	
	//loop_i.setI(I);
	
	// calcuate field for a single loop
	double BArr_i[N_P][3];
	
	// calculate field for a full mirror coil	
	for(int i = 0; i < N_rho; i++){
		//~ std::cout << i << std::endl;
		
		for(int j = 0; j < N_z; j++){
			// The radius and axial position of the loop is set, so the wires are distributed evenly over the whole volume
			// The current of the loop is set, so the currents in the loop add up to the total current
			Loop loop_i = Loop(innerRad + delta_rho/2.0 + i*delta_rho, I/N_wires, 0, 0, coilCentre[2] - 0.5*width_z + delta_z/2.0 + j*delta_z);
			
			SAMRhoWrapper(loop_i,phi,z,rho_min,rho_max,N_P,carPArr,BArr_i);
			//~ printVec(BArr_i[0],"B_i");
			
			for(int k = 0; k < N_P; k++){
				BArr[k][0] += BArr_i[k][0];
				BArr[k][1] += BArr_i[k][1];
				BArr[k][2] += BArr_i[k][2];
			}		
		}
	}
	
	writeBFieldToFile author(path,"A2");
	author.write(carPArr, BArr, N_P);

}

// Gaussian Quadrature models

void exactMirrorModel_GaussQuad_Loop_Line(std::string path, double carP1[3], double carP2[3]){	
	int N_rho = 4; 								// 8 layers AG, 4 A2
	int N_z = 30; 								// 90 windings per layer AG, 30 A2

	double I = 600; 							// 600 A A2 100 A AG
	
	double innerRad = 0.04125;					// 43.0658 mm AG, 41.25 mm A2
	double outerRad = 0.04637;					// 46.37 mm A2
	double width_rho = outerRad - innerRad;
	double width_z = 0.03468; 					// 34.6811 mm A2
			
	int N_P = 100; 								//number of points on line
		
	// define arrays to hold data
	double BArr_tot[N_P][3];
	double BArr_single[N_P][3];
	for(int k = 0; k < N_P; k++){
		BArr_tot[k][0] = 0.0;
		BArr_tot[k][1] = 0.0;
		BArr_tot[k][2] = 0.0;
		
		BArr_single[k][0] = 0.0;
		BArr_single[k][1] = 0.0;
		BArr_single[k][2] = 0.0;
	}
	double carPArr[N_P][3];
	//~ calcPointsLine3D(carP1,carP2,N_P,carPArr); // outcommented as the McDWrapper function defines the array itself
	
	int NGP_rho = 1; // number of gaussian quadrature points in the radial dimension
	int NGP_z = 3;
	
	double GPRhoValues[NGP_rho];
	double GPZValues[NGP_z];
	
	double GPRhoWeights[NGP_rho];
	double GPZWeights[NGP_z];
	
	const int NWiresRho = N_rho/2; // half the number of wires in each dimension
	const int NWiresZ = N_z/2;
	
	getGaussianQuadratureParams(NGP_rho,GPRhoValues,GPRhoWeights,width_rho/2.0,NWiresRho);
	getGaussianQuadratureParams(NGP_z,GPZValues,GPZWeights,width_z/2.0,NWiresZ);
	
	double centre_rho = innerRad + width_rho/2.0;
	
	for(int nGP_rho = 0; nGP_rho < NGP_rho; nGP_rho++){
		for(int nGP_z = 0; nGP_z < NGP_z; nGP_z++){
			//~ std::cout << "N Rho = " << nGP_rho << std::endl;
			//~ std::cout << "N Z = " << nGP_z << std::endl;
			Loop loop_i(centre_rho + GPRhoValues[nGP_rho],I,0,0,GPZValues[nGP_z]);		
			mcDonaldWrapper(loop_i, 7, carP1, carP2, N_P, carPArr, BArr_single);
			//~ mcDonaldZWrapper(loop_i,7,carP1[0],carP1[1],carP1[2],carP2[2],N_P,carPArr,BArr_single);
			//printVecArr(BArr_single,N_P);
			//SAMWrapper(loop_i, carP1,carP2, N_P, carPArr, BArr_single);
			
			//~ std::cout << "N wires in Z = " << GPZ_NWires[nGP_z] << std::endl;
			//~ std::cout << "Weights Z = " << GPZWeights[nGP_z] << std::endl;
			//~ std::cout << "N wires in Rho = " << GPRho_NWires[nGP_rho] << std::endl;
			//~ std::cout << "Weights Rho = " << GPRhoWeights[nGP_rho] << std::endl;
			
			double GFac = GPZWeights[nGP_z]*GPRhoWeights[nGP_rho];
			//~ std::cout << "Gaussian Quadrature factor = " << GFac << std::endl;
			
			for(int k = 0; k < N_P; k++){
				BArr_tot[k][0] += BArr_single[k][0]*GFac;
				BArr_tot[k][1] += BArr_single[k][1]*GFac;
				BArr_tot[k][2] += BArr_single[k][2]*GFac;
			}
			
			//~ printVecArr(BArr_single,N_P,"B");
			//~ printVec(BArr_single[0],"B1");
			//~ printVec(BArr_single[N_P-1],"B2");
			//~ printVec(carPArr[0],"r1");
			//~ printVec(carPArr[N_P-1],"r2");
			
		}
	}
	
	writeBFieldToFile author_loop(path,"gauss_34");
	author_loop.write(carPArr, BArr_tot, N_P);
	//~ printVecArr(BArr_tot,N_P,"BArr");
	
}

void exactMirrorModel_GaussQuad_Shell(std::string path, double carP1[3], double carP2[3]){	
	int N_rho = 4; 								// 8 layers AG, 4 A2
	int N_z = 30; //used to calc total current	// 90 windings per layer AG, 30 A2

	double i = 600; 							// current in a winding. 100 A AG, 600 A2
		
	double innerRad = 0.04125;					// 43.0658 mm AG, 41.25 mm A2
	double outerRad = 0.04637;					// 46.37 mm A2
	double width_rho = outerRad - innerRad;
	double L = 0.03468; 						// 34.6811 mm A2
	
	int N_P = 100;								//number of points on line
	
	// define arrays to hold data
	double BArr_tot[N_P][3];
	double BArr_single[N_P][3];
	for(int k = 0; k < N_P; k++){
		BArr_tot[k][0] = 0.0;
		BArr_tot[k][1] = 0.0;
		BArr_tot[k][2] = 0.0;
		
		BArr_single[k][0] = 0.0;
		BArr_single[k][1] = 0.0;
		BArr_single[k][2] = 0.0;
	}
	double carPArr[N_P][3];
	
	// calculate Gaussian Quadrature values
	int NGP_rho = 3; // number of gaussian quadrature points in the radial dimension	
	double GPRhoValues[NGP_rho];	
	double GPRhoWeights[NGP_rho];	
	
	const int NElementsRho = N_rho/2; // half the number of layers in each dimension	
	getGaussianQuadratureParams(NGP_rho,GPRhoValues,GPRhoWeights,width_rho/2.0,NElementsRho);
	
	double centre_rho = innerRad + width_rho/2.0;

	for(int nGP_rho = 0; nGP_rho < NGP_rho; nGP_rho++){
		//~ std::cout << "N Rho = " << nGP_rho << std::endl;
		Shell shell_i = Shell(centre_rho + GPRhoValues[nGP_rho],N_z,i,L,0,0,0);
		
		Conway1DWrapper(shell_i,carP1,carP2,N_P,carPArr,BArr_single);
		
		//~ std::cout << "N wires in Rho = " << GPRho_NWires[nGP_rho] << std::endl;
		//~ std::cout << "Weights Rho = " << GPRhoWeights[nGP_rho] << std::endl;
		
		double GFac = GPRhoWeights[nGP_rho];
		//~ std::cout << "Gaussian Quadrature factor = " << GFac << std::endl;
		
		for(int k = 0; k < N_P; k++){
			
			BArr_tot[k][0] += BArr_single[k][0]*GFac;
			BArr_tot[k][1] += BArr_single[k][1]*GFac;
			BArr_tot[k][2] += BArr_single[k][2]*GFac;
		}
		
		//~ printVecArr(BArr_single,N_P,"B");
		//~ printVec(BArr_single[0],"B1");
		//~ printVec(BArr_single[N_P-1],"B2");
		//~ printVec(carPArr[0],"r1");
		//~ printVec(carPArr[N_P-1],"r2");

	}
	writeBFieldToFile author_shell(path,"gauss_4");
	author_shell.write(carPArr, BArr_tot, N_P);
	//~ printVecArr(BArr_tot,N_P,"BArr");

}

void getGaussianQuadratureParams(const int N, double points[], double weights[], const double dimLength, const double NWires){
	if(N == 1){
		points[0] = 0.0;
		
		weights[0] = 2.0*NWires;
	}else if(N == 2){
		points[0] = -1.0/sqrt(3.0)*dimLength;
		points[1] = 1.0/sqrt(3.0)*dimLength;
		
		weights[0] = 1.0*NWires;
		weights[1] = 1.0*NWires;
	}else if(N == 3){
		points[0] = -sqrt(0.6)*dimLength; //sqrt(3/5)
		points[1] = 0;
		points[2] = sqrt(0.6)*dimLength;
		
		weights[0] = 5.0/9.0*NWires;
		weights[1] = 8.0/9.0*NWires;
		weights[2] = 5.0/9.0*NWires;
	}else if(N == 4){
		points[0] = -sqrt(3./7.+2./7.*sqrt(6./5.))*dimLength;
		points[1] = -sqrt(3./7.-2./7.*sqrt(6./5.))*dimLength;
		points[2] = sqrt(3./7.-2./7.*sqrt(6./5.))*dimLength;
		points[3] = sqrt(3./7.+2./7.*sqrt(6./5.))*dimLength;
		
		weights[0] = (18.-sqrt(30.))/36.*NWires;
		weights[1] = (18.+sqrt(30.))/36.*NWires;
		weights[2] = (18.+sqrt(30.))/36.*NWires;
		weights[3] = (18.-sqrt(30.))/36.*NWires;
	}else if(N == 5){
		points[0] = -1./3.*sqrt(5+2*sqrt(10./7.))*dimLength; 
		points[1] = -1./3.*sqrt(5-2*sqrt(10./7.))*dimLength; 
		points[2] = 0.0; 
		points[3] = 1./3.*sqrt(5-2*sqrt(10./7.))*dimLength; 
		points[4] = 1./3.*sqrt(5+2*sqrt(10./7.))*dimLength; 
		
		weights[0] = (322.-13.*sqrt(70))/900.*NWires;
		weights[1] = (322.+13.*sqrt(70))/900.*NWires;
		weights[2] = 128./225.*NWires;
		weights[3] = (322.+13.*sqrt(70))/900.*NWires;
		weights[4] = (322.-13.*sqrt(70))/900.*NWires;
	}else{
		printf("N is not a valid value (In function getGussianQuadratureParams)");
	}
}

// A detailed Biot-Savart model of a solenoid. Wire wound in helix shape

void exactMirrorModel_BiotSavart(std::string path, double carP1[3], double carP2[3]){
	// this model calculates the field along a line from carP1 to carP2 of N layers of helix-shaped wire using
	// the Biot-Savart method
	
	// define the line to calculate field on:	
	const int N_P = 100; //number of points on line
	double carPArr[N_P][3];	
	double BCarArr[N_P][3];
	calcPointsLine3D(carP1,carP2,N_P,carPArr); // calculates points along line
	for(int k = 0; k < N_P; k++){
		BCarArr[k][0] = 0.0;
		BCarArr[k][1] = 0.0;
		BCarArr[k][2] = 0.0;
	}
	
	writeToFile author_segments(path,"Segments");
	writeBFieldToFile author_field(path,"Spiral");
	
	for(int i = 0; i<N_P; i++){
		double carP_i[3];
		double BCar_i[3];
		for(int j=0; j<3; j++){
			carP_i[j] = carPArr[i][j];
			BCar_i[j] = BCarArr[i][j];
		}
		exactMirrorModel_BiotSavart(carP_i, BCar_i);
	}
	
	author_field.write(carPArr,BCarArr,N_P);
}

void exactMirrorModel_BiotSavart(double carP[3], double BCarVec[3]){
	// this model calculates the field at a point car of N layers of helix-shaped wire using
	// the Biot-Savart method
	const double coilCentre_z = 0.0;
	
	const int N_rho = 4; 								// 8 layers AG, 4 A2
	const int N_z = 30; 								// 90 windings per layer AG, 30 A2
	const int N_BS = 1000;							// Number of straight line segments per winding to give to the BS model
	
	assert(N_rho >> 0 && N_z >> 0 && N_BS >> 0); // all integers must be positive

	const double I = 600; 							// 600 A A2 100 A AG
	
	const double innerRad = 0.04125;					// inner coil radius.	43.0658 mm AG, 41.25 mm A2
	const double outerRad = 0.04637;					// outer coil radius. 	46.37 mm A2
	const double width_rho = outerRad - innerRad;		// coil thickness.
	const double width_z = 0.03468; 					// length of coil.		34.6811 mm A2
	
	// calculate the thickness of the wires assuming closely packed wires
	const double delta_rho = width_rho/N_rho;		// 533.2 um AG (table value), 1156 um A2 (derived value)
	const double delta_z = width_z/N_z; 			// 385.3 um AG (table value), 1280 um A2 (derived value)
	
	//~ std::cout << "tube centre = " << coilCentre_z << std::endl;
	//~ std::cout << "R1 = " << innerRad << std::endl;
	//~ std::cout << "R2 = " << outerRad << std::endl;
	//~ std::cout << "thickness = " << width_rho << std::endl;
	//~ std::cout << "L = " << width_z << std::endl;
	//~ std::cout << "delta_rho = " << delta_rho << std::endl;
	//~ std::cout << "delta_z = " << delta_z << std::endl;
	//~ std::cout << "N_rho = " << N_rho << std::endl;
	//~ std::cout << "N_z = " << N_z << std::endl;
	//~ std::cout << "N_BS = " << N_BS << std::endl;
	//~ std::cout << "I = " << I << std::endl;

	//~ std::string path = "";
	//~ writeToFile author_segments(path,"Segments");
	//~ writeBFieldToFile author_field(path,"Spiral");
	
	// calculate field for a full mirror coil	
	for(int n_rho = 0; n_rho < N_rho; n_rho++){				// loop over all layers
		
		double s0[3];
		double s1[3]{innerRad + delta_rho/2.0 + n_rho*delta_rho, 0.0, coilCentre_z + pow(-1,n_rho)*(-width_z/2.0) }; // initial value for s1 (to copy into s0 in the loop) 
		s0[0] = s1[0]; // for a given layer, the rho-coordinate stays the same for all segments
	
		for(int n_z = 0; n_z < N_z; n_z++){				// loop over all windings in a layer
			for(int n_BS = 0; n_BS < N_BS; n_BS++){		// loop over all straight line segments in a winding
								
				s0[1] = s1[1]; 
				s0[2] = s1[2]; // the last end point of a segment is now the start point
				s1[1] = 2.0*ALPHAPhysicalConstants::pi/N_BS*(n_BS+1);
				s1[2] = coilCentre_z + pow(-1,n_rho)*(-width_z/2.0 + delta_z/N_BS*(n_BS+1) + n_z*delta_z);
				
				double s0Car[3], s1Car[3];
				cylPToCarP(s0,s0Car); // convert the start and end point to cartesian coordinates
				cylPToCarP(s1,s1Car);
				//~ author_segments.write(s0Car,3); // write the start coordinate of the segments to a file
				
				biotSavart(s0Car,s1Car,carP,I,BCarVec);		
			}
		}
	}
	//~ author_field.write(carPArr,BArr,N_P);
	//~ printVecArr(BArr,N_P,"BArr");	
}

void testExactMirror(){
	
	double zm = 0.0;
	
	//~ mirrorA = new mirrorMcD_v1(zm);
	mirrorMcD_v1 *mirrorA = new mirrorMcD_v1(zm);
	
  
  
    CartesianVector B{0.,0.,0.};
    CartesianVector position{0.,0.,0.1};
  
  // Calculate B Field here.

  //five_mirrors.concat_Bfield(r, bsol, bf);
  //~ mirrorA->concat_Bfield(position, currents_[0], B);
  mirrorA->concat_Bfield(position, 1.0, B);
  std::cout << B << std::endl;
  
  //~ writeBField("B_GaussQuad_ExactMirror", path2, carPArr, BArr_tot, N_P);
  
}
