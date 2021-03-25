#include "mcDTubeSupportClass.h"

void McD_Tube_Support::getA0(const double z, const double R, double a0_array[]) const{
	// This function calculates the the a_n terms for a given z and R. 
	// The derivatives of a_0 have been precalculated by the class constructor, which
	// was instantiated to be able to go to a given order. The requested order has to
	// be below this
	// various outcommented debugging code has been left here
	
	//Constructor makes this check
	//assert(n_McD <= McD_order); // make sure that the requested order is not bigger then the max order the class is configured to calculate
	
	const int n = 2 + 2*McD_order;
	//~ std::cout << "number of a_n terms to calculate, n = " << n <<std::endl;
	for(int i=0; i<n; i++){
		a0_array[i] = 0; //make sure all elements are 0
	}	
	
	const double R2 = R*R;
	const double z2 = z*z;
	const double A=sqrt(R2+z2);	// sqrt(R²+z²)
	const double B = A + R;       // sqrt(R²+z²)+R
	const double b = 1/B;   
	const double a = 1/A;
	const double logB = log(B); 	// ln(sqrt(R²+z²)+z)
	
	a0_array[0] = z*logB; 			// z*ln(sqrt(R²+z²)+z)
	a0_array[1] = logB+z2*a*b;
	
	//~ std::cout << "z = " << z << std::endl;
	//~ std::cout << "A = " << A << std::endl;
	//~ std::cout << "B = " << B << std::endl;
	//~ std::cout << "z*logB = " << z*logB << std::endl;
		
	
	if(n>1){
	const int z_power_needed = 3 + (n-2);
	const int a_power_needed = 3 + 2*(n-2);
	const int b_power_needed = 2 + (n-2);
	
	double z_powers[z_power_needed];
	double a_powers[a_power_needed];
	double b_powers[b_power_needed];
	
	for(int i = 0; i<z_power_needed; i++){
		if(i==0){
			z_powers[i] = 1;
		}else{
			z_powers[i] = z_powers[i-1]*z;
		}		
	}
	for(int i = 0; i<a_power_needed; i++){
		if(i==0){
			a_powers[i] = 1;
		}else{
			a_powers[i] = a_powers[i-1]*a;
		}		
	}
	for(int i = 0; i<b_power_needed; i++){
		if(i==0){
			b_powers[i] = 1;
		}else{
			b_powers[i] = b_powers[i-1]*b;
		}		
	}
	
	//~ printStdVec(k_arr, "k");
	//~ printStdVec(lambda_arr, "lambda");
	//~ printStdVec(mu_arr, "mu");
	//~ printStdVec(nu_arr, "nu");
	//~ printStdVec(Ti_arr, "T_i");	

	//~ printArr(z_powers, z_power_needed, "z^n");
	//~ printArr(a_powers, a_power_needed, "a^n");
	//~ printArr(b_powers, b_power_needed, "b^n");
	
	//~ std::cout << "an_0 = " << a0_array[0] << std::endl;
	//~ std::cout << "an_1 = " << a0_array[1] << std::endl;
	
	
	
	//~ std::cout << "total number of terms, NT = " << NT << std::endl;
	for(int i=0; i<NT; i++){
		a0_array[Ti_arr[i]] += k_arr[i]*z_powers[lambda_arr[i]]*a_powers[mu_arr[i]]*b_powers[nu_arr[i]];
		//~ a0_array[Ti_arr.at(i)] += k_arr.at(i)*z_powers[lambda_arr.at(i)]*a_powers[mu_arr.at(i)]*b_powers[nu_arr.at(i)];	
		//~ std::cout << "T_i = " << k_arr.at(i)*z_powers[lambda_arr.at(i)]*a_powers[mu_arr.at(i)]*b_powers[nu_arr.at(i)] << std::endl;
	}
	
	//~ std::cout << "an_0 = " << a0_array[0] << std::endl;
	//~ std::cout << "an_1 = " << a0_array[1] << std::endl;
	//~ std::cout << "an_2 = " << a0_array[2] << std::endl;
	//~ std::cout << "an_3 = " << a0_array[3] << std::endl;
	//~ std::cout << "an_4 = " << a0_array[4] << std::endl;
	//~ std::cout << "an_5 = " << a0_array[5] << std::endl;
	//~ std::cout << "an_6 = " << a0_array[6] << std::endl;
	//~ std::cout << "an_7 = " << a0_array[7] << std::endl;
	//~ std::cout << " " << std::endl;
	
	}	
}

void McD_Tube_Support::getB(Tube tube, const double cylP[3], double BCylVec[3]) const{	
	/* Calculates the magnet field in BCylVec at the cylindrical coordinate cylVec for a finite solenoid of radius R and total current I using the method described by the McDonald model
	 * 
	 * @param tube 			the tube creating the magnetic field
	 * @param cylP 			the cylindrical coordinate of interest where the magnetfield is calculated
	 * @param BCylVec 		the magnetfield in cylP being calulated in cylindrical coordinates
	 * @param McDSupport	this function uses a class to precalculate the derivatives of the McD model. This is FASTER AND BETTER
	 */
	
	//~ std::cout <<"McD tube with support and Vectors in the class" << std::endl;
	//~ std::cout << "P = (" << cylP.GetRho() << ", " << cylP.GetPhi() << ", " << cylP.GetZ() << ")" << std::endl;
	//~ std::cout << "B = (" << BCylVec.GetRho() << ", " << BCylVec.GetPhi() << ", " << BCylVec.GetZ() << ")" << std::endl;
	
	// Coordinates
	const double rho = cylP[0];
	const double z   = cylP[2]-tube.getz();  	// Posistion the tube in the center
	const double R1  = tube.getR1();				// inner radius
	const double R2  = tube.getR2();				// outer radius
	const double I  = tube.getI();					// current
	//~ std::cout << "getB called with I = " << I << std::endl;
	const double L  = tube.getL();					// length of solenoid
	const double L_2 = 0.5*L;
	
	const double z1 = z + L_2;	// The distance from the lower point of the tube to the point of evaluation
	const double z2 = z - L_2;	// The distance from the upper point of the tube to the point of evaluation
	
	const int arraySize = 2*McD_order + 2; 	// the number of coefficients needed to calculate the requested order (nmax)
	double an11[arraySize];		// array to hold terms with z1 and R1
	double an12[arraySize];		// array to hold terms with z1 and R2
	double an21[arraySize];		// array to hold terms with z2 and R1
	double an22[arraySize];		// array to hold terms with z2 and R2
	
	this->getA0(z1,R1,an11);
	this->getA0(z1,R2,an12);
	this->getA0(z2,R1,an21);
	this->getA0(z2,R2,an22);
	
	//~ printArr(an11,arraySize);
	//~ printArr(an12,arraySize);
	//~ printArr(an21,arraySize);
	//~ printArr(an11,arraySize);
	
	// Preparing terms for loop
	double B_z = an12[0]-an11[0]-an22[0]+an21[0];
	//~ BCylVec.SetZ(an12[0]-an11[0]-an22[0]+an21[0]);
	double constZTerm = 1;
	
	double B_rho = -(rho/2)*(an12[1]-an11[1]-an22[1]+an21[1]);
	//~ BCylVec.SetRho(-(rho/2)*(an12[1]-an11[1]-an22[1]+an21[1]));
	double constRTerm = -(rho/2);
	
	// Looping over the series
	for (int n=1; n<=McD_order; n++){ 
			constZTerm *= (-1)*pow(rho/2,2)/pow(n,2);
			B_z+= constZTerm*(an12[2*n]-an11[2*n]-an22[2*n]+an21[2*n]);
			//~ BCylVec.AddZ(constZTerm*(an12[2*n]-an11[2*n]-an22[2*n]+an21[2*n]));
			constRTerm *= (-1)*pow(rho/2,2)*(n)/((n+1)*pow(n,2));
			B_rho += constRTerm*(an12[2*n+1]-an11[2*n+1]-an22[2*n+1]+an21[2*n+1]);
			//~ BCylVec.AddRho(constRTerm*(an12[2*n+1]-an11[2*n+1]-an22[2*n+1]+an21[2*n+1]));
			
			//~ std::cout << "n = " << n << std::endl;
			//~ std::cout << "B = (" << BCylVec.GetRho() << ", " << BCylVec.GetPhi() << ", " << BCylVec.GetZ() << ")" << std::endl;
			//~ std::cout << "constZTerm = " << constZTerm << std::endl;
			//~ std::cout << "constRTerm = " << constRTerm << std::endl;
	}
	
	// Preparing result
	const double C = PhysicsConstants::mu0*I/(2*L*(R2-R1));
	//~ std::cout << "I = " << I << std::endl;
	//~ std::cout << "L = " << L << std::endl;
	//~ std::cout << "R1 = " << R1 << std::endl;
	//~ std::cout << "R2 = " << R2 << std::endl;
	//~ std::cout << "C = " << C << std::endl;
	B_rho *= C;
	//~ BCylVec.MultRho(C);
	B_z *= C;
	//~ BCylVec.MultZ(C);
	
	//~ std::cout << "P = (" << cylP.GetRho() << ", " << cylP.GetPhi() << ", " << cylP.GetZ() << ")" << std::endl;
	//~ std::cout << "B = (" << BCylVec.GetRho() << ", " << BCylVec.GetPhi() << ", " << BCylVec.GetZ() << ")" << std::endl;
	
	// Result
	BCylVec[0]=B_rho;
	BCylVec[1]=0.0;
	BCylVec[2]=B_z;
}
