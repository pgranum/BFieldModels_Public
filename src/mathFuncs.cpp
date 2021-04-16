#include "mathFuncs.h"

// vector functions

void vecAdd(const double a[3], const double b[3], double c[3]){
/* This function adds two vectors 
* (a+b=c)
* 
* @param a		the left vector being added
* @param b 	 	the right vector being added
* @param c		the resulting vector of the operation
*/
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	c[2] = a[2] + b[2];
}

void vecAddOvrwrt(double a[3], const double b[3]){
/* This function adds two vectors 
* (a+b=a)
* 
* @param a		the left vector of the operation
* @param b 	 	the right vector of the operation 
*/
	a[0] += b[0];
	a[1] += b[1];
	a[2] += b[2];
}

void vecSub(const double a[3], const double b[3], double c[3]){
/* This function subtracts two vectors 
* (a-b=c)
* 
* @param a		the left vector of the operation
* @param b 	 	the right vector of the operation
* @param c		the resulting vector of the operation
*/
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}

void vecMult(const double a[3], const double b[3], double c[3]){
/* This function  multiply two vectors entrywise
* (a*b=c) (entrywise)
* 
* @param a		the left vector of the operation
* @param b 	 	the right vector of the operation
* @param c		the resulting vector of the operation
*/		
	c[0] = a[0] * b[0];
	c[1] = a[1] * b[1];
	c[2] = a[2] * b[2];
}

double vecDotP(const double a[3], const double b[3]){
/* This function finds the dot product of two vectors
* (a*b=return) (dot product)
* 
* @param a		the left vector of the operation
* @param b 	 	the right vector of the operation
* @return		the dot product
*/	
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void vecCrsP(const double  a[3], const double  b[3], double  c[3]){
/* This function  finds the crossproduct of two vectors
* (a*b=c) (cross-product)
* 
* @param a		the left vector of the operation
* @param b 	 	the right vector of the operation
* @param c		the resulting vector of the operation
*/		
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
}

double vecNorm(const double a[3]){
/* This function finds the norm of a vector
* 
* @param a		the vector finding the norm of 
* @return		the norm of vector a
*/	
	return std::sqrt(vecDotP(a,a));
}

void vecMultScal(const double a[3], const double k, double c[3]){
/* This function multiply the vector a with a scalar k
* (k*a=c)
* 
* @param a		the vector of the operation
* @param k 	 	the scalar of the operation
* @param c		the resulting vector of the operation
*/		

	c[0] = a[0]*k;
	c[1] = a[1]*k;
	c[2] = a[2]*k;
}

void vecMultScalOvrwrt(double a[3], const double k){
/* This function multiply the vector a with a scalar k and overwrites a with the resulting vector
* (k*a=a)
* 
* @param a		the vector of the operation and the result
* @param k 	 	the scalar of the operation
*/	
	a[0] *= k;
	a[1] *= k;
	a[2] *= k;
}

// Conversion of points

void cylPToCarP(const double cylP[3], double carP[3]){
	// converts a POINT in cylindrical coordinates to cartesian coordinates
	carP[0] = cylP[0]*cos(cylP[1]);
	carP[1] = cylP[0]*sin(cylP[1]);
	carP[2] = cylP[2];
}

void carPToCylP(const double carP[3], double cylP[3]){
	// converts a POINT in cartesian coordinates to cylindrical coordinates
	//~ ALPHAPhysicalConstants apc;

	cylP[0] = std::sqrt(carP[0]*carP[0] + carP[1]*carP[1]);
	cylP[2] = carP[2];
	
	//~ if (carP[0] < 0 && carP[1] < 0){
		//~ cylP[1] = 0;
	//~ }else if (carP[0] >= 0){
		//~ cylP[1] = atan(carP[1]/carP[0]);
	//~ }else if (carP[0] < 0){
		//~ cylP[1] = -atan(carP[1]/carP[0]) + ALPHAPhysicalConstants::pi;
	//~ }else{
		//~ std::cout << "WARNING: cannot determine phi" << std::endl;
		//~ cylP[1] = std::nan("1");
	//~ };
	if (carP[0] == 0 && carP[1] == 0){
		cylP[1] = 0;
	}else {
		cylP[1] = atan2(carP[1],carP[0]); // arctan(y/x)
	};
	
}

void sphPToCarP(const double sphP[3], double carP[3]){
	// converts a point in spherical coordinates to cartesian coordinates
	
	const double rsintheta = sphP[0]*sin(sphP[1]);
	carP[0]=rsintheta*cos(sphP[2]);
	carP[1]=rsintheta*sin(sphP[2]);
	carP[2]=sphP[0]*cos(sphP[1]);	
}

void carPToSphP(const double carP[3], double sphP[3]){
	// converts a POINT in cartesian coordinates to cylindrical coordinates
	
	const double s2 = carP[0]*carP[0] + carP[1]*carP[1]; 
	
	sphP[0] = std::sqrt(s2 + carP[2]*carP[2]);
	
	if (s2 == 0 && carP[2] == 0){
		sphP[1] = 0;
	}else if (s2 == 0 && carP[2] > 0){
		sphP[1] = 0;
	}else if (s2 == 0 && carP[2] < 0){
		sphP[1] = PhysicsConstants::pi;
	}else if (carP[2] == 0){
		sphP[1] = PhysicsConstants::pi/2.0;
	}else{
		sphP[1] = atan2(sqrt(s2),carP[2]); // arctan(sqrt(y^2 + x^2)/z)
	};
	
	if (carP[0] == 0 && carP[1] == 0){
		sphP[2] = 0;
	}else {
		sphP[2] = atan2(carP[1],carP[0]); // arctan(y/x)
	};
	
}

// Conversion of vectors

void carVecToCylVec(const double carVec[3], const double carP[3], double cylVec[3]){
	// converts a cartesian vector to a cylindrical vector
	
	double cylP[3];
	carPToCylP(carP,cylP);
	
	double cosPhi = cos(cylP[1]);
	double sinPhi = sin(cylP[1]);
	
	cylVec[0]= cosPhi*carVec[0]+sinPhi*carVec[1];
	cylVec[1]=-sinPhi*carVec[0]+cosPhi*carVec[1];
	cylVec[2]= carVec[2];
}

void carVecToSphVec(const double carVec[3], const double carP[3], double sphVec[3]){
	// converst a cartesian vector at the azimuthal angle phi and polar angle theta 
	// to a spherical vector
	double sphP[3];
	carPToSphP(carP,sphP);
	
	double sinTheta=sin(sphP[1]);
	double sinPhi  =sin(sphP[2]);
	double cosTheta=cos(sphP[1]);
	double cosPhi  =cos(sphP[2]);
	sphVec[0]=sinTheta*cosPhi*carVec[0]+sinTheta*sinPhi*carVec[1]+cosTheta*carVec[2];
	sphVec[1]=cosTheta*cosPhi*carVec[0]+cosTheta*sinPhi*carVec[1]-sinTheta*carVec[2];
	sphVec[2]=-sinPhi*carVec[0]+cosPhi*carVec[1];	
}

void cylVecToCarVec(const double cylVec[3], const double cylP[3], double carVec[3]){
	// converts a cylindrical VECTOR at point P to a cartesian vector
	// P should be in cylindrical coordinates
	double cphi = cos(cylP[1]);
	double sphi = sin(cylP[1]);
	
	carVec[0] = cphi*cylVec[0] - sphi*cylVec[1];
	carVec[1] = sphi*cylVec[0] + cphi*cylVec[1];
	carVec[2] = cylVec[2];
}

void sphVecToCarVec(const double sphVec[3], const double sphP[3], double carVec[3]){
	// converts a spherical VECTOR at point P to a cartesian vector
	// P should be in spherical coordinates
	double ctheta = cos(sphP[1]);
	double stheta = sin(sphP[1]);
	double cphi = cos(sphP[2]);
	double sphi = sin(sphP[2]);	
	
	carVec[0] = sphVec[0]*stheta*cphi + sphVec[1]*ctheta*cphi - sphVec[2]*sphi;
	carVec[1] = sphVec[0]*stheta*sphi + sphVec[1]*ctheta*sphi + sphVec[2]*cphi;
	carVec[2] = sphVec[0]*ctheta - sphVec[1]*stheta;
}

// integral function

double I_010(const double R, const double zDiff, const double cylP[3]){
	// A Bessel-Laplace integral from the "Exact Solution..." paper by J. T. Conway
	// Radius of cylindrical shell, R
	// Axial distance between observation point and the centre of the cylindrical shell
	// observation point cylP = (r, phi, z)
	
	//~ std::cout << "z_src = " <<z_src << std::endl;
	
	const double r = cylP[0];
	
	const double rSum = r+R;
	const double rSum2 = rSum*rSum;
	const double rDiff = r-R;
	const double rDiff2 = rDiff*rDiff;
	const double zDiff2 = zDiff*zDiff;
	
	const double k0 = 2.0/sqrt(rSum2 + zDiff2);
	const double k = sqrt(r*R)*k0;
	const double beta = std::asin( zDiff / (sqrt( rDiff2 + zDiff2 )) );
	const double K_comp = std::tr1::comp_ellint_1(k); 
	
	//~ std::cout << "k 010 = " << k << std::endl;
	//~ std::cout << "beta 010 = " << beta << std::endl;
	//~ std::cout << "K_comp 010 = " << K_comp << std::endl;
	
	//////// THE SOLUTION FOR BETA OR K = 0 SHOULD REDUCE TO SOMETHING SIMPLE, BUT THE RESULTS SEEMS TO BE OFF BY ABOUT 1% ///////////
	//~ if(zDiff == 0){
		//~ std::cout << "zDiff = 0 => beta = 0, 010" << std::endl;
		//~ if(r < R){
			//~ return 1.0/R;
		//~ }
		//~ if(r > R){
			//~ return 0.0;
		//~ }
	//~ }
	

	//~ if(r == 0){
		//~ std::cout << "r = 0, 010" << std::endl;
		//~ if(r < R){
			//~ return 1.0/R * (1 - fabs(zDiff)/(2.0*sqrt(R*R+zDiff2)) - fabs(beta)/2.0);
		//~ }
		//~ if(r > R){
			//~ return 1.0/R * (-fabs(zDiff)/(2.0*sqrt(R*R+zDiff2)) - fabs(beta)/2.0);
		//~ }
	//~ }
	
	//~ if(beta == 0){
		//~ std::cout << "beta = 0, 010" << std::endl;
		//~ if(r < R){
			//~ return 1.0/R;
		//~ }
		//~ if(r > R){
			//~ return 0.0;
		//~ }
	//~ }
		
	if(r < R){		
		return 1.0/R * ( 1 - fabs(zDiff)*k0*K_comp/(2*PhysicsConstants::pi) - HeumansLambda(fabs(beta),k)/2.0 );		
	}else if(r > R){		
		return 1.0/R * ( -fabs(zDiff)*k0*K_comp/(2*PhysicsConstants::pi) + HeumansLambda(fabs(beta),k)/2.0 );		
	}else{		
		std::cout << "Error in I_010 function. Cannot evaluate integraol for r=R" << std::endl;
		return 0;
	}
	
}

double I_011(const double R, const double z_src, const double cylP[3]){
	// A Bessel-Laplace integral from the "Exact Solution..." paper by J. T. Conway
	// Radius of cylindrical shell, R
	// Axial distance between observation point and the centre of the cylindrical shell
	// observation point cylP = (r, phi, z)
	
	const double r = cylP[0];
	
	if(r==0){
		return 0.0;
	}
	
	const double rSum = r+R;
	const double rSum2 = rSum*rSum;
	const double zDiff = z_src;
	const double zDiff2 = zDiff*zDiff;
	
	const double k2 = 4*r*R / (rSum2 + zDiff2);
	const double k = sqrt(k2);
	
	//~ std::cout << "k 011 = " << k << std::endl;
	
	const double E_comp = std::tr1::comp_ellint_2(k); // complete elliptic integral of second kind
	const double K_comp = std::tr1::comp_ellint_1(k); // complete elliptic integral of first kind
	
	//~ std::cout << "E_comp 011 = " << E_comp << std::endl;
	//~ std::cout << "K_comp 011 = " << K_comp << std::endl;
	
	return 1.0/(PhysicsConstants::pi*k*sqrt(r*R)) * ( (2-k2)*K_comp - 2*E_comp );
}

double HeumansLambda(const double beta, const double k){
	// Heumans Lambda function. Used in the "Exact Solution..." paper by J. T. Conway
	
	//////// THE SOLUTION FOR BETA OR K = 0 SHOULD REDUCE TO SOMETHING SIMPLE, BUT THE RESULTS SEEMS TO BE OFF BY ABOUT 1% ///////////
	//~ if(beta == 0){
		//~ return 0.0;
	//~ }
	//~ if(k == 0){
		//~ return beta;
	//~ }	
	
	const double kPrime = sqrt(1.0 - k*k);
	const double E_comp = std::tr1::comp_ellint_2(k); // complete elliptic integral of second kind
	const double K_comp = std::tr1::comp_ellint_1(k); // complete elliptic integral of first kind
	const double E = std::tr1::ellint_2(kPrime,beta); // elliptic integral of second kind
	const double F = std::tr1::ellint_1(kPrime,beta); // elliptic integral of first kind
	
	//~ std::cout << "k = " << k << std::endl;
	//~ std::cout << "beta = " << beta << std::endl;
	//~ std::cout << "E = " << E << std::endl;
	//~ std::cout << "F = " << F << std::endl;
	
	return 2.0/PhysicsConstants::pi * (E_comp*F + K_comp*E - K_comp*F);
}
