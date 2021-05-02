#include "mathFuncs.h"

// vector functions

//~ void vecAdd(const double a[3], const double b[3], double c[3]){
//~ /* This function adds two vectors 
//~ * (a+b=c)
//~ * 
//~ * @param a		the left vector being added
//~ * @param b 	 	the right vector being added
//~ * @param c		the resulting vector of the operation
//~ */
	//~ c[0] = a[0] + b[0];
	//~ c[1] = a[1] + b[1];
	//~ c[2] = a[2] + b[2];
//~ }

//~ void vecAddOvrwrt(double a[3], const double b[3]){
//~ /* This function adds two vectors 
//~ * (a+b=a)
//~ * 
//~ * @param a		the left vector of the operation
//~ * @param b 	 	the right vector of the operation 
//~ */
	//~ a[0] += b[0];
	//~ a[1] += b[1];
	//~ a[2] += b[2];
//~ }

//~ void vecSub(const double a[3], const double b[3], double c[3]){
//~ /* This function subtracts two vectors 
//~ * (a-b=c)
//~ * 
//~ * @param a		the left vector of the operation
//~ * @param b 	 	the right vector of the operation
//~ * @param c		the resulting vector of the operation
//~ */
	//~ c[0] = a[0] - b[0];
	//~ c[1] = a[1] - b[1];
	//~ c[2] = a[2] - b[2];
//~ }

//~ void vecMult(const double a[3], const double b[3], double c[3]){
//~ /* This function  multiply two vectors entrywise
//~ * (a*b=c) (entrywise)
//~ * 
//~ * @param a		the left vector of the operation
//~ * @param b 	 	the right vector of the operation
//~ * @param c		the resulting vector of the operation
//~ */		
	//~ c[0] = a[0] * b[0];
	//~ c[1] = a[1] * b[1];
	//~ c[2] = a[2] * b[2];
//~ }

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

//~ void vecCrsP(const double  a[3], const double  b[3], double  c[3]){
//~ /* This function  finds the crossproduct of two vectors
//~ * (a*b=c) (cross-product)
//~ * 
//~ * @param a		the left vector of the operation
//~ * @param b 	 	the right vector of the operation
//~ * @param c		the resulting vector of the operation
//~ */		
	//~ c[0] = a[1]*b[2] - a[2]*b[1];
	//~ c[1] = a[2]*b[0] - a[0]*b[2];
	//~ c[2] = a[0]*b[1] - a[1]*b[0];
//~ }

double vecNorm(const double a[3]){
/* This function finds the norm of a vector
* 
* @param a		the vector finding the norm of 
* @return		the norm of vector a
*/	
	return std::sqrt(vecDotP(a,a));
}

//~ void vecMultScal(const double a[3], const double k, double c[3]){
//~ /* This function multiply the vector a with a scalar k
//~ * (k*a=c)
//~ * 
//~ * @param a		the vector of the operation
//~ * @param k 	 	the scalar of the operation
//~ * @param c		the resulting vector of the operation
//~ */		

	//~ c[0] = a[0]*k;
	//~ c[1] = a[1]*k;
	//~ c[2] = a[2]*k;
//~ }

//~ void vecMultScalOvrwrt(double a[3], const double k){
//~ /* This function multiply the vector a with a scalar k and overwrites a with the resulting vector
//~ * (k*a=a)
//~ * 
//~ * @param a		the vector of the operation and the result
//~ * @param k 	 	the scalar of the operation
//~ */	
	//~ a[0] *= k;
	//~ a[1] *= k;
	//~ a[2] *= k;
//~ }

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


