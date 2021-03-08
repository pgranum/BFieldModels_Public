#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1 // needed for elliptic integral functions
#include "BFieldModels.h"
#include <iostream>
void biotSavart(double s0[3], double s1[3], double r[3], const double I, double B_vec[3]){

	//~ printVec(s0,"s0");
	//~ printVec(s1,"s1");
	//~ printVec(r,"r");
	//~ printVec(B_vec,"B in");
	
	double i[3]={0.}; double k0[3]={0.}; double k1[3]={0.};
	vecSub(s1,s0,i);
	vecSub(s0,r,k0);
	vecSub(s1,r,k1);
	
	//~ double i_norm; double k0_norm; double k1_norm;
	double i_norm = vecNorm(i);
	double k0_norm = vecNorm(k0);
	double k1_norm = vecNorm(k1);
	
	
	double i_unit[3]={0.}; double k0_unit[3]={0.}; double k1_unit[3]={0.};
	vecMultScal(i, 1.0/i_norm, i_unit);
	vecMultScal(k0,1.0/k0_norm,k0_unit);
	vecMultScal(k1,1.0/k1_norm,k1_unit);
	
	//~ printVec(k0,"k0");
	//~ printVec(k1,"k1");

	double A[3];
	vecCrsP(k0,i_unit,A);
	
	double A_norm;
	A_norm = vecNorm(A);
	vecMultScalOvrwrt(A,1/(A_norm*A_norm));
	
	//~ printVec(A,"A");
	
	double a0; double dk[3];
	vecSub(k1_unit,k0_unit,dk);
	
	a0 = vecDotP(dk,i_unit);
		
	const double C = ALPHAPhysicalConstants::mu0/(4*ALPHAPhysicalConstants::pi)*I*a0;
	vecMultScalOvrwrt(A,C);
	
	//~ std::cout << "C = " << C << std::endl;
	//~ printVec(A,"A*C");
	
	vecAddOvrwrt(B_vec,A);
	
	//~ printVec(B_vec,"B out");
}

void biotSavart_vecClass(const CartesianVector& s0, const CartesianVector& s1, const CartesianVector& r, double I, CartesianVector& B_vec){

	//double b[3]=CartesianVector& a;
	//CartesianVector* a=&b[0]
	//double b[3] = a.GetXYZptr();
	//CartesianVector* a=(CartesianVector*) &b[0];
	
	CartesianVector i=s1-s0;
	CartesianVector k0=s0-r;
	CartesianVector k1=s1-r;
	
	//things like this are also possible
	//double whatever=(s1-s0).GetVectorNorm();
	double i_norm  = i.GetVectorNorm();
	double k0_norm = k0.GetVectorNorm();

	double k1_norm = k1.GetVectorNorm();
	
	
	CartesianVector i_unit = i*(1/i_norm);
	CartesianVector k0_unit = k0*(1/k0_norm);
	CartesianVector k1_unit = k1*(1/k1_norm);

	CartesianVector A=k0.cross(i_unit);

	double A_norm;
	A_norm = A.GetVectorNorm();
	A*=(1/(A_norm*A_norm));
	
	CartesianVector dk=k1_unit-k0_unit;
	
	double a0 = dk.dot(i_unit);
	
	double C = ALPHAPhysicalConstants::mu0/(4*ALPHAPhysicalConstants::pi)*I*a0;
	A *= C;
	
	B_vec += A;
}

void currentLoopFieldOnAxisAnalytic_vecClass(const double R, const double I, const double z, CartesianVector& B_vec){
	// calculates the B-field on the z-axis produced by a current loop in the x-y-plane  using the analytical expression.
	// radius R, zloop is the position of the current loop on the zaxis , current I, zobs is observations point on the z axis
	double k = ALPHAPhysicalConstants::mu0/(4*ALPHAPhysicalConstants::pi)*(2*ALPHAPhysicalConstants::pi*pow(R,2)*I)/(pow(pow(z,2)+pow(R,2),3.0/2.0));
	B_vec.SetXYZ(0.,0.,k);
}

void lineSegmentsOct(double R, double L, double z, double dPhi, const int N, double s_arr[][3]){
	// splits an octupole into straight line segments, which are returned in the array s_arr
	// the octupole has radius R, length L, is centered at z, amd is rotated dPhi
	// each end turn is split into N segments
	// the calculations are done in cylindrical coordinates before being converted to cartesian at the end
	// cylindrical vector = [r, phi, z]
	
	double pi = ALPHAPhysicalConstants::pi;
	double LHalf = L/2.0;
	double piFourths = 2.0*pi/8.0;
	
	// angles at which the straight wire sections are located:
	double phiM[8] = {0.0};
	calcPointsLine1D(0.0+dPhi, 7*piFourths+dPhi, 8, phiM);

	
	
	double sCyl[3] = {R, 0.0, 0.0};
	double sCar[3] = {0.}; // placeholder for cartesian point
	
	double theta_vec[N]; // angles to describe end turn
	
	for(int i = 0; i<8; i++){
		double phi1 = phiM[i]; // end turn start azimuthal angle
		double phi2 = phi1 + piFourths; // end turn end azimuthal angle
		calcPointsLine1D(pi,0,N,theta_vec);
		printArr(theta_vec,N,"theta");
		double turnRad = (phi2 - phi1)*R/2.0; // the radius of the end turn
		
		for(int j = 0; j<N; j++){
			sCyl[1] = cos(theta_vec[j]) * (phi2-phi1)/2.0 + phi1 + (phi2-phi1)/2.0;
			
			if(i%2 == 1){
				// if mod(i,2) = 1, so if i is odd, the loop is curving one way
				sCyl[2] = sin(theta_vec[j])*turnRad + z + LHalf;
			}else if (i%2 == 0){
				// if mod(i,2) = 0, so if i is even, the loop curves the other way
				sCyl[2] = -sin(theta_vec[j])*turnRad + z - LHalf;
			}
			
			cylPToCarP(sCyl,sCar);
			s_arr[i*N+j][0] = sCar[0];
			s_arr[i*N+j][1] = sCar[1];
			s_arr[i*N+j][2] = sCar[2];
		}
	} 
}

void mirrorBFieldLine_BS(Loop loop, double r0[3], double r1[3]){
	
	// calculates the B-field of a mirror along the line from r0 to r1
	// does not return result for now
	
	const int N = 10; // number of observations points
	const int NSegments = 64; // number of BS segments
	double rVecArr[N][3]={0.0};
	double BVecArr[N][3]={0.0};
	calcPointsLine3D(r0,r1,N,rVecArr);
	printVecArr(rVecArr,N,"rVecArr");
	
	for(int i = 0; i<N; i++){
		loopBiotSavart(loop,NSegments,rVecArr[i],BVecArr[i]);
	}
	
	//~ writeBField("test","/eos/user/p/penielse/Test folder",rVecArr,BVecArr,N);
	
	printVecArr(BVecArr,N,"B field array");
	
}

// Current Loop

void loopExactSAM(Loop loop, const double cylP[3], double BCylVec[3]){
/* Calculates the magnet field in BCylVec at the cylindrical coordinate cylVec for a loop of radius R and current I	using the method described by the SAM model
 * 
 * @param loop the loop creating the magnetic field
 * @param cylP the cylindrical coordinate of interest where the magnetfield is calculated
 * @param BCylVec the magnetfield in cylP being calulated in cylindrical coordinates
*/
	
	// Unpacking Loop 
	const double R = loop.getR();

	// Coordinates
	const double rho = cylP[0];
	const double z   = cylP[2]-loop.getz();
	
	// Powers
	const double rho2 = rho*rho;
	const double z2 = z*z; 
	const double R2 = R*R;
	
	// Model
	const double alpha2 = R2 + rho2 + z2 - 2*R*rho;
	const double beta2 = R2 + rho2 + z2 + 2*R*rho;
	const double k = std::sqrt(1 - alpha2/beta2);
	const double beta = std::sqrt(beta2);

	const double E = std::tr1::comp_ellint_2(k);
	const double K = std::tr1::comp_ellint_1(k);
	
	//Preparing Result
	const double C = ALPHAPhysicalConstants::mu0*loop.getI()/ALPHAPhysicalConstants::pi;
	const double denom_rho = 2*alpha2*beta*rho; // Denominator
	const double denom_z = 2*alpha2*beta; // Denominator
	
	double B_rho;
	double B_z;

	if (denom_rho != 0){
		B_rho=C*z/(denom_rho)*((R2+rho2+z2)*E-alpha2*K);
	} else {
		B_rho = 0;
	}
	if (denom_z != 0){
		B_z = C/(denom_z)*((R2-rho2-z2)*E+alpha2*K);
	} else {
		B_z = 0;
	} 	 
	// Result			 
	BCylVec[0]=B_rho;
	BCylVec[1]=0.0;
	BCylVec[2]=B_z;
}

const CartesianVector SAM(Loop loop, CartesianVector& carP){
/* Calculates the magnet field in BCarVec at the cartesian coordinate carVec for a loop of radius R and current I	using the method described by the SAM model
 * 
 * @param loop the loop creating the magnetic field
 * @param cylP the cylindrical coordinate of interest where the magnetfield is calculated
 * @return the magnetfield of the loop at carP being calulated in cartesian coordinates
*/
	
	// Unpacking Loop 
	double R = loop.getR();

	// Coordinates
	double x = carP.GetX();
	double y = carP.GetY();
	double z = carP.GetZ() - loop.getz();
	
		
	// Powers
	double x2 = x*x;
	double y2 = y*y;
	double z2 = z*z;
	double r2 = x2+y2+z2;
	double rho2 = x2+y2;
	double rho = sqrt(rho2);
	double R2 = R*R;
	
	// Model
	double alpha2 = R2 + r2 - 2*R*rho;
	double beta2 = R2 + r2 + 2*R*rho;
	
	double k = std::sqrt(1 - alpha2/beta2);
	double beta = std::sqrt(beta2);

	double E = std::tr1::comp_ellint_2(k);
	double K = std::tr1::comp_ellint_1(k);
	
	//Preparing Result
	double C = ALPHAPhysicalConstants::mu0*loop.getI()/ALPHAPhysicalConstants::pi;
	double denom_x = 2*alpha2*beta*rho2; // Denominator
	double denom_z = 2*alpha2*beta; // Denominator
	
	assert(denom_x != 0);  // To not divide by zero
	assert(denom_z != 0); 
	
	double B_x = C*x*z*denom_x*((R2+r2)*E-alpha2*K);
	double B_y = y/x*B_x;
	double B_z = C*denom_z*((R2-r2)*E+alpha2*K);
	
	return CartesianVector(B_x,B_y,B_z);
}

void loopExactSAM_vecClass(Loop loop, CylindricalVector& cylP, CylindricalVector& BCylVec){
/* Calculates the magnet field in BCylVec at the cylindrical coordinate cylVec for a loop of radius R and current I	using the method described by the SAM model
 * 
 * @param loop the loop creating the magnetic field
 * @param cylP the cylindrical coordinate of interest where the magnetfield is calculated
 * @param BCylVec the magnetfield in cylP being calulated in cylindrical coordinates
*/
	
	//~ CylindricalVector* cylP = (CylindricalVector*)&cylP_[0];
	//~ CylindricalVector* BCylVec = (CylindricalVector*)&BCylVec_[0];
	
	// Unpacking Loop 
	double R = loop.getR();

	// Coordinates
	double rho = cylP.GetRho();
	double z   = cylP.GetZ() - loop.getz();
	
	// Powers
	double rho2 = rho*rho;
	double z2 = z*z; 
	double R2 = R*R;
	
	// Model
	double alpha2 = R2 + rho2 + z2 - 2*R*rho;
	double beta2 = R2 + rho2 + z2 + 2*R*rho;
	double k = std::sqrt(1 - alpha2/beta2);
	double beta = std::sqrt(beta2);

	double E = std::tr1::comp_ellint_2(k);
	double K = std::tr1::comp_ellint_1(k);
	
	//Preparing Result
	double C = ALPHAPhysicalConstants::mu0*loop.getI()/ALPHAPhysicalConstants::pi;
	double denom_rho = 2*alpha2*beta*rho; // Denominator
	double denom_z = 2*alpha2*beta; // Denominator
	
	assert(denom_rho != 0);  // To not divide by zero
	assert(denom_z != 0); 
	
	double B_rho=C*z/(denom_rho)*((R2+rho2+z2)*E-alpha2*K);
	double B_z=C/(denom_z)*((R2-rho2-z2)*E+alpha2*K);	 
	
	BCylVec.SetRhoPhiZ(B_rho,0.0,B_z);
}

// Three McDonalds Model for Loop, Shell, Tube

void loopApproxMcDonald_vecClass(double z, const double I, const double R, const int nmax, const CylindricalVector& cylP, CylindricalVector& BCylVec){	
/* Calculates the magnet field in BCylVec at the cylindrical coordinate cylVec for a loop of radius R and current I	using the method described by the McDonald model
 * 
 * @param loop the loop creating the magnetic field
 * @param nmax an integer [0,7] that the denotes the order of the McDonald calculation
 * @param cylP the cylindrical coordinate of interest where the magnetfield is calculated
 * @param BCylVec the magnetfield in cylP being calulated in cylindrical coordinates
 */

	// Coordinates
	double rho = cylP.GetRho();
	z   = cylP.GetZ() - z;
	
	
	double an[2*nmax+2];
	mcDonaldLoopSupFunc(nmax,z,R,an);
	
	// Preparing terms for loop
	//~ double B_z = an[0];
	BCylVec.SetRhoPhiZ( -(rho/2)*an[1], 0.0, an[0] );
	double constZTerm = 1;
	
	//~ double B_rho = -(rho/2)*an[1];
	double constRTerm = -(rho/2);
	
	// Looping over the series
	for (int n=1; n<=nmax; n++){ 
			constZTerm *= (-1)*pow(rho/2,2)/pow(n,2);
			//~ B_z+= constZTerm*an[2*n];
			BCylVec.AddZ(constZTerm*an[2*n]);
			constRTerm *= (-1)*pow(rho/2,2)*(n)/((n+1)*pow(n,2));
			//~ B_rho += constRTerm*an[2*n+1];
			BCylVec.AddRho(constRTerm*an[2*n+1]);
	}
	
	// Preparing result
	double C = ALPHAPhysicalConstants::mu0*I*R*R/2;
	//~ B_rho *= C;
	BCylVec.MultRho(C);
	//~ B_z *= C;
	BCylVec.MultZ(C);
	
	// Result
	//~ BCylVec[0]=B_rho;
	//~ BCylVec[1]=0.0;
	//~ BCylVec[2]=B_z;
}

void mcDonald(Loop loop, const int nmax, const double cylP[3], double BCylVec[3]){	
/* Calculates the magnet field in BCylVec at the cylindrical coordinate cylVec for a loop of radius R and current I	using the method described by the McDonald model
 * 
 * @param loop the loop creating the magnetic field
 * @param nmax an integer [0,7] that the denotes the order of the McDonald calculation
 * @param cylP the cylindrical coordinate of interest where the magnetfield is calculated
 * @param BCylVec the magnetfield in cylP being calulated in cylindrical coordinates
 */

	// Coordinates
	const double rho = cylP[0];
	const double z   = cylP[2]-loop.getz();
	const double R  = loop.getR();
	const double I  = loop.getI();
	
	double an[2*nmax+2];
	mcDonaldLoopSupFunc(nmax,z,R,an);
	
	// Preparing terms for loop
	double B_z = an[0]; 
	double constZTerm = 1;
	
	const double rho_2 = 0.5*rho;
	double B_rho = -rho_2*an[1];
	double constRTerm = -rho_2;
	
	double rho2_4 = rho_2*rho_2; //(rho/2)²
	// Looping over the series
	for (int n=1; n<=nmax; n++){ 
			constZTerm *= -rho2_4/(n*n); // -(rho/2)²/n²
			B_z+= constZTerm*an[2*n];	
			constRTerm *= -rho2_4/((n+1)*n); //-(rho/2)²n/((n+1)n²)
			B_rho += constRTerm*an[2*n+1];
	}
	
	// Preparing result
	const double C = ALPHAPhysicalConstants::mu0*I*R*R/2;
	B_rho *= C;
	B_z *= C;
	
	// Result
	BCylVec[0]=B_rho;
	BCylVec[1]=0.0;
	BCylVec[2]=B_z;
}

void mcDonaldLoopSupFunc(const int n, const double z, const double R, double an[]){
	
/* This function is a support function to loopApproxMcDonald it generates  the an[] list depending on the order of n to reduce calculation time.
 * 
 * @param loop 	the loop creating the magnetic field
 * @param z		the z-coordinate where the field is being calculated
 * @param R 	the radius of the loop
 * @param an	the list of derivatives of the axial magnetic field
 */
	
	// Powers of z
	const double z2=z*z;  

	// Powers of R
	const double R2 = R*R;
		
	// Powers of d
	const double d = 1/(z2+R2); 
	const double d3_2 = sqrt(d)*d;
	const double d5_2 = d3_2*d;
	
	an[0]=d3_2;
	an[1]=-3*z*d5_2;
if (n > 0){ 
	const double z3=z2*z;
			
	const double d7_2=d5_2*d;
	const double d9_2=d7_2*d;
			
	an[2]=(12*z2-3*R2)*d7_2;
	an[3]=(-60*z3+45*R2*z)*d9_2;
		
if (n > 1){
	const double z4=z3*z;
	const double z5=z4*z;

	const double R4 = R2*R2;
				
	const double d11_2=d9_2*d;
	const double d13_2=d11_2*d;
			
	an[4]=(360*z4-540*R2*z2+45*R4)*d11_2;
	an[5]=(-2520*z5+6300*R2*z3-1575*R4*z)*d13_2;
if (n > 2){
	const double z6=z5*z;
	const double z7=z6*z;

	const double R6 = R4*R2;
		
	const double d15_2=d13_2*d;
	const double d17_2=d15_2*d;


	an[6]=(20160*z6-75600*R2*z4+37800*R4*z2-1575*R6)*d15_2;
	an[7]=(-181440*z7+952560*R2*z5-793800*R4*z3+99225*R6*z)*d17_2;

if (n > 3){
	const double z8=z7*z;
	const double z9=z8*z;
	
	const double R8 = R6*R2;
		
	const double d19_2=d17_2*d;
	const double d21_2=d19_2*d;
	
	an[8]=(1814400*z8-12700800*R2*z6+15876000*R4*z4-3969000*R6*z2+99225*R8)*d19_2;
	an[9]=-(19958400*z9-179625600*R2*z7+314344800*R4*z5-130977000*R6*z3+9823275*R8*z)*d21_2;
if (n > 4){
	const double z10=z9*z;
	const double z11=z10*z;

	const double R10 = R8*R2;

	const double d23_2=d21_2*d;
	const double d25_2=d23_2*d;

	an[10]=(239500800*z10-2694384000*R2*z8+6286896000*R4*z6-3929310000*R6*z4+589396500*R8*z2-9823275*R10)*d23_2;
	an[11]=-(3113510400*z11-42810768000*R2*z9+128432304000*R4*z7-112378266000*R6*z5+28094566500*R8*z3-1404728325*R10*z)*d25_2;
if (n> 5){

	const double z12=z11*z;
	const double z13= z12*z;

	const double R12 = R10*R2;
			
	const double d27_2=d25_2*d;
	const double d29_2=d27_2*d;

	an[12]=(43589145600*z12-719220902400*R2*z10+2697078384000*R4*z8-3146591448000*R6*z6+1179971793000*R8*z4-117997179300*R10*z2+1404728325*R12)*d27_2;
	an[13]=-(653837184000*z13-12749825088000*R2*z11+58436698320000*R4*z9-87655047480000*R6*z7+46018899927000*R8*z5-7669816654500*R10*z3+273922023375*R12*z)*d29_2;
if (n>6){
	const double z14=z13*z;
	const double z15= z14*z;

	const double R14 = R12*R2;
			
	const double d31_2=d29_2*d;
	const double d33_2=d31_2*d;

	an[14]=(10461394944000*z14-237996734976000*R2*z12+1308982042368000*R4*z10-2454341329440000*R6*z8+1718038930608000*R8*z6-429509732652000*R10*z4+30679266618000*R12*z2-273922023375*R14)*d31_2;
	an[15]=-(177843714048000*z15-4668397493760000*R2*z13+30344583709440000*R4*z11-69539671000800000*R6*z9+62585703900720000*R8*z7-21904996365252000*R10*z5+2607737662530000*R12*z3-69850115960625*R14*z)*d33_2;
}
}	
}
}
}
}
}
}

void mcDonald(Shell shell, const int nmax, const double cylP[3], double BCylVec[3]){	
/* Calculates the magnet field in BCylVec at the cylindrical coordinate cylVec for a shell solenoid of radius R and total current I	using the method described by the McDonald model
 * 
 * @param shell the shell creating the magnetic field
 * @param nmax an integer [0,7] that the denotes the order of the McDonald calculation
 * @param cylP the cylindrical coordinate of interest where the magnetfield is calculated
 * @param BCylVec the magnetfield in cylP being calulated in cylindrical coordinates
 */
	
	// Coordinates
	const double rho = cylP[0];
	const double z   = cylP[2]-shell.getz();
	const double R  = shell.getR();
	const double I  = shell.getI();
	const double L  = shell.getL();
	const double L_2  = 0.5*L;
	
	const double z1 = z + L_2;
	const double z2 = z - L_2;
	
	double an1[2*nmax+2];
	double an2[2*nmax+2];
	mcDonaldShellSupFunc(nmax,z1,R,an1);
	mcDonaldShellSupFunc(nmax,z2,R,an2);
	
	// Preparing terms for loop
	double B_z = an1[0]-an2[0]; 
	double constZTerm = 1;
	
	double B_rho = -(rho/2)*(an1[1]-an2[1]);
	double constRTerm = -(rho/2);
	
	// Looping over the series
	for (int n=1; n<=nmax; n++){ 
			constZTerm *= (-1)*pow(rho/2,2)/pow(n,2);
			B_z+= constZTerm*(an1[2*n]-an2[2*n]);
			constRTerm *= (-1)*pow(rho/2,2)*(n)/((n+1)*pow(n,2));
			B_rho += constRTerm*(an1[2*n+1]-an2[2*n+1]);
	}
	
	// Preparing result
	const double C = ALPHAPhysicalConstants::mu0*I/(2*L);
	B_rho *= C;
	B_z *= C;
	
	// Result
	BCylVec[0]=B_rho;
	BCylVec[1]=0.0;
	BCylVec[2]=B_z;
}

void mcDonaldShellSupFunc(const int n, const double z, const double R, double an[]){	
/* This function is a support function to loopApproxMcDonald it generates  the an[] list depending on the order of n to reduce calculation time.
 * 
 * @param n		the order of the McDonald model required
 * @param z		the z-coordinate where the field is being calculated
 * @param R 	the radius of the loop
 * @param an	the list of derivatives of the axial magnetic field
 */
	
	// Powers of z
	const double z2=z*z;  

	// Powers of R
	const double R2 = R*R;
		
	// Powers of d
	const double d = 1/(z2+R2);
	const double d1_2 = sqrt(d); 
	const double d3_2 = R2*d1_2*d;

	an[0]=z*d1_2;
	an[1]=d3_2;
if (n > 0){ 
	const double z2 = z*z;
	
	const double d5_2=d3_2*d;
	const double d7_2=d5_2*d;
	an[2]=-3*z*d5_2;
	an[3]=(12*z2-3*R2)*d7_2;
		
if (n > 1){
	const double z3=z2*z;
	const double z4=z3*z;

	const double R4 = R2*R2;
				
	const double d9_2=d7_2*d;
	const double d11_2=d9_2*d;
	an[4]=(-60*z3+45*R2*z)*d9_2;			
	an[5]=(360*z4-540*R2*z2+45*R4)*d11_2;

if (n > 2){
	const double z5=z4*z;
	const double z6=z5*z;

	const double R6 = R4*R2;
		
	const double d13_2=d11_2*d;
	const double d15_2=d13_2*d;

	an[6]=(-2520*z5+6300*R2*z3-1575*R4*z)*d13_2;
	an[7]=(20160*z6-75600*R2*z4+37800*R4*z2-1575*R6)*d15_2;

if (n > 3){
	const double z7=z6*z;
	const double z8=z7*z;
	
	const double R8 = R6*R2;
		
	const double d17_2=d15_2*d;
	const double d19_2=d17_2*d;
	
	an[8]=(-181440*z7+952560*R2*z5-793800*R4*z3+99225*R6*z)*d17_2;
	an[9]=(1814400*z8-12700800*R2*z6+15876000*R4*z4-3969000*R6*z2+99225*R8)*d19_2;

if (n > 4){
	const double z9=z8*z;
	const double z10=z9*z;

	const double R10 = R8*R2;

	const double d21_2=d19_2*d;
	const double d23_2=d21_2*d;
	
	an[10]=-(19958400*z9-179625600*R2*z7+314344800*R4*z5-130977000*R6*z3+9823275*R8*z)*d21_2;
	an[11]=(239500800*z10-2694384000*R2*z8+6286896000*R4*z6-3929310000*R6*z4+589396500*R8*z2-9823275*R10)*d23_2;
if (n> 5){
	const double z11=z10*z;
	const double z12=z11*z;

	const double R12 = R10*R2;
			
	const double d25_2 = d23_2*d;	
	const double d27_2 = d25_2*d;
	
	an[12]=-(3113510400*z11-42810768000*R2*z9+128432304000*R4*z7-112378266000*R6*z5+28094566500*R8*z3-1404728325*R10*z)*d25_2;
	an[13]=(43589145600*z12-719220902400*R2*z10+2697078384000*R4*z8-3146591448000*R6*z6+1179971793000*R8*z4-117997179300*R10*z2+1404728325*R12)*d27_2;

if (n>6){
	const double z13 = z12*z;
	const double z14 = z13*z;
	
	const double R14 = R12*R2;

	const double d29_2 = d27_2*d;			
	const double d31_2 = d29_2*d;

	an[14]=-(653837184000*z13-12749825088000*R2*z11+58436698320000*R4*z9-87655047480000*R6*z7+46018899927000*R8*z5-7669816654500*R10*z3+273922023375*R12*z)*d29_2;
	an[15]=(10461394944000*z14-237996734976000*R2*z12+1308982042368000*R4*z10-2454341329440000*R6*z8+1718038930608000*R8*z6-429509732652000*R10*z4+30679266618000*R12*z2-273922023375*R14)*d31_2;
}
}	
}
}
}
}
}
}

void mcDonald(Tube tube, const int nmax, const double cylP[3], double BCylVec[3]){	
/* Calculates the magnet field in BCylVec at the cylindrical coordinate cylVec for a finite solenoid of radius R and total current I using the method described by the McDonald model
 * 
 * @param tube 		the tube creating the magnetic field
 * @param nmax 		an integer [0,7] that the denotes the order of the McDonald calculation
 * @param cylP 		the cylindrical coordinate of interest where the magnetfield is calculated
 * @param BCylVec 	the magnetfield in cylP being calulated in cylindrical coordinates
 */

	// Coordinates
	double rho = cylP[0];
	double z   = cylP[2]-tube.getz();  	// Posistion the tube in the center
	double R1  = tube.getR1();			// inner radius
	double R2  = tube.getR2();			// outer radius
	double I  = tube.getI();			// current
	double L  = tube.getL();			// length of solenoid
	double L_2 = 0.5*L;
	
	double z1 = z + L_2;	// The distance from the lower point of the tube to the point of evaluation
	double z2 = z - L_2;	// The distance from the upper point of the tube to the point of evaluation
	
	int arraySize = 2*nmax + 2; // the number of coefficients needed to calculate the requested order (nmax)
	double an11[arraySize];		// array to hold terms with z1 and R1
	double an12[arraySize];		// array to hold terms with z1 and R2
	double an21[arraySize];		// array to hold terms with z2 and R1
	double an22[arraySize];		// array to hold terms with z2 and R2
	
	mcDonaldTubeSupFunc(nmax,z1,R1,an11);
	mcDonaldTubeSupFunc(nmax,z1,R2,an12);
	mcDonaldTubeSupFunc(nmax,z2,R1,an21);
	mcDonaldTubeSupFunc(nmax,z2,R2,an22);
	
	// Preparing terms for loop
	double B_z = an12[0]-an11[0]-an22[0]+an21[0];
	double constZTerm = 1;
	
	double B_rho = -(rho/2)*(an12[1]-an11[1]-an22[1]+an21[1]);
	double constRTerm = -(rho/2);
	
	// Looping over the series
	for (int n=1; n<=nmax; n++){ 
			constZTerm *= (-1)*pow(rho/2,2)/pow(n,2);
			B_z+= constZTerm*(an12[2*n]-an11[2*n]-an22[2*n]+an21[2*n]);
			constRTerm *= (-1)*pow(rho/2,2)*(n)/((n+1)*pow(n,2));
			B_rho += constRTerm*(an12[2*n+1]-an11[2*n+1]-an22[2*n+1]+an21[2*n+1]);
			
			//~ std::cout << "n = " << n << std::endl;
			//~ std::cout << "constZTerm = " << constZTerm << std::endl;
			//~ std::cout << "constRTerm = " << constRTerm << std::endl;
	}
	
	// Preparing result
	double C = ALPHAPhysicalConstants::mu0*I/(2*L*(R2-R1));
	B_rho *= C;
	B_z *= C;
	
	// Result
	BCylVec[0]=B_rho;
	BCylVec[1]=0.0;
	BCylVec[2]=B_z;
}

void mcDonald(Tube tube, const int nmax, const double cylP[3], double BCylVec[3], McD_Tube_Support& McDSupport){	
/* Calculates the magnet field in BCylVec at the cylindrical coordinate cylVec for a finite solenoid of radius R and total current I using the method described by the McDonald model
 * 
 * @param tube 			the tube creating the magnetic field
 * @param nmax 			an integer [0,7] that the denotes the order of the McDonald calculation
 * @param cylP 			the cylindrical coordinate of interest where the magnetfield is calculated
 * @param BCylVec 		the magnetfield in cylP being calulated in cylindrical coordinates
 * @param McDSupport	this function uses a class to precalculate the derivatives of the McD model. This is FASTER AND BETTER
 */

	//~ std::cout <<"McD tube with support and arrays" << std::endl;
	
	// Coordinates
	double rho = cylP[0];
	double z   = cylP[2]-tube.getz();  	// Posistion the tube in the center
	double R1  = tube.getR1();			// inner radius
	double R2  = tube.getR2();			// outer radius
	double I  = tube.getI();			// current
	double L  = tube.getL();			// length of solenoid
	double L_2 = 0.5*L;
	
	double z1 = z + L_2;	// The distance from the lower point of the tube to the point of evaluation
	double z2 = z - L_2;	// The distance from the upper point of the tube to the point of evaluation
	
	int arraySize = 2*nmax + 2; // the number of coefficients needed to calculate the requested order (nmax)
	double an11[arraySize];		// array to hold terms with z1 and R1
	double an12[arraySize];		// array to hold terms with z1 and R2
	double an21[arraySize];		// array to hold terms with z2 and R1
	double an22[arraySize];		// array to hold terms with z2 and R2
	
	McDSupport.getA0(z1,R1,nmax,an11);
	McDSupport.getA0(z1,R2,nmax,an12);
	McDSupport.getA0(z2,R1,nmax,an21);
	McDSupport.getA0(z2,R2,nmax,an22);
	
	//~ printArr(an11,arraySize);
	//~ printArr(an12,arraySize);
	//~ printArr(an21,arraySize);
	//~ printArr(an11,arraySize);
	
	// Preparing terms for loop
	double B_z = an12[0]-an11[0]-an22[0]+an21[0];
	double constZTerm = 1;
	
	double B_rho = -(rho/2)*(an12[1]-an11[1]-an22[1]+an21[1]);
	double constRTerm = -(rho/2);
	
	// Looping over the series
	for (int n=1; n<=nmax; n++){ 
			constZTerm *= (-1)*pow(rho/2,2)/pow(n,2);
			B_z+= constZTerm*(an12[2*n]-an11[2*n]-an22[2*n]+an21[2*n]);
			constRTerm *= (-1)*pow(rho/2,2)*(n)/((n+1)*pow(n,2));
			B_rho += constRTerm*(an12[2*n+1]-an11[2*n+1]-an22[2*n+1]+an21[2*n+1]);
			
			//~ std::cout << "n = " << n << std::endl;
			//~ std::cout << "B = (" << B_rho << ", " << 0 << ", " << B_z << ")" << std::endl;
			//~ std::cout << "constZTerm = " << constZTerm << std::endl;
			//~ std::cout << "constRTerm = " << constRTerm << std::endl;
	}
	
	// Preparing result
	double C = ALPHAPhysicalConstants::mu0*I/(2*L*(R2-R1));
	//~ std::cout << "I = " << I << std::endl;
	//~ std::cout << "L = " << L << std::endl;
	//~ std::cout << "R1 = " << R1 << std::endl;
	//~ std::cout << "R2 = " << R2 << std::endl;
	//~ std::cout << "C = " << C << std::endl;
	B_rho *= C;
	B_z *= C;
	
	// Result
	BCylVec[0]=B_rho;
	BCylVec[1]=0.0;
	BCylVec[2]=B_z;
	
	//~ printVec(BCylVec);
}

void mcDonaldTubeSupFunc(const int n, const double z, const double R, double an[]){	
/* This function is a support function to loopApproxMcDonald it generates  the an[] list depending on the order of n to reduce calculation time. n=3 IS THE MAX
 * 
 * @param n		the order of the McDonald model required
 * @param z		the z-coordinate where the field is being calculated
 * @param R 	the radius of the loop
 * @param an	the list of derivatives of the axial magnetic field
 */
	double R2 = R*R;
	double z2 = z*z;
	double A=sqrt(R2+z2);	// sqrt(R²+z²)
	double B = A + R;       // sqrt(R²+z²)+R
	double b = 1/B;   
	
	double logB = log(B); // ln(sqrt(R²+z²)+z)
	an[0]=z*logB; // z*ln(sqrt(R²+z²)+z)
	double a = 1/A;
	an[1]=logB+z2*a*b;
	
	//~ std::cout << "an_0 = " << an[0] << std::endl;
	//~ std::cout << "an_1 = " << an[1] << std::endl;
	
	if (n > 0){
		double a2 = a*a;
		double a3 = a2*a;
		double b2 = b*b;
		double z3 = z2*z;
		an[2]=3*z*a*b-z3*a3*b-z3*a2*b2;
		
		double a4 = a3*a;
		double a5 = a4*a;
		double b3 = b2*b;
		double z4 = z3*z;
		an[3]=3*a*b-6*z2*a3*b+3*z4*a5*b-6*z2*a2*b2+3*z4*a4*b2+2*z4*a3*b3;
		
		//~ std::cout << "an_2 = " << an[2] << std::endl;
		//~ std::cout << "an_3 = " << an[3] << std::endl;
	if (n > 1){
		double a6 = a5*a;
		double a7 = a6*a;
		double b3 = b2*b;
		double b4 = b3*b;
		double z5 = z4*z;
		an[4]=-15*z*a3*b+30*z3*a5*b-15*z5*a7*b-15*z*a2*b2+30*z3*a4*b2-15*z5*a6*b2+20*z3*a3*b3-12*z5*a5*b3-6*z5*a4*b4;
		double z6 = z5*z;
		double a8 = a7*a;
		double a9 = a8*a;
		double b5 = b4*b;
		an[5]= -15*a3*b+135*z2*a5*b-225*z4*a7*b+105*z6*a9*b-15*a2*b2+135*z2*a4*b2-225*z4*a6*b2+105*z6*a8*b2+90*z2*a3*b3-180*z4*a5*b3+90*z6*a7*b3-90*z4*a4*b4+60*z6*a6*b4+24*z6*a5*b5;
		
		//~ std::cout << "an_4 = " << an[4] << std::endl;
		//~ std::cout << "an_5 = " << an[5] << std::endl;
	if (n > 2){
		double z7 = z6*z;
		double a10 = a9*a;
		double a11 = a10*a;
		double b6 = b5*b;
		an[6]=315*z*a5*b-1575*z3*a7*b+2205*z5*a9*b-945*z7*a11*b+315*z*a4*b2-1575*z3*a6*b2+2205*z5*a8*b2-945*z7*a10*b2+210*z*a3*b3-1260*z3*a5*b3
			+1890*z5*a7*b3-840*z7*a9*b3-630*z3*a4*b4+1260*z5*a6*b4-630*z7*a8*b4+504*z5*a5*b5-360*z7*a7*b5-120*z7*a6*b6;
		double z8 = z7*z;
		double a12 = a11*a;
		double a13 = a12*a;
		double b7 = b6*b;
		an[7]=315*a5*b-6300*z2*a7*b+22050*z4*a9*b-26460*z6*a11*b+10395*z8*a13*b+315*a4*b2-6300*z2*a6*b2+22050*z4*a8*b2-26460*z6*a10*b2+10395*z8*a12*b2
			+210*a3*b3-5040*z2*a5*b3+18900*z4*a7*b3-23520*z6*a9*b3+9450*z8*a11*b3-2520*z2*a4*b4+12600*z4*a6*b4-17640*z6*a8*b4+7560*z8*a10*b4+5040*z4*a5*b5
			-10080*z6*a7*b5+5040*z8*a9*b5-3360*z6*a6*b6+2520*z8*a8*b6+720*z8*a7*b7;
	
		//~ std::cout << "an_6 = " << an[6] << std::endl;
		//~ std::cout << "an_7 = " << an[7] << std::endl;
	}
	}
	}
	//~ std::cout << " " << std::endl;
}	

// Others

void loopBiotSavart(Loop loop, const int NSegements, double carP[3], double BCarVec[3]){
	// calculates the B-field of a current loop in the x-y plane by using the BS method
	// radius R, z position z, current I, number of segments N, obs point carP
	
	// Unfolding the loop
	const double R = loop.getR();
	const double I = loop.getI();
	const double z_loop = loop.getz();

	double B = vecNorm(BCarVec);
	if(B != 0){
		//print("Warning, B is nonzero in loopBiotSavart. Setting B_vec = 0");
		BCarVec[0] = 0;
		BCarVec[1] = 0;
		BCarVec[2] = 0;
	}

	// setup start and end coordinate for straight line segment
	double phis[NSegements];
	calcPointsLine1D(0,2*ALPHAPhysicalConstants::pi,NSegements,phis);
	double x0 = cos(phis[0])*R;
	double y0 = sin(phis[0])*R;
	double x1;
	double y1;
	
	for(int i=0; i<NSegements-1; i++){
		x1 = cos(phis[i+1])*R;
		y1 = sin(phis[i+1])*R;
		double r0[3]{x0,y0,z_loop};
		double r1[3]{x1,y1,z_loop};
			
		biotSavart(r0,r1,carP,I,BCarVec);
		
		x0 = x1;
		y0 = y1;
	}	
}

void loopExactJack510(Loop loop, double cylP[3], double BCylVec[3]){
/* Calculates the magnet field in cylindrical coordinates BCylVec at the cylindircal coordinate cylP for a loop of radius R and current Ibased on the exact solution of exercise 5.10 in Jackson. 
 * 
 * @param R the radius of the loop
 * @param z_loop the z position of the loop
 * @param I the current in the coil loop
 * @param cylP the cylindrical coordinate of interest where the magnetfield is calculated
 * @param BCylVec the magnetic field in cylP calulated in cylindrical coordinates
*/
	// Coordinates
	double rho = cylP[0];
	double z   = cylP[2] - loop.getz();
	
	double I   = loop.getI();
	double R   = loop.getR();
	
	// Initilizing components of the magnetic field
	double B_rho;
	double B_z;
	
	int N = 10000;
	double kArr[N];
	double kMax = 100;
	double kMin = kMax/N;
	linspace(kMin,kMax, N,kArr);
	//printVec(kArr,N,"kArr");

		if(rho == R && z==0){
			std::cout<<"ERROR: cannot evaluate mirror field on mirror"<<std::endl;
			B_rho = 0.0;
			B_z = 0.0;
		} else if(rho < R){
			double integrandRhoArr[N];
			double integrandZArr[N];
			for (int i = 0; i < N-1; i++){
				double kArri=kArr[i];
				double cyl_bessel_k_1_kArri_R = std::tr1::cyl_bessel_k(1,kArri*R);
				//std::cout << "cyl_bessel"<< cyl_bessel_k_1_kArri_R << std:: endl;
				integrandRhoArr[i] = kArri*sin(kArri*z)*std::tr1::cyl_bessel_i(1,kArri*rho)*cyl_bessel_k_1_kArri_R; 
				integrandZArr[i] = kArri*cos(kArri*z)*std::tr1::cyl_bessel_i(0,kArri*rho)*cyl_bessel_k_1_kArri_R;
			}
			B_rho = trapz(N,kArr,integrandRhoArr);
			B_z = trapz(N,kArr,integrandZArr);
		} else if(rho > R){
			double integrandRhoArr[N];
			double integrandZArr[N];
			for (int i = 0; i < N-1; i++){
				double kArri=kArr[i];
				double cyl_bessel_i_1_kArri_R = std::tr1::cyl_bessel_i(1,kArri*R);
				integrandRhoArr[i] = kArri*sin(kArri*z)*std::tr1::cyl_bessel_k(1,kArri*rho)*cyl_bessel_i_1_kArri_R; 
				integrandZArr[i] = kArri*cos(kArri*z)*std::tr1::cyl_bessel_k(0,kArri*rho)*cyl_bessel_i_1_kArri_R;
			}
			B_rho = trapz(N,kArr,integrandRhoArr);
			B_z = trapz(N,kArr,integrandZArr);
		} else{
			std::cout<<"ERROR: cannot evaluate r coordinate"<<std::endl;
			B_rho = 0.0;
			B_z = 0.0;
		};
	
	// Preparing result
	double C = ALPHAPhysicalConstants::mu0*I*R/ALPHAPhysicalConstants::pi;
	B_rho *= C;
	B_z *= C;
	

	// Result
	BCylVec[0] = B_rho;
	BCylVec[1] = 0.0;
	BCylVec[2] = B_z;
}

void loopApproxFrancis(Loop loop, const double lambda, const double sphP[3], double BSphVec[3]){
	
	// an approximation of the B-field of a current loop
	// based on eq A.2 in https://iopscience.iop.org/article/10.1088/1367-2630/14/1/015010
	// input: radius R, mirror z position z, current I, model parameter lambda,
	// observation point r, placeholder for B in cartesian coor B_vec
	
	const double R = loop.getR();
	const double r = sphP[0];
	const double theta = sphP[1];

	// Constants
	const double C = 0.125*ALPHAPhysicalConstants::mu0*loop.getI()*R/lambda;
	const double r2 = r*r;
	const double R2 = R*R;
	
	double B_r;
	double B_theta;

	const double sintheta = sin(theta);
	const double crossterm = 2*R*lambda*r*sintheta;

	// There are some special cases that need to be handled:
	// For theta = 0 or pi, s = 0 and we devide by 0
	// In practice, m12m-m12p is seen to be rounded off to 0 for sin(theta) < 1e-8, which
	// leads to B being wrong by a factor of 1/2
	// Hence, we treat  sin(theta) < 1e-7 to be  sin(theta) = 0
	// r has a similar problem. The result starts to deviate on the 1e-5 level when r < 1e-12

	if(r < 1e-12){
		//~ std::cout << "in centre" << std::endl;
		B_r = ALPHAPhysicalConstants::mu0*loop.getI()/(2*R);		
		B_theta = 0;
	}else if(sin(theta) < 1e-7){ 
		//~ std::cout << "on axis" << std::endl;
		const double m12 = 1/std::sqrt(R2 + r2);
		const double m32 = m12*m12*m12;
		B_r = C*(4*m32*R*lambda)*cos(theta); // the cos(theta) will either be 1 or -1
		//~ std::cout << "theta = " << theta << std::endl;
		//~ std::cout << "cos(theta) = " << cos(theta) << std::endl;
		B_theta = 0;
	} else {
		const double m12m = 1/std::sqrt(R2 + r2 - crossterm);
		const double m12p = 1/std::sqrt(R2 + r2 + crossterm);
		const double m32m = m12m*m12m*m12m;
		const double m32p = m12p*m12p*m12p;
		//~ std::cout << "m12m = " << m12m << std::endl;
		//~ std::cout << "m12p = " << m12p << std::endl;
		//~ std::cout << "m32m = " << m32m << std::endl;
		//~ std::cout << "m32p = " << m32p << std::endl;
		//~ std::cout << "cos(theta) = " << cos(theta) << std::endl;
		//~ std::cout << "sin(theta) = " << sintheta << std::endl;
		//~ std::cout << "r = " << r << std::endl;
		//~ std::cout << "C = " << C << std::endl;
		//~ std::cout << "R = " << R << std::endl;
		//~ std::cout << "lambda = " << lambda << std::endl;
		//~ std::cout << "m12m-m12p = " << m12m-m12p << std::endl;
		//~ std::cout << "T1 = " << (m12m-m12p)/(r*sintheta) << std::endl;
		//~ std::cout << "T2 = " << (m32m+m32p)*R*lambda << std::endl;
		//~ std::cout << "Br = " << C*((m12m-m12p)/(r*sintheta) + (m32m+m32p)*R*lambda)*cos(theta) << std::endl;
		
		B_r = C*((m12m-m12p)/(r*sintheta) + (m32m+m32p)*R*lambda)*cos(theta);
		B_theta = -C*((m12m-m12p)/r -(r-R*lambda*sintheta)*m32m +(r+R*lambda*sintheta)*m32p);
	}
		   
	// Result
	BSphVec[0] = B_r;
	BSphVec[1] = B_theta;
	BSphVec[2] = 0.;
}

void loopApproxJackson(Loop loop, double sphP[3], double BSphVec[3]){
/* Calculates the magnet field in spherical coordinates BSphVec at the spherical coordinate sphP for a loop of radius R and current I based on approximation (5.40) in Jackson Third Edition
 * 
 * @param loop
 * @param sphP the spherical coordinate of interest where the magnetfield is calculated
 * @param BSphVec the magnetic field in cylP calulated in spherical coordinates
*/
	// Shift 
	double carP[3];
	sphPToCarP(sphP,carP);
	carP[2]=carP[2]-loop.getz();
	carPToSphP(carP,sphP);
	
	double I = loop.getI();
	double R = loop.getR();
	
	// Coordinates
	double r = sphP[0];
	double theta = sphP[1];
	
	// Constants
	double sintheta  = sin(theta);
	double costheta  = cos(theta);
	double sintheta2 = sintheta*sintheta;
	double R2 = R*R;
	double r2 = r*r;
	double d = 1/(R2+r2);
	double d2 = d*d;
	double d1_2 = sqrt(d);
	double d3_2 = d1_2*d;
	double d5_2 = d3_2*d;
	
	// Model
	double B_r = 1+3.75*R2*r2*sintheta2*d2;
	double B_z = 2*R2-r2+1.875*R2*r2*sintheta2*(4*R2-3*r2)*d2;
	
	// Preparing Result
	double C_r = 0.5*ALPHAPhysicalConstants::mu0*I*R2*costheta*d3_2;
	B_r *= C_r;
	
	double C_z = -0.25*ALPHAPhysicalConstants::mu0*I*R2*sintheta*d5_2;
	B_z *= C_z;
	
	// Result 
	BSphVec[0]=B_r;
	BSphVec[1]=0.0;
	BSphVec[2]=B_z;
}

void Conway1D(Shell shell, const double cylP[3], double BCylVec[3]){
	// From the "Exact Solution..." paper by J. T. Conway
	
	assert(shell.getx() == 0 && shell.gety() == 0); // the solenoid has to be centered around the axis
	
	const double Z1 = shell.getz() - shell.getL()/2.0;
	const double Z2 = shell.getz() + shell.getL()/2.0;
	const double z = cylP[2];
	const double R = shell.getR();
	
	//~ std::cout << "Z1 = " << Z1 << std::endl;
	//~ std::cout << "Z2 = " << Z2 << std::endl;
	//~ std::cout << "z = " << z << std::endl;
	
	const double IDist = shell.getI()/shell.getL(); // the current density [A/m]
	const double C = ALPHAPhysicalConstants::mu0*IDist*R/2.0;	
	
	if(z < Z1){
		BCylVec[2] = C * ( I_010(R, Z1-z, cylP)
						 - I_010(R, Z2-z, cylP) );
						 
		//~ std::cout << "I_010(Z1-z) =  " << I_010(R, Z1-z, cylP) << std::endl;
		//~ std::cout << "I_010(Z2-z) =  " << I_010(R, Z2-z, cylP) << std::endl;
	}else if(z >= Z1 && z <= Z2){
		BCylVec[2] = C * ( 2* I_010(R, 0.0 , cylP)
							- I_010(R, Z1-z, cylP)
							- I_010(R, Z2-z, cylP) );
							
		//~ std::cout << "I_010(0.0) =  " << I_010(R, 0.0 , cylP) << std::endl;
		//~ std::cout << "I_010(Z1-z) = " << I_010(R, Z1-z, cylP) << std::endl;
		//~ std::cout << "I_010(Z2-z) = " << I_010(R, Z2-z, cylP) << std::endl;
	}else if(z > Z2){
		BCylVec[2] = C * ( I_010(R, Z2-z, cylP)
						 - I_010(R, Z1-z, cylP) );
						 
		//~ std::cout << "I_010(Z2-z) =  " << I_010(R, Z2-z, cylP) << std::endl;
		//~ std::cout << "I_010(Z1-z) =  " << I_010(R, Z1-z, cylP) << std::endl;
	}else{
		std::cout << "Error, cannot determine z in Conway1D" << std::endl;
	}
		
	BCylVec[0] = C * ( I_011(R, Z2-z, cylP)
					 - I_011(R, Z1-z, cylP) );
	
	//~ std::cout << "I_011(Z2-z) = " << I_011(R, Z2-z, cylP) << std::endl;
	//~ std::cout << "I_011(Z1-z) = " << I_011(R, Z1-z, cylP) << std::endl;
}
