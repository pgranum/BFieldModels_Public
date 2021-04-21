#include "BFieldModels.h"

// Utility Functions

void BiotSavart::getB(const double carP[3], const double I, double BCarVec[3]) const {

	//~ printVec(carS0,"carS0");
	//~ printVec(carS1,"carS1");
	//~ printVec(carP,"carP");
	//~ printVec(BCarVec,"B in");
	
	double i[3]={carS1[0]-carS0[0], carS1[1]-carS0[1], carS1[2]-carS0[2]}; 
	double k0[3]={carS0[0]-carP[0], carS0[1]-carP[1], carS0[2]-carP[2]}; 
	double k1[3]={carS1[0]-carP[0], carS1[1]-carP[1], carS1[2]-carP[2]}; 
	
	//~ double k0[3]={0.}; 
	//~ double k1[3]={0.};
	//~ vecSub(carS1,carS0,i);
	//~ vecSub(carS0,carP,k0);
	//~ vecSub(carS1,carP,k1);
	
	//~ double i_norm; double k0_norm; double k1_norm;
	const double i_norm = sqrt(i[0]*i[0] + i[1]*i[1] + i[2]*i[2]);
	const double k0_norm = sqrt(k0[0]*k0[0] + k0[1]*k0[1] + k0[2]*k0[2]);
	const double k1_norm = sqrt(k1[0]*k1[0] + k1[1]*k1[1] + k1[2]*k1[2]);
	
	// the inverse norm
	const double i_inorm = 1.0/i_norm;
	const double k0_inorm = 1.0/k0_norm;
	const double k1_inorm = 1.0/k1_norm;	
	
	const double i_unit[3] = {i[0]*i_inorm, i[1]*i_inorm, i[2]*i_inorm}; 
	const double k0_unit[3] = {k0[0]*k0_inorm, k0[1]*k0_inorm, k0[2]*k0_inorm};
	const double k1_unit[3] = {k1[0]*k1_inorm, k1[1]*k1_inorm, k1[2]*k1_inorm};
	//~ vecMultScal(i, 1.0/i_norm, i_unit);
	//~ vecMultScal(k0,1.0/k0_norm,k0_unit);
	//~ vecMultScal(k1,1.0/k1_norm,k1_unit);
	
	//~ printVec(k0,"k0");
	//~ printVec(k1,"k1");

	double A[3] = {k0[1]*i_unit[2] - k0[2]*i_unit[1],
				   k0[2]*i_unit[0] - k0[0]*i_unit[2],
				   k0[0]*i_unit[1] - k0[1]*i_unit[0]};
	//~ vecCrsP(k0,i_unit,A);
	
	const double A_inorm = 1.0/sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
	const double A_inorm2 = A_inorm*A_inorm;
	//~ A_norm = vecNorm(A);
	//~ vecMultScalOvrwrt(A,1/(A_norm*A_norm));
	A[0] = A[0]*A_inorm2;
	A[1] = A[1]*A_inorm2;
	A[2] = A[2]*A_inorm2;
	
	//~ printVec(A,"A");
	
	
	const double dk[3] = {k1_unit[0]-k0_unit[0], k1_unit[1]-k0_unit[1], k1_unit[2]-k0_unit[2]};
	//~ vecSub(k1_unit,k0_unit,dk);
	
	const double a0 = dk[0]*i_unit[0] + dk[1]*i_unit[1] + dk[2]*i_unit[2]; 
	//~ a0 = vecDotP(dk,i_unit);
		
	const double C = PhysicsConstants::mu0/(4*PhysicsConstants::pi)*I*a0;
	BCarVec[0] = A[0]*C;
	BCarVec[1] = A[1]*C;
	BCarVec[2] = A[2]*C;
	//~ vecMultScalOvrwrt(A,C);
	
	//~ std::cout << "C = " << C << std::endl;
	//~ printVec(A,"A*C");
	
	//~ vecAddOvrwrt(BCarVec,A);
	
	//~ printVec(BCarVec,"B out");
}

void GQ_Support::getGaussianQuadratureParams(const int N, std::vector<double> &points, std::vector<double> &weights, const double dimLength, const double NWires) const {
	if(N == 1){
		points.push_back(0.0);
		
		weights.push_back(2.0*NWires);
	}else if(N == 2){
		points.push_back(-1.0/sqrt(3.0)*dimLength);
		points.push_back(1.0/sqrt(3.0)*dimLength);
		
		weights.push_back(1.0*NWires);
		weights.push_back(1.0*NWires);
	}else if(N == 3){
		points.push_back(-sqrt(0.6)*dimLength); //sqrt(3/5)
		points.push_back(0);
		points.push_back(sqrt(0.6)*dimLength);
		
		weights.push_back(5.0/9.0*NWires);
		weights.push_back(8.0/9.0*NWires);
		weights.push_back(5.0/9.0*NWires);
	}else if(N == 4){
		points.push_back(-sqrt(3./7.+2./7.*sqrt(6./5.))*dimLength);
		points.push_back(-sqrt(3./7.-2./7.*sqrt(6./5.))*dimLength);
		points.push_back(sqrt(3./7.-2./7.*sqrt(6./5.))*dimLength);
		points.push_back(sqrt(3./7.+2./7.*sqrt(6./5.))*dimLength);
		
		weights.push_back((18.-sqrt(30.))/36.*NWires);
		weights.push_back((18.+sqrt(30.))/36.*NWires);
		weights.push_back((18.+sqrt(30.))/36.*NWires);
		weights.push_back((18.-sqrt(30.))/36.*NWires);
	}else if(N == 5){
		points.push_back(-1./3.*sqrt(5+2*sqrt(10./7.))*dimLength); 
		points.push_back(-1./3.*sqrt(5-2*sqrt(10./7.))*dimLength); 
		points.push_back(0.0); 
		points.push_back(1./3.*sqrt(5-2*sqrt(10./7.))*dimLength); 
		points.push_back(1./3.*sqrt(5+2*sqrt(10./7.))*dimLength); 
		
		weights.push_back((322.-13.*sqrt(70))/900.*NWires);
		weights.push_back((322.+13.*sqrt(70))/900.*NWires);
		weights.push_back(128./225.*NWires);
		weights.push_back((322.+13.*sqrt(70))/900.*NWires);
		weights.push_back((322.-13.*sqrt(70))/900.*NWires);
	}else{
		std::cout << "N is not a valid value (In function getGussianQuadratureParams)" << std::endl;
	}
}

// Current Loop

void SimpleAnalyticModel::getB(const double cylP[3], double BCylVec[3]) const {
/* Calculates the magnet field in BCylVec at the cylindrical coordinate cylVec for a loop of radius R and current I	using the method described by the SAM model
 * 
 * @param loop the loop creating the magnetic field
 * @param cylP the cylindrical coordinate of interest where the magnetfield is calculated
 * @param BCylVec the magnetfield in cylP being calulated in cylindrical coordinates
*/
	// Coordinates
	const double rho = cylP[0];
	const double z   = cylP[2]-Loop::z;
	
	// Powers
	const double rho2 = rho*rho;
	const double z2 = z*z; 
	
	// Model
	const double alpha2 = R2 + rho2 + z2 - 2*R*rho;
	const double beta2 = R2 + rho2 + z2 + 2*R*rho;
	const double k = std::sqrt(1 - alpha2/beta2);
	const double beta = std::sqrt(beta2);

	const double E = std::tr1::comp_ellint_2(k);
	const double K = std::tr1::comp_ellint_1(k);
	
	//Preparing Result
	const double C = PhysicsConstants::mu0*I/PhysicsConstants::pi;
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

void BiotSavart_Loop::getB(const double carP[3], double BCarVec[3]) const {
	// calculates the B-field of a current loop in the x-y plane by using the BS method
	// radius R, z position z, current I, number of segments N, obs point carP

	// make sure that the placeholder for the magnetic field is 0
	BCarVec[0] = 0;
	BCarVec[1] = 0;
	BCarVec[2] = 0;

	//~ double x0 = xs[0];
	//~ double y0 = ys[0];
	//~ double x1;
	//~ double y1;
	//~ printVec(BCarVec,"B"); 
	
	for(int n_BS=0; n_BS<N_BS; n_BS++){
		//~ x1 = xs[i+1];
		//~ y1 = ys[i+1];
		//~ double r0[3]{x0,y0,Loop::z};
		//~ double r1[3]{x1,y1,Loop::z};
			
		//~ this->BSSegment(r0,r1,carP,I,BCarVec);
		segments[n_BS].getB(carP, I, BCarVec);
		
		//~ printVec(BCarVec,"B"); 
		
		//~ x0 = x1;
		//~ y0 = y1;
	}	
}

void McD_Loop::getB(const double cylP[3], double BCylVec[3]) const {	
/* Calculates the magnet field in BCylVec at the cylindrical coordinate cylVec for a loop of radius R and current I	using the method described by the McDonald model
 * 
 * @param loop the loop creating the magnetic field
 * @param McDOrder an integer [0,7] that the denotes the order of the McDonald calculation
 * @param cylP the cylindrical coordinate of interest where the magnetfield is calculated
 * @param BCylVec the magnetfield in cylP being calulated in cylindrical coordinates
 */

	// Coordinates
	const double rho = cylP[0];
	const double z   = cylP[2]-Loop::z;
	
	double an[2*McDOrder+2];
	mcDonaldLoopSupFunc(z,an);
	
	// Preparing terms for loop
	double B_z = an[0]; 
	double constZTerm = 1;
	
	const double rho_2 = 0.5*rho;
	double B_rho = -rho_2*an[1];
	double constRTerm = -rho_2;
	
	//~ std::cout << "B_z_i_Acc = " << B_z << std::endl;
	//~ std::cout << "B_rho_i_Acc = " << B_rho << std::endl;
	
	double rho2_4 = rho_2*rho_2; //(rho/2)²
	// Looping over the series
	for (int n=1; n<=McDOrder; n++){ 
		constZTerm *= -rho2_4/(n*n); // -(rho/2)²/n²
		B_z+= constZTerm*an[2*n];	
		constRTerm *= -rho2_4/((n+1)*n); //-(rho/2)²n/((n+1)n²)
		B_rho += constRTerm*an[2*n+1];
		
		//~ std::cout << "B_z_i = " << constZTerm*an[2*n] << std::endl;
		//~ std::cout << "B_rho_i = " << constRTerm*an[2*n+1] << std::endl;
		//~ std::cout << "B_z_i_Acc = " << B_z << std::endl;
		//~ std::cout << "B_rho_i_Acc = " << B_rho << std::endl;
	}
	
	// Preparing result
	//~ std::cout << "C = " << C << std::endl;
	B_rho *= C;
	B_z *= C;
	
	// Result
	BCylVec[0]=B_rho;
	BCylVec[1]=0.0;
	BCylVec[2]=B_z;
}

void McD_Loop::mcDonaldLoopSupFunc(const double z, double an[]) const {
	
/* This function is a support function to loopApproxMcDonald it generates  the an[] list depending on the order of n to reduce calculation time.
 * 
 * @param loop 	the loop creating the magnetic field
 * @param z		the z-coordinate where the field is being calculated
 * @param R 	the radius of the loop
 * @param an	the list of derivatives of the axial magnetic field
 */
	
	// Powers of z
	const double z2=z*z;  
		
	// Powers of d
	const double d = 1/(z2+R2); 
	const double d3_2 = sqrt(d)*d;
	const double d5_2 = d3_2*d;
	
	an[0]=d3_2;
	an[1]=-3*z*d5_2;
if (McDOrder > 0){ 
	const double z3=z2*z;
			
	const double d7_2=d5_2*d;
	const double d9_2=d7_2*d;
			
	an[2]=(12*z2-3*R2)*d7_2;
	an[3]=(-60*z3+45*R2*z)*d9_2;
		
if (McDOrder > 1){
	const double z4=z3*z;
	const double z5=z4*z;

	const double R4 = R2*R2;
				
	const double d11_2=d9_2*d;
	const double d13_2=d11_2*d;
			
	an[4]=(360*z4-540*R2*z2+45*R4)*d11_2;
	an[5]=(-2520*z5+6300*R2*z3-1575*R4*z)*d13_2;
if (McDOrder > 2){
	const double z6=z5*z;
	const double z7=z6*z;

	const double R6 = R4*R2;
		
	const double d15_2=d13_2*d;
	const double d17_2=d15_2*d;


	an[6]=(20160*z6-75600*R2*z4+37800*R4*z2-1575*R6)*d15_2;
	an[7]=(-181440*z7+952560*R2*z5-793800*R4*z3+99225*R6*z)*d17_2;

if (McDOrder > 3){
	const double z8=z7*z;
	const double z9=z8*z;
	
	const double R8 = R6*R2;
		
	const double d19_2=d17_2*d;
	const double d21_2=d19_2*d;
	
	an[8]=(1814400*z8-12700800*R2*z6+15876000*R4*z4-3969000*R6*z2+99225*R8)*d19_2;
	an[9]=-(19958400*z9-179625600*R2*z7+314344800*R4*z5-130977000*R6*z3+9823275*R8*z)*d21_2;
if (McDOrder > 4){
	const double z10=z9*z;
	const double z11=z10*z;

	const double R10 = R8*R2;

	const double d23_2=d21_2*d;
	const double d25_2=d23_2*d;

	an[10]=(239500800*z10-2694384000*R2*z8+6286896000*R4*z6-3929310000*R6*z4+589396500*R8*z2-9823275*R10)*d23_2;
	an[11]=-(3113510400*z11-42810768000*R2*z9+128432304000*R4*z7-112378266000*R6*z5+28094566500*R8*z3-1404728325*R10*z)*d25_2;
if (McDOrder> 5){

	const double z12=z11*z;
	const double z13= z12*z;

	const double R12 = R10*R2;
			
	const double d27_2=d25_2*d;
	const double d29_2=d27_2*d;

	an[12]=(43589145600*z12-719220902400*R2*z10+2697078384000*R4*z8-3146591448000*R6*z6+1179971793000*R8*z4-117997179300*R10*z2+1404728325*R12)*d27_2;
	an[13]=-(653837184000*z13-12749825088000*R2*z11+58436698320000*R4*z9-87655047480000*R6*z7+46018899927000*R8*z5-7669816654500*R10*z3+273922023375*R12*z)*d29_2;
if (McDOrder>6){
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

// Shell

void Conway1D::getB(const double cylP[3], double BCylVec[3]) const {
	// From the "Exact Solution..." paper by J. T. Conway
	
	const double z = cylP[2];		
	//~ std::cout << "Z1 = " << Z1 << std::endl;
	//~ std::cout << "Z2 = " << Z2 << std::endl;
	//~ std::cout << "z = " << z << std::endl;
	
	if(z < Z1){
		BCylVec[2] = C * ( I_010(Z1-z, cylP) - I_010(Z2-z, cylP) );						 
		//~ std::cout << "I_010(Z1-z) =  " << I_010(Z1-z, cylP) << std::endl;
		//~ std::cout << "I_010(Z2-z) =  " << I_010(Z2-z, cylP) << std::endl;
	}else if(z >= Z1 && z <= Z2){
		BCylVec[2] = C * ( 2* I_010(0.0 , cylP)
							- I_010(Z1-z, cylP)
							- I_010(Z2-z, cylP) );							
		//~ std::cout << "I_010(0.0) =  " << I_010(0.0 , cylP) << std::endl;
		//~ std::cout << "I_010(Z1-z) = " << I_010(Z1-z, cylP) << std::endl;
		//~ std::cout << "I_010(Z2-z) = " << I_010(Z2-z, cylP) << std::endl;
	}else if(z > Z2){
		BCylVec[2] = C * ( I_010(Z2-z, cylP) - I_010(Z1-z, cylP) );						 
		//~ std::cout << "I_010(Z2-z) =  " << I_010(Z2-z, cylP) << std::endl;
		//~ std::cout << "I_010(Z1-z) =  " << I_010(Z1-z, cylP) << std::endl;
	}else{
		std::cout << "Error, cannot determine z in Conway1D" << std::endl;
	}
		
	BCylVec[0] = C * ( I_011(Z2-z, cylP) - I_011(Z1-z, cylP) );	
	//~ std::cout << "I_011(Z2-z) = " << I_011(Z2-z, cylP) << std::endl;
	//~ std::cout << "I_011(Z1-z) = " << I_011(Z1-z, cylP) << std::endl;
}

double Conway1D::I_010(const double zDiff, const double cylP[3]) const {
	// A Bessel-Laplace integral from the "Exact Solution..." paper by J. T. Conway
	// Radius of cylindrical shell, R
	// Axial distance between observation point and the centre of the cylindrical shell
	// observation point cylP = (r, phi, z)
	
	//~ std::cout << "zDiff = " << zDiff << std::endl;
	
	const double rho = cylP[0];
	//~ std::cout << "rho = " << zDiff << std::endl;	
	
	const double rSum = rho+R;
	const double rSum2 = rSum*rSum;
	const double rDiff = rho-R;
	const double rDiff2 = rDiff*rDiff;
	const double zDiff2 = zDiff*zDiff;
	
	const double k0 = 2.0/sqrt(rSum2 + zDiff2);
	const double k = sqrt(rho*R)*k0;
	const double beta = std::asin( zDiff / (sqrt( rDiff2 + zDiff2 )) );
	const double K_comp = std::tr1::comp_ellint_1(k); 
	
	//~ std::cout << "k 010 = " << k << std::endl;
	//~ std::cout << "beta 010 = " << beta << std::endl;
	//~ std::cout << "K_comp 010 = " << K_comp << std::endl;

	if(zDiff == 0 || beta == 0){
		//~ std::cout << "zDiff = 0 => beta = 0, 010" << std::endl;
		if(rho < R){
			//~ std::cout << "I_010 = " << 1.0/R << std::endl;
			return 1.0/R;
		}
		if(rho > R){
			return 0.0;
		}
	}

	if(rho == 0){
		//~ std::cout << "rho = 0, 010" << std::endl;
		if(rho < R){
			//~ std::cout << "I_010 = " << 1.0/R * (1 - fabs(zDiff)/(2.0*sqrt(R*R+zDiff2)) - std::tr1::ellint_2(1,fabs(beta))/2.0) << std::endl;
			return 1.0/R * (1 - fabs(zDiff)/(2.0*sqrt(R*R+zDiff2)) - std::tr1::ellint_2(1,fabs(beta))/2.0);
		}
		if(rho > R){
			return 1.0/R * (-fabs(zDiff)/(2.0*sqrt(R*R+zDiff2)) - std::tr1::ellint_2(1,fabs(beta))/2.0);
		}
	}
		
	if(rho < R){
		//~ std::cout << "I_010 = " << 1.0/R * ( 1 - fabs(zDiff)*k0*K_comp/(2*PhysicsConstants::pi) - HeumansLambda(fabs(beta),k)/2.0 ) << std::endl;
		return 1.0/R * ( 1 - fabs(zDiff)*k0*K_comp/(2*PhysicsConstants::pi) - HeumansLambda(fabs(beta),k)/2.0 );		
	}else if(rho > R){		
		return 1.0/R * ( -fabs(zDiff)*k0*K_comp/(2*PhysicsConstants::pi) + HeumansLambda(fabs(beta),k)/2.0 );		
	}else{		
		std::cout << "Error in I_010 function. Cannot evaluate integraol for rho=R" << std::endl;
		return 0;
	}
	
}

double Conway1D::I_011(const double z_src, const double cylP[3]) const {
	// A Bessel-Laplace integral from the "Exact Solution..." paper by J. T. Conway
	// Radius of cylindrical shell, R
	// Axial distance between observation point and the centre of the cylindrical shell
	// observation point cylP = (rho, phi, z)
	
	const double rho = cylP[0];
	
	if(rho==0){
		//~ std::cout << "rho = 0 => I_011 = 0" << std::endl; 
		return 0.0;
	}
	
	const double rSum = rho+R;
	const double rSum2 = rSum*rSum;
	const double zDiff = z_src;
	const double zDiff2 = zDiff*zDiff;
	
	const double k2 = 4*rho*R / (rSum2 + zDiff2);
	const double k = sqrt(k2);
	
	//~ std::cout << "k 011 = " << k << std::endl;
	
	const double E_comp = std::tr1::comp_ellint_2(k); // complete elliptic integral of second kind
	const double K_comp = std::tr1::comp_ellint_1(k); // complete elliptic integral of first kind
	
	//~ std::cout << "E_comp 011 = " << E_comp << std::endl;
	//~ std::cout << "K_comp 011 = " << K_comp << std::endl;
	
	//~ std::cout << "pi*k*sqrt(rho*R) = " << PhysicsConstants::pi*k*sqrt(rho*R)<< std::endl;
	
	return 1.0/(PhysicsConstants::pi*k*sqrt(rho*R)) * ( (2-k2)*K_comp - 2*E_comp );
}

double Conway1D::HeumansLambda(const double beta, const double k) const {
	// Heumans Lambda function. Used in the "Exact Solution..." paper by J. T. Conway
	
	if(beta == 0){
		return 0.0;
	}	
	
	const double kPrime = sqrt(1.0 - k*k);
	const double E_comp = std::tr1::comp_ellint_2(k); // complete elliptic integral of second kind
	const double K_comp = std::tr1::comp_ellint_1(k); // complete elliptic integral of first kind
	const double E = std::tr1::ellint_2(kPrime,beta); // elliptic integral of second kind
	const double F = std::tr1::ellint_1(kPrime,beta); // elliptic integral of first kind
	
	if(k == 0){
		return E;
	}
	
	//~ std::cout << "k = " << k << std::endl;
	//~ std::cout << "beta = " << beta << std::endl;
	//~ std::cout << "E = " << E << std::endl;
	//~ std::cout << "F = " << F << std::endl;
	
	return 2.0/PhysicsConstants::pi * (E_comp*F + K_comp*E - K_comp*F);
}

void McD_Shell::getB(const double cylP[3], double BCylVec[3]) const {
/* Calculates the magnet field in BCylVec at the cylindrical coordinate cylVec for a shell solenoid of radius R and total current I	using the method described by the McDonald model
 * 
 * @param shell the shell creating the magnetic field
 * @param McDOrder an integer [0,7] that the denotes the order of the McDonald calculation
 * @param cylP the cylindrical coordinate of interest where the magnetfield is calculated
 * @param BCylVec the magnetfield in cylP being calulated in cylindrical coordinates
 */
	
	// Coordinates
	const double rho = cylP[0];
	const double z1 = cylP[2] - Z1;
	const double z2 = cylP[2] - Z2;
	
	double an1[2*McDOrder+2];
	double an2[2*McDOrder+2];
	mcDonaldShellSupFunc(z1,an1);
	mcDonaldShellSupFunc(z2,an2);
	
	// Preparing terms for loop
	double B_z = an1[0]-an2[0]; 
	double constZTerm = 1;
	
	double B_rho = -(rho/2)*(an1[1]-an2[1]);
	double constRTerm = -(rho/2);
	
	// Looping over the series
	for (int n=1; n<=McDOrder; n++){ 
		constZTerm *= (-1)*pow(rho/2,2)/pow(n,2);
		B_z+= constZTerm*(an1[2*n]-an2[2*n]);
		constRTerm *= (-1)*pow(rho/2,2)*(n)/((n+1)*pow(n,2));
		B_rho += constRTerm*(an1[2*n+1]-an2[2*n+1]);
	}
	
	// Preparing result	
	B_rho *= C;
	B_z *= C;
	
	// Result
	BCylVec[0]=B_rho;
	BCylVec[1]=0.0;
	BCylVec[2]=B_z;
}

void McD_Shell::mcDonaldShellSupFunc(const double z, double an[]) const {	
/* This function is a support function to loopApproxMcDonald it generates  the an[] list depending on the order of n to reduce calculation time.
 * 
 * @param n		the order of the McDonald model required
 * @param z		the z-coordinate where the field is being calculated
 * @param R 	the radius of the loop
 * @param an	the list of derivatives of the axial magnetic field
 */
	
	// Powers of z
	const double z2=z*z;  
		
	// Powers of d
	const double d = 1/(z2+R2);
	const double d1_2 = sqrt(d); 
	const double d3_2 = R2*d1_2*d;

	an[0]=z*d1_2;
	an[1]=d3_2;
if (McDOrder > 0){ 
	const double z2 = z*z;
	
	const double d5_2=d3_2*d;
	const double d7_2=d5_2*d;
	an[2]=-3*z*d5_2;
	an[3]=(12*z2-3*R2)*d7_2;
		
if (McDOrder > 1){
	const double z3=z2*z;
	const double z4=z3*z;

	const double R4 = R2*R2;
				
	const double d9_2=d7_2*d;
	const double d11_2=d9_2*d;
	an[4]=(-60*z3+45*R2*z)*d9_2;			
	an[5]=(360*z4-540*R2*z2+45*R4)*d11_2;

if (McDOrder > 2){
	const double z5=z4*z;
	const double z6=z5*z;

	const double R6 = R4*R2;
		
	const double d13_2=d11_2*d;
	const double d15_2=d13_2*d;

	an[6]=(-2520*z5+6300*R2*z3-1575*R4*z)*d13_2;
	an[7]=(20160*z6-75600*R2*z4+37800*R4*z2-1575*R6)*d15_2;

if (McDOrder > 3){
	const double z7=z6*z;
	const double z8=z7*z;
	
	const double R8 = R6*R2;
		
	const double d17_2=d15_2*d;
	const double d19_2=d17_2*d;
	
	an[8]=(-181440*z7+952560*R2*z5-793800*R4*z3+99225*R6*z)*d17_2;
	an[9]=(1814400*z8-12700800*R2*z6+15876000*R4*z4-3969000*R6*z2+99225*R8)*d19_2;

if (McDOrder > 4){
	const double z9=z8*z;
	const double z10=z9*z;

	const double R10 = R8*R2;

	const double d21_2=d19_2*d;
	const double d23_2=d21_2*d;
	
	an[10]=-(19958400*z9-179625600*R2*z7+314344800*R4*z5-130977000*R6*z3+9823275*R8*z)*d21_2;
	an[11]=(239500800*z10-2694384000*R2*z8+6286896000*R4*z6-3929310000*R6*z4+589396500*R8*z2-9823275*R10)*d23_2;
if (McDOrder> 5){
	const double z11=z10*z;
	const double z12=z11*z;

	const double R12 = R10*R2;
			
	const double d25_2 = d23_2*d;	
	const double d27_2 = d25_2*d;
	
	an[12]=-(3113510400*z11-42810768000*R2*z9+128432304000*R4*z7-112378266000*R6*z5+28094566500*R8*z3-1404728325*R10*z)*d25_2;
	an[13]=(43589145600*z12-719220902400*R2*z10+2697078384000*R4*z8-3146591448000*R6*z6+1179971793000*R8*z4-117997179300*R10*z2+1404728325*R12)*d27_2;

if (McDOrder>6){
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

void NWire_Shell::getB(const double cylP[3], double BCylVec[3]) const {
	// this model represents the magnet with an NxM wire grid, and calculates the total field
	// as the sum of the individual wire contributions. The field of a wire is calculated
	// with the SAM or McD method (see what method is commented in)	
	
	// make sure that the placeholder for the magnetic field is 0
	BCylVec[0] = 0;
	BCylVec[1] = 0;
	BCylVec[2] = 0;

    double BCylVec_i[3];
    
	for(int n_z = 0; n_z < N_z; n_z++){
		loops[n_z].getB(cylP,BCylVec_i);
		
		BCylVec[0] += BCylVec_i[0];
		BCylVec[1] += BCylVec_i[1];
		BCylVec[2] += BCylVec_i[2];		
	}
}

void GaussianQuadratureLoops_Shell::getB(const double cylP[3], double BCylVec[3]) const {	
    // make sure that the placeholder for the magnetic field is 0
	BCylVec[0] = 0;
	BCylVec[1] = 0;
	BCylVec[2] = 0;
 
    double BCylVec_i[3];
    
	for(int nGP_z = 0; nGP_z < NGP_z; nGP_z++){
		loops[nGP_z].getB(cylP,BCylVec_i);		
		//~ printVec(BCylVec_i,"BCyl_i");
		
		const double GFac = GPZWeights[nGP_z];
		
		BCylVec[0] += BCylVec_i[0]*GFac;
		BCylVec[1] += BCylVec_i[1]*GFac;
		BCylVec[2] += BCylVec_i[2]*GFac;
	}
}

// Tube

void McD_Tube::getB(const double cylP[3], double BCylVec[3]) const{	
	/* Calculates the magnet field in BCylVec at the cylindrical coordinate cylVec for a finite solenoid of radius R and total current I using the method described by the McDonald model
	 * 
	 * @param tube 			the tube creating the magnetic field
	 * @param cylP 			the cylindrical coordinate of interest where the magnetfield is calculated
	 * @param BCylVec 		the magnetfield in cylP being calulated in cylindrical coordinates
	 * @param McDSupport	this function uses a class to precalculate the derivatives of the McD model. This is FASTER AND BETTER
	 */
	
	// Coordinates
	const double rho = cylP[0];
	const double z1 = cylP[2] - Z1; // define coordinates relative to the position of the solenoid
	const double z2 = cylP[2] - Z2;	
	
	double an11[Na] = {0};		// array to hold terms with z1 and R1
	double an12[Na] = {0};		// array to hold terms with z1 and R2
	double an21[Na] = {0};		// array to hold terms with z2 and R1
	double an22[Na] = {0};		// array to hold terms with z2 and R2
	
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
	double constZTerm = 1;	
	double B_rho = -(rho/2)*(an12[1]-an11[1]-an22[1]+an21[1]);
	double constRTerm = -(rho/2);
	
	// Looping over the series
	for (int n=1; n<=McDOrder; n++){ 
		constZTerm *= (-1)*pow(rho/2,2)/pow(n,2);
		B_z+= constZTerm*(an12[2*n]-an11[2*n]-an22[2*n]+an21[2*n]);
		constRTerm *= (-1)*pow(rho/2,2)*(n)/((n+1)*pow(n,2));
		B_rho += constRTerm*(an12[2*n+1]-an11[2*n+1]-an22[2*n+1]+an21[2*n+1]);
	}
	
	// Preparing result
	B_rho *= C;
	B_z *= C;
	
	// Result
	BCylVec[0]=B_rho;
	BCylVec[1]=0.0;
	BCylVec[2]=B_z;
}

void McD_Tube::getA0(const double z, const double R, double a0_array[]) const{
	// This function calculates the the a_n terms for a given z and R. 
	// The derivatives of a_0 have been precalculated by the class constructor, which
	// was instantiated to be able to go to a given order. The requested order has to
	// be below this
	// various outcommented debugging code has been left here

	//~ std::cout << "number of a_n terms to calculate, Na = " << Na <<std::endl;
	//~ for(int i=0; i<Na; i++){
		//~ a0_array[i] = 0; //make sure all elements are 0
	//~ }	
	
	const double z2 = z*z;
	const double R2 = R*R;
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
	
	if(Na>1){
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

void TAVP::getB(const double sphP[3], double BSphVec[3]) const {	
	// an approximation of the B-field of a current loop
	// based on eq A.2 in https://iopscience.iop.org/article/10.1088/1367-2630/14/1/015010
	// input: radius R, mirror z position z, current I, model parameter lambda,
	// observation point r, placeholder for B in cartesian coor B_vec
	
	//~ printVec(sphP,"sphP");
	
	// Converting the spherical point to cartesian coordinates, and makes it realtive to the centre of the magnet
	const double rsintheta = sphP[0]*sin(sphP[1]);
	const double carP[3] = {rsintheta*cos(sphP[2]), rsintheta*sin(sphP[2]), sphP[0]*cos(sphP[1]) - Loop::getz()};
	// Convert back to sphericla coordinates. The transformation has not affected phi = sphP[2]
	const double s2 = carP[0]*carP[0] + carP[1]*carP[1]; 
	const double r = std::sqrt(s2 + carP[2]*carP[2]);
	double theta;
	if (s2 == 0 && carP[2] == 0){
		theta = 0;
	}else if (s2 == 0 && carP[2] > 0){
		theta = 0;
	}else if (s2 == 0 && carP[2] < 0){
		theta = PhysicsConstants::pi;
	}else if (carP[2] == 0){
		theta = PhysicsConstants::pi/2.0;
	}else{
		theta = atan2(sqrt(s2),carP[2]); // arctan(sqrt(y^2 + x^2)/z)
	};


	double B_r;
	double B_theta;
	
	const double r2 = r*r;
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
		B_r = PhysicsConstants::mu0*I/(2*R);		
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

void Helix::getB(const double carP[3], double BCarVec[3]) const {
	// make sure that the placeholder for the magnetic field is 0
	BCarVec[0] = 0;
	BCarVec[1] = 0;
	BCarVec[2] = 0;		
	
	//~ std::cout << "tube centre = " << z << std::endl;
	//~ std::cout << "R1 = " << R1 << std::endl;
	//~ std::cout << "R2 = " << outerRadius << std::endl;
	//~ std::cout << "thickness = " << width_rho << std::endl;
	//~ std::cout << "L = " << width_z << std::endl;
	//~ std::cout << "delta_rho = " << delta_rho << std::endl;
	//~ std::cout << "delta_z = " << delta_z << std::endl;
	//~ std::cout << "N_rho = " << N_rho << std::endl;
	//~ std::cout << "N_z = " << N_z << std::endl;
	//~ std::cout << "N_BS = " << N_BS << std::endl;
	//~ std::cout << "I = " << i << std::endl;	
	
	double BCarVec_i[3]{0.,0.,0.};
	
	for(int n=0; n<N_rho*N_z*N_BS; n++){
		segments[n].getB(carP,i,BCarVec_i);
		//~ std::cout << "i = " << n_rho*N_z*N_BS + n_z*N_BS + n_BS << std::endl;
			
		BCarVec[0] += BCarVec_i[0];
		BCarVec[1] += BCarVec_i[1];
		BCarVec[2] += BCarVec_i[2];	
		
	}
}

void NWire_Tube::getB(const double cylP[3], double BCylVec[3]) const {
	// this model represents the magnet with an NxM wire grid, and calculates the total field
	// as the sum of the individual wire contributions. The field of a wire is calculated
	// with the SAM or McD method (see what method is commented in)	
	
	// make sure that the placeholder for the magnetic field is 0
	BCylVec[0] = 0;
	BCylVec[1] = 0;
	BCylVec[2] = 0;

    double BCylVec_i[3];
    
    for(int n=0; n<N_rho*N_z; n++){
		loops[n].getB(cylP, BCylVec_i);
		
		BCylVec[0] += BCylVec_i[0];
		BCylVec[1] += BCylVec_i[1];
		BCylVec[2] += BCylVec_i[2];	
	}
}

void GaussianQuadratureLoops_Tube::getB(const double cylP[3], double BCylVec[3]) const {	
    // make sure that the placeholder for the magnetic field is 0
	BCylVec[0] = 0;
	BCylVec[1] = 0;
	BCylVec[2] = 0;
    
    double BCylVec_i[3];
    
	for(int n=0; n<NGP_rho*NGP_z; n++){
		loops[n].getB(cylP, BCylVec_i);
		
		const double GFac = GFacs[n];
		
		BCylVec[0] += BCylVec_i[0]*GFac;
		BCylVec[1] += BCylVec_i[1]*GFac;
		BCylVec[2] += BCylVec_i[2]*GFac;	
	}
}

void GaussianQuadratureShells_Tube::getB(const double cylP[3], double BCylVec[3]) const {	
    // make sure that the placeholder for the magnetic field is 0
	BCylVec[0] = 0;
	BCylVec[1] = 0;
	BCylVec[2] = 0;
    
    double BCylVec_i[3] = {0.,0.,0.};
    
	for(int nGP_rho = 0; nGP_rho < NGP_rho; nGP_rho++){		
		shells[nGP_rho].getB(cylP, BCylVec_i);
				
		const double GFac = GPRhoWeights[nGP_rho];
			
		BCylVec[0] += BCylVec_i[0]*GFac;
		BCylVec[1] += BCylVec_i[1]*GFac;
		BCylVec[2] += BCylVec_i[2]*GFac;	
	}
}
