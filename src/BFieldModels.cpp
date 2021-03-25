//~ #define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1 // needed for elliptic integral functions
#include "BFieldModels.h"

// Basic functions

void biotSavart(const double s0[3], const double s1[3], const double r[3], const double I, double B_vec[3]){

	//~ printVec(s0,"s0");
	//~ printVec(s1,"s1");
	//~ printVec(r,"r");
	//~ printVec(B_vec,"B in");
	
	double i[3]={0.}; double k0[3]={0.}; double k1[3]={0.};
	vecSub(s1,s0,i);
	vecSub(s0,r,k0);
	vecSub(s1,r,k1);
	
	//~ double i_norm; double k0_norm; double k1_norm;
	const double i_norm = vecNorm(i);
	const double k0_norm = vecNorm(k0);
	const double k1_norm = vecNorm(k1);
	
	
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
		
	const double C = PhysicsConstants::mu0/(4*PhysicsConstants::pi)*I*a0;
	vecMultScalOvrwrt(A,C);
	
	//~ std::cout << "C = " << C << std::endl;
	//~ printVec(A,"A*C");
	
	vecAddOvrwrt(B_vec,A);
	
	//~ printVec(B_vec,"B out");
}

// Current Loop

void loopExactSAM(const Loop& loop, const double cylP[3], double BCylVec[3]){
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
	const double C = PhysicsConstants::mu0*loop.getI()/PhysicsConstants::pi;
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

void loopBiotSavart(const Loop& loop, const int NSegments, const double carP[3], double BCarVec[3]){
	// calculates the B-field of a current loop in the x-y plane by using the BS method
	// radius R, z position z, current I, number of segments N, obs point carP

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
	double phis[NSegments];
	double dPhi = 2*PhysicsConstants::pi/(NSegments-1);
	for(int i=0; i<NSegments; i++){
		phis[i] = dPhi*i;
	}
	double x0 = cos(phis[0])*R;
	double y0 = sin(phis[0])*R;
	double x1;
	double y1;
	
	for(int i=0; i<NSegments-1; i++){
		x1 = cos(phis[i+1])*R;
		y1 = sin(phis[i+1])*R;
		double r0[3]{x0,y0,z_loop};
		double r1[3]{x1,y1,z_loop};
			
		biotSavart(r0,r1,carP,I,BCarVec);
		
		x0 = x1;
		y0 = y1;
	}	
}

void mcDonald(const Loop& loop, const int nmax, const double cylP[3], double BCylVec[3]){	
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
	const double C = PhysicsConstants::mu0*I*R*R/2;
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


// Shell

void Conway1D(const Shell& shell, const double cylP[3], double BCylVec[3]){
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
	const double C = PhysicsConstants::mu0*IDist*R/2.0;	
	
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

void mcDonald(const Shell& shell, const int nmax, const double cylP[3], double BCylVec[3]){	
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
	const double C = PhysicsConstants::mu0*I/(2*L);
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

void NWire(const Shell& shell, const int N_z, const double cylP[3], double BCylVec[3]){
	// this model represents the magnet with an NxM wire grid, and calculates the total field
	// as the sum of the individual wire contributions. The field of a wire is calculated
	// with the SAM or McD method (see what method is commented in)	
		
    // Tube Center
    const double x = shell.getx();
    const double y = shell.gety();
    const double z = shell.getz();

	// Radii
	const double radius = shell.getR();

    // Length
	const double width_z = shell.getL(); 	

    // Spacing
	const double delta_z = width_z/N_z; 

    // Current			
	const double I = shell.getI();					
	
    double BCylVec_i[3];
	for(int j = 0; j < N_z; j++){
		const Loop loop = Loop(radius, I/N_z, x, y, z - 0.5*width_z + delta_z/2.0 + j*delta_z);
				
		loopExactSAM(loop, cylP, BCylVec_i);
		BCylVec[0] += BCylVec_i[0];
		BCylVec[1] += BCylVec_i[1];
		BCylVec[2] += BCylVec_i[2];		
	}
}

void GaussianQuadratureLoop(const Shell& shell, const int N_z, const int NG_z, const double cylP[3], double BCylVec[3]){	
    // Tube Center
    const double x = shell.getx();
    const double y = shell.gety();
    const double z = shell.getz();

	// Radii
	const double radius = shell.getR();					
	
    // Dimension of Tube
	const double width_z = shell.getL();

    // Current			
	const double I = shell.getI()/(N_z);	// divide the current of the shell on N_z loops

	double GPZValues[NG_z];

	double GPZWeights[NG_z];

	const int NWiresZ = N_z/2;
	
	getGaussianQuadratureParams(NG_z,GPZValues,GPZWeights,width_z/2.0,NWiresZ);
	
    double BCylVec_i[3];
	for(int nGP_z = 0; nGP_z < NG_z; nGP_z++){
		const Loop loop(radius,I,x,y,z+GPZValues[nGP_z]);		
		loopExactSAM(loop, cylP, BCylVec_i);
		//~ printVec(BCylVec_i,"BCyl_i");
		
		const double GFac = GPZWeights[nGP_z];
		
		BCylVec[0] += BCylVec_i[0]*GFac;
		BCylVec[1] += BCylVec_i[1]*GFac;
		BCylVec[2] += BCylVec_i[2]*GFac;
	}
	
	//~ printVec(BCylVec,"BCyl");
}

// Tube

void mcDonald(const Tube& tube, const int nmax, const double cylP[3], double BCylVec[3], const McD_Tube_Support& McDSupport){	
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
	const double rho = cylP[0];
	const double z   = cylP[2]-tube.getz();  	// Posistion the tube in the center
	const double R1  = tube.getR1();			// inner radius
	const double R2  = tube.getR2();			// outer radius
	const double I  = tube.getI();			// current
	const double L  = tube.getL();			// length of solenoid
	const double L_2 = 0.5*L;
	
	const double z1 = z + L_2;	// The distance from the lower point of the tube to the point of evaluation
	const double z2 = z - L_2;	// The distance from the upper point of the tube to the point of evaluation
	
	const int arraySize = 2*nmax + 2; // the number of coefficients needed to calculate the requested order (nmax)
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
	const double C = PhysicsConstants::mu0*I/(2*L*(R2-R1));
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

void loopApproxFrancis(const Loop& loop, const double lambda, const double sphP[3], double BSphVec[3]){
	
	// an approximation of the B-field of a current loop
	// based on eq A.2 in https://iopscience.iop.org/article/10.1088/1367-2630/14/1/015010
	// input: radius R, mirror z position z, current I, model parameter lambda,
	// observation point r, placeholder for B in cartesian coor B_vec
	
	const double R = loop.getR();
	const double r = sphP[0];
	const double theta = sphP[1];

	// Constants
	const double C = 0.125*PhysicsConstants::mu0*loop.getI()*R/lambda;
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
		B_r = PhysicsConstants::mu0*loop.getI()/(2*R);		
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

void Helix(const Tube& tube, const int N_rho, const int N_z, const int N_BS, const double carP[3], double BCarVec[3]){
		
    // Tube Center
    const double z = tube.getz();

	// Radii
	const double innerRadius = tube.getR1();			
	const double outerRadius = tube.getR2();				
	const double width_rho = outerRadius - innerRadius;
	
    // Length
	const double width_z = tube.getL(); 	

    // Spacing			
	const double delta_rho = width_rho/N_rho; 	
	const double delta_z = width_z/N_z; 

    // Current			
	const double i = tube.geti();				
	
	//~ std::cout << "tube centre = " << z << std::endl;
	//~ std::cout << "R1 = " << innerRadius << std::endl;
	//~ std::cout << "R2 = " << outerRadius << std::endl;
	//~ std::cout << "thickness = " << width_rho << std::endl;
	//~ std::cout << "L = " << width_z << std::endl;
	//~ std::cout << "delta_rho = " << delta_rho << std::endl;
	//~ std::cout << "delta_z = " << delta_z << std::endl;
	//~ std::cout << "N_rho = " << N_rho << std::endl;
	//~ std::cout << "N_z = " << N_z << std::endl;
	//~ std::cout << "N_BS = " << N_BS << std::endl;
	//~ std::cout << "I = " << i << std::endl;
	

	for(int n_rho = 0; n_rho < N_rho; n_rho++){				// loop over all layers
		double cylS0[3];
		double cylS1[3]{innerRadius + delta_rho/2.0 + n_rho*delta_rho, 0.0, z + pow(-1,n_rho)*(-width_z/2.0) }; // initial value for cylS1 (to copy into cylS0 in the loop) 
		cylS0[0] = cylS1[0]; // for a given layer, the rho-coordinate stays the same for all segments
	
		for(int n_z = 0; n_z < N_z; n_z++){				// loop over all windings in a layer
			for(int n_BS = 0; n_BS < N_BS; n_BS++){		// loop over all straight line segments in a winding
				
                double BCarVec_i[3]{0.,0.,0.};				
				cylS0[1] = cylS1[1]; 
				cylS0[2] = cylS1[2]; // the last end point of a segment is now the start point
				cylS1[1] = 2.0*PhysicsConstants::pi/N_BS*(n_BS+1);
				cylS1[2] = z + pow(-1,n_rho)*(-width_z/2.0 + delta_z/N_BS*(n_BS+1) + n_z*delta_z);
				
				double carS0[3], carS1[3];
				cylPToCarP(cylS0,carS0); // convert the start and end point to cartesian coordinates
				cylPToCarP(cylS1,carS1);				

				//~ if(n_z == 0 && n_BS == 0){
				//~ std::cout << "BS params:" << std::endl;
				//~ std::cout << "n_rho = " << n_rho << " n_z = " << n_z << " n_BS = " << n_BS << std::endl;
				//~ printVec(carS0,"s0");
				//~ printVec(carS1,"s1");
				//~ printVec(carP,"r_i");
				//~ std::cout << "I = " << i << std::endl;
				//~ printVec(BCarVec_i,"B_i");
				//~ }

				biotSavart(carS0,carS1,carP,i,BCarVec_i);	
				BCarVec[0] += BCarVec_i[0];
				BCarVec[1] += BCarVec_i[1];
				BCarVec[2] += BCarVec_i[2];					
							
			}
		}
	}
}

void NWire(const Tube& tube, const int N_rho, const int N_z, const double cylP[3], double BCylVec[3]){
	// this model represents the magnet with an NxM wire grid, and calculates the total field
	// as the sum of the individual wire contributions. The field of a wire is calculated
	// with the SAM or McD method (see what method is commented in)	
		
    // Tube Center
    const double x = tube.getx();
    const double y = tube.gety();
    const double z = tube.getz();

	// Radii
	const double innerRadius = tube.getR1();			
	const double outerRadius = tube.getR2();				
	const double width_rho = outerRadius - innerRadius;
	
    // Number of wires
	const int N_wires = N_rho*N_z;

    // Length
	const double width_z = tube.getL(); 	

    // Spacing			
	const double delta_rho = width_rho/N_rho; 	
	const double delta_z = width_z/N_z; 

    // Current			
	const double I = tube.getI();					
	
    double BCylVec_i[3];
	for(int i = 0; i < N_rho; i++){
		for(int j = 0; j < N_z; j++){
			//~ std::cout << " i, j = " << i << ", " << j << std::endl;
			Loop loop = Loop(innerRadius + delta_rho/2.0 + i*delta_rho, I/N_wires, x, y, z - 0.5*width_z + delta_z/2.0 + j*delta_z);
					
			loopExactSAM(loop, cylP, BCylVec_i);

			BCylVec[0] += BCylVec_i[0];
			BCylVec[1] += BCylVec_i[1];
			BCylVec[2] += BCylVec_i[2];		
		}
	}
}

void GaussianQuadratureLoop(const Tube& tube, const int N_rho, const int N_z, const int NG_rho, const int NG_z, const double cylP[3], double BCylVec[3]){	
    // Tube Center
    const double x = tube.getx();
    const double y = tube.gety();
    const double z = tube.getz();

	// Radii
	const double innerRadius = tube.getR1();			
	const double outerRadius = tube.getR2();				
	
    // Dimension of Tube
	const double width_z = tube.getL(); 	
    const double width_rho = outerRadius - innerRadius;

    // Current			
	const double I = tube.getI()/(N_rho*N_z);	// divide the current of the tube on N_rho*N_z loops

	
	double GPRhoValues[NG_rho];
	double GPZValues[NG_z];
	
	double GPRhoWeights[NG_rho];
	double GPZWeights[NG_z];
	
	const int NWiresRho = N_rho/2; // half the number of wires in each dimension
	const int NWiresZ = N_z/2;
	
	getGaussianQuadratureParams(NG_rho,GPRhoValues,GPRhoWeights,width_rho/2.0,NWiresRho);
	getGaussianQuadratureParams(NG_z,GPZValues,GPZWeights,width_z/2.0,NWiresZ);
	
	const double centre_rho = innerRadius + width_rho/2.0;
	
    double BCylVec_i[3];
	for(int nG_rho = 0; nG_rho < NG_rho; nG_rho++){
		for(int nGP_z = 0; nGP_z < NG_z; nGP_z++){
			const Loop loop(centre_rho + GPRhoValues[nG_rho],I,x,y,z+GPZValues[nGP_z]);		
			loopExactSAM(loop, cylP, BCylVec_i);
			//~ printVec(BCylVec_i,"BCyl_i");
			
			const double GFac = GPZWeights[nGP_z]*GPRhoWeights[nG_rho];
			
            BCylVec[0] += BCylVec_i[0]*GFac;
			BCylVec[1] += BCylVec_i[1]*GFac;
			BCylVec[2] += BCylVec_i[2]*GFac;
		}
	}
	
	//~ printVec(BCylVec,"BCyl");
}

void GaussianQuadratureShell(const Tube& tube, const int N_rho, const int NG_rho, const double cylP[3], double BCylVec[3]){	
    // Tube Center
    const double x = tube.getx();
    const double y = tube.gety();
    const double z = tube.getz();

	// Radii
	const double innerRadius = tube.getR1();			
	const double outerRadius = tube.getR2();				
	
    // Dimension of Tube
	const double width_z = tube.getL(); 	
    const double width_rho = outerRadius - innerRadius;

    // Current			
	const double I = tube.getI()/N_rho; // devide the current of the tube on N_rho shells

	double GPRhoValues[NG_rho];
	double GPRhoWeights[NG_rho];

	const int NElementsRho = N_rho/2; // half the number of layers in each dimension	
	getGaussianQuadratureParams(NG_rho,GPRhoValues,GPRhoWeights,width_rho/2.0,NElementsRho);
	
	const double centre_rho = innerRadius + width_rho/2.0;

    double BCylVec_i[3];
	for(int nG_rho = 0; nG_rho < NG_rho; nG_rho++){
		//~ std::cout << "nG_Rho = " << nG_rho << std::endl;
		const Shell shell = Shell(centre_rho + GPRhoValues[nG_rho],I,width_z,x,y,z);
		
		Conway1D(shell, cylP, BCylVec_i);
		
		const double GFac = GPRhoWeights[nG_rho];
			
		BCylVec[0] += BCylVec_i[0]*GFac;
		BCylVec[1] += BCylVec_i[1]*GFac;
		BCylVec[2] += BCylVec_i[2]*GFac;
	
	}
	
	//~ printVec(BCylVec,"BCyl");
}

// Utils

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
