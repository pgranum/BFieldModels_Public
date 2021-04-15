#ifndef BFIELDMODELS_H
#define BFIELDMODELS_H

#include <iostream>
#include <cmath>
#include <math.h>
#include <tr1/cmath>
#include <tgmath.h>
#include <chrono>
#include <cassert>

#include "Loop.h"
#include "Shell.h"
#include "Tube.h"
#include "utils.h"
#include "PhysicsConstants.hpp"
#include "mathFuncs.h"

// UTILS

class BiotSavart {
	private:
	
	public:
	BiotSavart(){} // end of constructor
	
	~BiotSavart(){} // end of destructor
	
	void BSSegment(const double s0[3], const double s1[3], const double r[3], const double I, double B_vec[3]) const ;
};

class GQ_Support {
	private:
	
	public:
	GQ_Support(){
	
	} // end of constructor
	
	~GQ_Support(){
	} // end of destructor
	
	//~ void getGaussianQuadratureParams(const int N, double points[], double weights[], const double dimLength, const double NWires) const ;
	void getGaussianQuadratureParams(const int N, std::vector<double> &points, std::vector<double> &weights, const double dimLength, const double NWires) const ;
};

// LOOP

class SimpleAnalyticModel : public Loop{
	private:
	double R2;
	
	public:
	SimpleAnalyticModel(const double R, const double I, const double x, const double y, const double z) : Loop(R,I,x,y,z){
		R2 = R*R;
	} // end of constructor
	
	~SimpleAnalyticModel(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ;
};

class BiotSavart_Loop : public Loop, public BiotSavart{
	private:
	const int N_BS;
	std::vector<double> xs;
	std::vector<double> ys;
	
	public:
	BiotSavart_Loop(const int N_BS, const double R, const double I, const double x, const double y, const double z) : Loop(R,I,x,y,z), N_BS(N_BS){
		// setting up the straght line segments used by the Biot-Savart method
		const double dPhi = 2*PhysicsConstants::pi/(N_BS-1);
		for(int i=0; i<N_BS; i++){
			const double phi = i*dPhi;
			xs.push_back(cos(phi)*R);
			ys.push_back(sin(phi)*R);
		}
	} // end of constructor
	
	~BiotSavart_Loop(){
	} // end of destructor
	
	void getB(const double carP[3], double BCarVec[3]) const ;
};

class McD_Loop : public Loop{
	private:
	const int McDOrder;
	const double C;
	
	public:
	McD_Loop(const int McDOrder, const double R, const double I, const double x, const double y, const double z) : 
		Loop(R,I,x,y,z), McDOrder(McDOrder),
		C(PhysicsConstants::mu0*I*R*R*0.5)
	{
	} // end of constructor
	
	~McD_Loop(){
	} // end of destructor
	
	void getB(const double carP[3], double BCarVec[3]) const ;
	void mcDonaldLoopSupFunc(const int n, const double z, const double R, double an[]) const ;
};

// SHELL

class Conway1D : public Shell{
	private:
	const double IDens;
	const double C;
	
	public:
	Conway1D(const double R, const int N, const double i, const double L, const double x, const double y, const double z) : 
	Shell(R,N,i,L,x,y,z),
	IDens(i*N/L),
	C(PhysicsConstants::mu0*IDens*R*0.5)
	{
		// R	shell radius
		// N 	number of wire turns in the shell/number of loops representing the shell
		// i 	current in wire/current per loop
		// L	length of the shell
		// x	x position of shell centre
		// y	y position of shell centre
		// z	z position of shell centre
	
	} // end of constructor
	
	~Conway1D(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ;
};

class McD_Shell : public Shell {
	private:
	const int McDOrder;
	const double L_2 = 0.5*this->getL();
	const double C;
	
	public:
	McD_Shell(const int McDOrder, const double R, const int N, const double i, const double L, 
			  const double x, const double y, const double z): 
	Shell(R,N,i,L,x,y,z), 
	McDOrder(McDOrder),
	C(PhysicsConstants::mu0*i*N/(2*L))
	{
	} // end of constructor
	
	~McD_Shell(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ; 
	void mcDonaldShellSupFunc(const int n, const double z, double an[]) const ;
};

class NWire_Shell : public Shell{
	private:
	const int N_z; // the number of wires in the shell
	std::vector<SimpleAnalyticModel> loops;
	
	public:
	NWire_Shell(const int N_z, const double R, const int N, const double i, const double L, 
				const double x, const double y, const double z) : Shell(R,N,i,L,x,y,z), N_z(N_z){
		
		// Spacing
		const double delta_z = L/N_z; // the spcing between the individual loops
		const double I = i*N/N_z;
		
		for(int n_z = 0; n_z < N_z; n_z++){
			loops.push_back(SimpleAnalyticModel(R,I,x,y,z - 0.5*L + delta_z*0.5 + n_z*delta_z));
		}
	} // end of constructor
	
	~NWire_Shell(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ;
};

class GaussianQuadratureLoops_Shell : public Shell, public GQ_Support{
	private:
	const int N_z; // the number of wires making up the Shell
	const int NGP_z; // the number of wires used to represent the shell in the Gaussian Quadrature
	std::vector<double> GPZValues;
	std::vector<double> GPZWeights;
	std::vector<SimpleAnalyticModel> loops;
	
	public:
	GaussianQuadratureLoops_Shell(const int N_z, const int NGP_z, const double R, const int N, const double i, const double L,
								  const double x, const double y, const double z) : Shell(R,N,i,L,x,y,z), N_z(N_z), NGP_z(NGP_z){
		assert(N_z % 2 == 0); // N_z has to be even, as the GQ alogorithm assumes the same number of wires in each half of the magnet
		
		const double width_z = this->getL(); 
		
		getGaussianQuadratureParams(NGP_z,GPZValues,GPZWeights,width_z*0.5,N_z/2);
		
		const double I = i*N/(N_z);	// divide the current of the shell on N_z loops
		
		for(int nGP_z = 0; nGP_z < NGP_z; nGP_z++){
			loops.push_back(SimpleAnalyticModel(R,I,x,y,z+GPZValues[nGP_z]));
		}
	} // end of constructor
	
	~GaussianQuadratureLoops_Shell(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ;
};

// TUBE

class McD_Tube : public Tube{
	private:
	const int McDOrder;
	const int Na; 	// the number of derivatives needed to calculate the requested order (McDOrder)
	const double L_2;
	const double C;
	
	std::vector<int> k_arr = {3, -1, -1};
	std::vector<int> lambda_arr = {1, 3, 3};
	std::vector<int> mu_arr = {1, 2, 3};
	std::vector<int> nu_arr = {1, 2, 1};
	std::vector<int> Ti_arr = {2, 2, 2};
	int NT; // the total number of terms in a_n
	
	public:
	int GetMcDOrder() const {
        return McDOrder;
	}
	
	McD_Tube(const int McDOrder, const double R1, const double R2, 
			 const int N, const double i, const double L, const double x, const double y, const double z) : 
	Tube(R1,R2,N,i,L,x,y,z),
	McDOrder(McDOrder),
	Na(2*McDOrder + 2),
	L_2(0.5*L),
	C(PhysicsConstants::mu0*i*N/(2*L*(R2-R1)))
	{	
		//~ std::cout << "Constructing McD Class..." << std::endl;
		int n = 2 + 2*McDOrder; // To calc the (n_McD)'th order in the McD expansion, (2+2*n_McD) orders are needed for a_n
		// for n = 2 and higher there is a system to an[]. Each terms will be of the form k*pow(z,lambda)*pow(a,mu)*pow(b,nu)
		if(n > 1){
			
		NT = 0; // the total number of terms contributed by the 2nd to the n'th order
		for(int i=2; i<n; i++){
			NT += pow(3,i-1); // the number of terms contributed by the n'th order
		}
		
		// Each term gives rise to another 3 terms for the next order. Calculate the remaining parameters:
		int NTAcc_i = 0; 				// placeholder for the accumulated number of terms to the (n_i)'th order (for n_i > 1)
		for (int n_i = 2; n_i < n; n_i++){ 	// for the 2nd and the remaining orders
			const int NT_im1 = pow(3,n_i-2);		// the number of terms associated with the (n_i-1)'th order
			const int NT_i = pow(3,n_i-1); 		// the number of terms associated with the (n_i)'th order			
			//~ std::cout << "n_i = " << n_i << std::endl;
			
			if (n_i > 2){					// The parameters for the 2nd order have already been calculated, so skip this step if n_i = 2
			for (int i = 0; i < NT_im1; i++){ 		// Loop over all terms in the previous order (parent terms) and calc the terms derived from them
				const int P_i = NTAcc_i - NT_im1 + i; 	// Parent index. The index of the term that gives rise to the next 3 terms	
				
				k_arr.push_back( 		k_arr[P_i]*lambda_arr[P_i] 	);	
				lambda_arr.push_back( 	lambda_arr[P_i] - 1			); 	
				mu_arr.push_back( 		mu_arr[P_i]					); 	
				nu_arr.push_back( 		nu_arr[P_i]					);
				Ti_arr.push_back( 		n_i							);
				
				k_arr.push_back( 		-k_arr[P_i]*nu_arr[P_i]		);	
				lambda_arr.push_back( 	lambda_arr[P_i] + 1			); 	
				mu_arr.push_back( 		mu_arr[P_i] + 1				); 	
				nu_arr.push_back( 		nu_arr[P_i] + 1				);
				Ti_arr.push_back( 		n_i							);
				
				k_arr.push_back( 		-k_arr[P_i]*mu_arr[P_i]		);	
				lambda_arr.push_back( 	lambda_arr[P_i] + 1			); 	
				mu_arr.push_back(		mu_arr[P_i] + 2				); 
				nu_arr.push_back( 		nu_arr[P_i]					);
				Ti_arr.push_back( 		n_i							);				
			}
			}
			// The parameters of the terms of the (n_i)'th order are now calculated

			NTAcc_i += NT_i;
		}
		//~ std::cout << "Calculation of terms done" << std::endl;
		// now that all the terms have been calculated it is time to group the terms together
		
		//~ std::cout << "k size = " << k_arr.size() << std::endl;
		//~ std::cout << "lambda size = " << lambda_arr.size() << std::endl;
		//~ std::cout << "mu size = " << mu_arr.size() << std::endl;
		//~ std::cout << "nu size = " << nu_arr.size() << std::endl;
		assert(k_arr.size() == lambda_arr.size() && lambda_arr.size() == mu_arr.size() && mu_arr.size() == nu_arr.size());
		
		unsigned int i = 0;
		while (i < k_arr.size()){
			//~ std::cout << "i = " << i << "/" << k_arr.size() << std::endl;
			const int k_i = k_arr[i];
			const int lambda_i = lambda_arr[i];
			const int mu_i = mu_arr[i];
			const int nu_i = nu_arr[i];
			const int Ti_i = Ti_arr[i];
			
			if(k_i == 0){
				k_arr.erase(k_arr.begin()+i); 	// delete this term since the coefficient is 0
				lambda_arr.erase(lambda_arr.begin()+i);
				mu_arr.erase(mu_arr.begin()+i);
				nu_arr.erase(nu_arr.begin()+i);
				Ti_arr.erase(Ti_arr.begin()+i);			
			} else {
				unsigned int j = i + 1; // an index to scan throught the remaining elements
				while (j < k_arr.size()){
					if(lambda_arr[j] == lambda_i && mu_arr[j] == mu_i && nu_arr[j] == nu_i && Ti_arr[j] == Ti_i){
						k_arr[i] += k_arr[j]; 	// the two terms are the same and belong to the same order. Add their front factors together...
						k_arr.erase(k_arr.begin()+j); 	// ... and delete the matched element
						lambda_arr.erase(lambda_arr.begin()+j);
						mu_arr.erase(mu_arr.begin()+j);
						nu_arr.erase(nu_arr.begin()+j);
						Ti_arr.erase(Ti_arr.begin()+j);
					}else{
						j += 1;
					}				
				}			
				i += 1;
			}
		}
		//~ std::cout << "Grouping of terms done" << std::endl;
		
		//~ std::cout << "k size = " << k_arr.size() << std::endl;
		//~ std::cout << "lambda size = " << lambda_arr.size() << std::endl;
		//~ std::cout << "mu size = " << mu_arr.size() << std::endl;
		//~ std::cout << "nu size = " << nu_arr.size() << std::endl;
		assert(k_arr.size() == lambda_arr.size() && lambda_arr.size() == mu_arr.size() && mu_arr.size() == nu_arr.size());

		NT = k_arr.size();			
		}
		//~ std::cout << "Constructing McD Class done" << std::endl;
	} // end of constructor
	
	~McD_Tube(){
	} // end of destructor
	
	void getA0(const double z, const double R, double a0_array[]) const;
	void getB(const double cylP[3], double BCylVec[3]) const ;
};

class TAVP : public Loop{
	private:
	const double lambda;
	
	public:
	TAVP(const double lambda, const double R, const double I, const double x, const double y, const double z) : Loop(R,I,x,y,z), lambda(lambda){
	
	} // end of constructor
	
	~TAVP(){
	} // end of destructor
	
	void getB(double sphP[3], double BSphVec[3]) const ;
};

class Helix : public Tube, public BiotSavart{
	private:
	const int N_z; // the number of wires making up the Tube
	const int N_rho; // the number of layers making up the Tube
	const int N_BS; // the number of straigt line segments per turn
	const double delta_rho; // the spacing between layers
	const double delta_z; // the spapring between wires in a layer	
	
	public:
	Helix(const int N_z, const int N_rho, const int N_BS, const double R1, const double R2, 
			   const int N, const double i, const double L, const double x, const double y, const double z) : 
	Tube(R1,R2,N,i,L,x,y,z),
	N_z(N_z), 
	N_rho(N_rho), 
	N_BS(N_BS),
	delta_rho(Tube::getThickness()/N_rho),
	delta_z(L/N_z)
	{
	
	} // end of constructor
	
	~Helix(){
	} // end of destructor
	
	void getB(const double carP[3], double BCarVec[3]) const ;
};

class NWire_Tube : public Tube{
	private:
	const int N_z; // the number of wires making up the Tube
	const int N_rho; // the number of layers making up the Tube
	const int NGP_z; // the number of wires used to represent the Tube in the Gaussian Quadrature
	const int NGP_rho; // the number of layers used to represent the Tube in the Gaussian Quadrature
	std::vector<std::vector<SimpleAnalyticModel>> loops;
	
	
	public:
	NWire_Tube(const int N_z, const int N_rho, const int NGP_z, const int NGP_rho, const double R1, const double R2, 
			   const int N, const double i, const double L, const double x, const double y, const double z) : Tube(R1,R2,N,i,L,x,y,z),
			   N_z(N_z), N_rho(N_rho), NGP_z(NGP_z), NGP_rho(NGP_rho){
		
		// Dimension of Tube			
		const double width_rho = this->getThickness();
		
	    // Number of wires
		const int N_wires = N_rho*N_z;
		
		const double I = i*N/N_wires; // total current in tube distributed across the N wires
		
	    // Length
		const double width_z = this->getL(); 	
	
	    // Spacing			
		const double delta_rho = width_rho/N_rho; 	
		const double delta_z = width_z/N_z; 
				
		for(int n_rho = 0; n_rho < N_rho; n_rho++){
			loops.push_back(std::vector<SimpleAnalyticModel>());
			for(int n_z = 0; n_z < N_z; n_z++){
				loops[n_rho].push_back(SimpleAnalyticModel(R1 + delta_rho*0.5 + n_rho*delta_rho, I, x, y, z - 0.5*width_z + delta_z*0.5 + n_z*delta_z));
			}
		}
	} // end of constructor
	
	~NWire_Tube(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ;
};

class GaussianQuadratureLoops_Tube : public Tube, public GQ_Support{
	private:
	const int N_z; // the number of wires making up the Tube
	const int N_rho; // the number of layers making up the Tube
	const int NGP_z; // the number of wires used to represent the Tube in the Gaussian Quadrature
	const int NGP_rho; // the number of layers used to represent the Tube in the Gaussian Quadrature
	std::vector<double> GPZValues;
	std::vector<double> GPZWeights;
	std::vector<double> GPRhoValues;
	std::vector<double> GPRhoWeights;
	std::vector<std::vector<SimpleAnalyticModel>> loops;
	
	public:
	GaussianQuadratureLoops_Tube(const int N_z, const int N_rho, const int NGP_z, const int NGP_rho, const double R1, const double R2, 
								 const int N, const double i, const double L, const double x, const double y, const double z) : Tube(R1,R2,N,i,L,x,y,z),
								 N_z(N_z), N_rho(N_rho), NGP_z(NGP_z), NGP_rho(NGP_rho){
		assert(N_z % 2 == 0); // N_z has to be even, as the GQ alogorithm assumes the same number of wires in each half of the magnet
		assert(N_rho % 2 == 0); // N_rho has to be even, as the GQ alogorithm assumes the same number of wires in each half of the magnet
		
		const double width_z = this->getL(); 	
		const double width_rho = this->getThickness();
		
		const int NWiresRho = N_rho/2; // half the number of wires in each dimension
		const int NWiresZ = N_z/2;		
		
		getGaussianQuadratureParams(NGP_rho,GPRhoValues,GPRhoWeights,width_rho*0.5,NWiresRho);
		getGaussianQuadratureParams(NGP_z,GPZValues,GPZWeights,width_z*0.5,NWiresZ);
		
		const double I_temp = this->getI()/(N_rho*N_z);	// divide the current of the tube between the loops use for GQ  (N_rho*N_z loops)
		const double centre_rho = R1 + width_rho*0.5;
		
		for(int nGP_rho = 0; nGP_rho < NGP_rho; nGP_rho++){
			loops.push_back(std::vector<SimpleAnalyticModel>());
			for(int nGP_z = 0; nGP_z < NGP_z; nGP_z++){
				loops[nGP_rho].push_back(SimpleAnalyticModel(centre_rho + GPRhoValues[nGP_rho],I_temp,x,y,z+GPZValues[nGP_z]));
			}
		}
		
	} // end of constructor
	
	~GaussianQuadratureLoops_Tube(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ;
};

class GaussianQuadratureShells_Tube : public Tube, public GQ_Support{
	private:
	const int N_rho; // the number of layers making up the Tube
	const int NGP_rho; // the number of shells used to represent the Tube in the Gaussian Quadrature
	std::vector<double> GPRhoValues;
	std::vector<double> GPRhoWeights;
	
	std::vector<Conway1D> shells;
	
	public:
	GaussianQuadratureShells_Tube(const int N_rho, const int NGP_rho, const double R1, const double R2, const int N, const double i, const double L,
								  const double x, const double y, const double z) : Tube(R1,R2,N,i,L,x,y,z), N_rho(N_rho), NGP_rho(NGP_rho) {
		
		const double width_rho = this->getThickness();
		const int NWiresRho = N_rho/2; // half the number of wires in each dimension
		getGaussianQuadratureParams(NGP_rho,GPRhoValues,GPRhoWeights,width_rho*0.5,NWiresRho);
		
		const double I = this->getI()/N_rho; // devide the current of the tube on N_rho shells
		//~ const double innerRadius = this->getR1();
		//~ const double width_z = this->getL(); 	
		const double centre_rho = R1 + width_rho*0.5;
		
		for(int nGP_rho = 0; nGP_rho < NGP_rho; nGP_rho++){
			shells.push_back(Conway1D(centre_rho + GPRhoValues[nGP_rho],1,I,L,x,y,z));
		}
	
	} // end of constructor
	
	~GaussianQuadratureShells_Tube(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ;
};


#endif
