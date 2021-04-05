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
	
	void getGaussianQuadratureParams(const int N, double points[], double weights[], const double dimLength, const double NWires) const ;
};

// LOOP

class SimpleAnalyticModel : public Loop{
	private:
	
	public:
	SimpleAnalyticModel(const double R, const double I, const double x, const double y, const double z) : Loop(R,I,x,y,z){
	
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
	
	public:
	McD_Loop(const int McDOrder, const double R, const double I, const double x, const double y, const double z) : Loop(R,I,x,y,z), McDOrder(McDOrder){
	
	} // end of constructor
	
	~McD_Loop(){
	} // end of destructor
	
	void getB(const double carP[3], double BCarVec[3]) const ;
	void mcDonaldLoopSupFunc(const int n, const double z, const double R, double an[]) const ;
};

// SHELL

class Conway1D : public Shell{
	private:
	
	public:
	Conway1D(const double R, const int N, const double i, const double L, const double x, const double y, const double z) : Shell(R,N,i,L,x,y,z){
	
	} // end of constructor
	
	~Conway1D(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ;
};

class McD_Shell : public Shell {
	private:
	const int McDOrder;
	const double L_2 = 0.5*this->getL();
	
	public:
	McD_Shell(const int McDOrder, const double R, const int N, const double i, const double L, 
			  const double x, const double y, const double z): Shell(R,N,i,L,x,y,z), McDOrder(McDOrder){
	} // end of constructor
	
	~McD_Shell(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ; 
	void mcDonaldShellSupFunc(const int n, const double z, double an[]) const ;
};

class NWire_Shell : public Shell{
	private:
	const int N_z; // the number of wires in the shell
	
	public:
	NWire_Shell(const int N_z, const double R, const int N, const double i, const double L, 
				const double x, const double y, const double z) : Shell(R,N,i,L,x,y,z), N_z(N_z){
	
	} // end of constructor
	
	~NWire_Shell(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ;
};

class GaussianQuadratureLoops_Shell : public Shell, public GQ_Support{
	private:
	const int N_z; // the number of wires making up the Shell
	const int NG_z; // the number of wires used to represent the shell in the Gaussian Quadrature
	
	public:
	GaussianQuadratureLoops_Shell(const int N_z, const int NG_z, const double R, const int N, const double i, const double L,
								  const double x, const double y, const double z) : Shell(R,N,i,L,x,y,z), N_z(N_z), NG_z(NG_z){
	
	} // end of constructor
	
	~GaussianQuadratureLoops_Shell(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ;
};

// TUBE

class McD_Tube : public Tube{
	private:
	const int McDOrder;
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
			 const int N, const double i, const double L, const double x, const double y, const double z) : Tube(R1,R2,N,i,L,x,y,z),
			 McDOrder(McDOrder){	
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
	
	void getB(const double sphP[3], double BSphVec[3]) const ;
};

class Helix : public Tube, public BiotSavart{
	private:
	const int N_z; // the number of wires making up the Tube
	const int N_rho; // the number of layers making up the Tube
	const int N_BS; // the number of straigt line segments per turn
	
	public:
	Helix(const int N_z, const int N_rho, const int N_BS, const double R1, const double R2, 
			   const int N, const double i, const double L, const double x, const double y, const double z) : Tube(R1,R2,N,i,L,x,y,z),
			   N_z(N_z), N_rho(N_rho), N_BS(N_BS){
	
	} // end of constructor
	
	~Helix(){
	} // end of destructor
	
	void getB(const double carP[3], double BCarVec[3]) const ;
};

class NWire_Tube : public Tube{
	private:
	const int N_z; // the number of wires making up the Tube
	const int N_rho; // the number of layers making up the Tube
	const int NG_z; // the number of wires used to represent the Tube in the Gaussian Quadrature
	const int NG_rho; // the number of layers used to represent the Tube in the Gaussian Quadrature
	
	public:
	NWire_Tube(const int N_z, const int N_rho, const int NG_z, const int NG_rho, const double R1, const double R2, 
			   const int N, const double i, const double L, const double x, const double y, const double z) : Tube(R1,R2,N,i,L,x,y,z),
			   N_z(N_z), N_rho(N_rho), NG_z(NG_z), NG_rho(NG_rho){
	
	} // end of constructor
	
	~NWire_Tube(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ;
};

class GaussianQuadratureLoops_Tube : public Tube, public GQ_Support{
	private:
	const int N_z; // the number of wires making up the Tube
	const int N_rho; // the number of layers making up the Tube
	const int NG_z; // the number of wires used to represent the Tube in the Gaussian Quadrature
	const int NG_rho; // the number of layers used to represent the Tube in the Gaussian Quadrature
	
	public:
	GaussianQuadratureLoops_Tube(const int N_z, const int N_rho, const int NG_z, const int NG_rho, const double R1, const double R2, 
								 const int N, const double i, const double L, const double x, const double y, const double z) : Tube(R1,R2,N,i,L,x,y,z),
								 N_z(N_z), N_rho(N_rho), NG_z(NG_z), NG_rho(NG_rho){
	
	} // end of constructor
	
	~GaussianQuadratureLoops_Tube(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ;
};

class GaussianQuadratureShells_Tube : public Tube, public GQ_Support{
	private:
	const int N_rho; // the number of layers making up the Tube
	const int NG_rho; // the number of shells used to represent the Tube in the Gaussian Quadrature
	
	public:
	GaussianQuadratureShells_Tube(const int N_rho, const int NG_rho, const double R1, const double R2, const int N, const double i, const double L,
								  const double x, const double y, const double z) : Tube(R1,R2,N,i,L,x,y,z), N_rho(N_rho), NG_rho(NG_rho){
	
	} // end of constructor
	
	~GaussianQuadratureShells_Tube(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]) const ;
};



#endif
