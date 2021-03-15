// A class to hold all the coefficients that can be precalculated for the McDonald expansion for
// a finite solenoid (a solenoid with finite length and  thickness). It also holds a method 
// to calculate the magnetic field

#ifndef MCD_TUBE_SUPPORT_H
#define MCD_TUBE_SUPPORT_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <cstdlib>

#include "Tube.h"
//~ #include "printFuncs.h"

class McD_Tube_Support {
private:
	const int McD_order;
	std::vector<int> k_arr = {3, -1, -1};
	std::vector<int> lambda_arr = {1, 3, 3};
	std::vector<int> mu_arr = {1, 2, 3};
	std::vector<int> nu_arr = {1, 2, 1};
	std::vector<int> Ti_arr = {2, 2, 2};
	int NT; // the total number of terms in a_n


public:
	McD_Tube_Support(const int McDOrder):McD_order(McDOrder)
	{	
		int n = 2 + 2*McD_order; // To calc the (n_McD)'th order in the McD expansion, (2+2*n_McD) orders are needed for a_n
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
		
		// now that all the terms have been calculated it is time to group the terms together
		
		//~ std::cout << "k size = " << k_arr.size() << std::endl;
		//~ std::cout << "lambda size = " << lambda_arr.size() << std::endl;
		//~ std::cout << "mu size = " << mu_arr.size() << std::endl;
		//~ std::cout << "nu size = " << nu_arr.size() << std::endl;
		assert(k_arr.size() == lambda_arr.size() && lambda_arr.size() == mu_arr.size() && mu_arr.size() == nu_arr.size());
		
		unsigned int i = 0;
		while (i < k_arr.size()){
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
		
		//~ std::cout << "k size = " << k_arr.size() << std::endl;
		//~ std::cout << "lambda size = " << lambda_arr.size() << std::endl;
		//~ std::cout << "mu size = " << mu_arr.size() << std::endl;
		//~ std::cout << "nu size = " << nu_arr.size() << std::endl;
		assert(k_arr.size() == lambda_arr.size() && lambda_arr.size() == mu_arr.size() && mu_arr.size() == nu_arr.size());

		NT = k_arr.size();			
		}
	} // end of constructor
	
	~McD_Tube_Support(){
	} // end of destructor
	
	void getA0(const double z, const double R, const int n_McD, double a0_array[]);
	void getB(Tube tube, const int nmax, const CylindricalVector& cylP, CylindricalVector& BCylVec);
	
};



#endif
