#include <iostream>
#include <cmath>
#include <limits>

#include "BFieldModels.h"
#include "mathFuncs.h"
//~ #include "Loop.h"
//~ #include "Shell.h"
//~ #include "Tube.h"
#include "utils.h"

void gradDesc(double a_n[], int N, double a_final[]);
double f(double xVec[2]);
void dfdx(const double xVec[], const int N, double dfdx[]);

void gradDesc(double a_n[], int N, double a_final[]){
	
	double gamma = 1;
	double err = 1;
	double eps = std::numeric_limits<double>::epsilon();
	
	// predefine a_n-1 and a_n+1
	double a_np1[N];
	double a_nm1[N];
	for(int i=0; i<N; i++){
		a_nm1[i] = a_n[i] - 100*eps;
		a_np1[i] = 0;
	}
	
	// placeholders for gradients
	double grad_nm1[N];
	double grad_n[N];
	double grad_np1[N];
	for (int i = 0; i < N; i++)
	{
		grad_nm1[i] = 0;
		grad_n[i] = 0;
		grad_np1[i] = 0;
	}
	// precalc gradients for a_n-1 and a_n
	dfdx(a_nm1,N,grad_nm1);
	dfdx(a_n,N,grad_n);
	return;
	int n = 0; // iteration counter
	while(fabs(err) > 1e-12 && n < 10){
		//~ printArr(a_min,N,"a_min");

		
		err = 0;
		gamma = 0; // reset gamma
		double denom = 0;
		for(int i=0; i<N; i++){
			gamma += (a_n[i]-a_nm1[i])*(grad_n[i]-grad_nm1[i]);
			denom += pow(grad_n[i]-grad_nm1[i],2);
		}
		gamma /= denom;
		//~ gamma = 0.01;
		std::cout << "gamma = " << gamma << std::endl;
		std::cout << "deltaA = " << gamma*grad_n[0] << std::endl;
		std::cout << "a_n = " << a_n[0] << std::endl;
		std::cout << "a_n+1 = " << a_n[0] - gamma*grad_n[0] << std::endl;
		
		
		// calc a_n+1
		for(int i=0; i<N; i++){
			a_np1[i] = a_n[i] - gamma*grad_n[i];
		}
		
		printArr(a_np1,N,"a_n+1");
		
		// find grad at a_n+1
		dfdx(a_np1,N,grad_np1);
		
		err = f(a_n) - f(a_np1); // ERROR SHOULD BE THE DIFFERENCE BETWEEN THE FUNCTION VALUES
		std::cout << "err = " << err << std::endl;
		
		// update a_n to be a_n+1 and a_nm1 to be a_n
		for(int i=0; i<N; i++){
			a_nm1[i] = a_n[i];
			grad_nm1[i] = grad_n[i];
			a_n[i] = a_np1[i];
			grad_n[i] = grad_np1[i];
			
		}
		
		
		
		n += 1;		
		if(fabs(err) > 1e-12 && n < 10){std::cout << "we continue" << std::endl;}
		if(!(fabs(err) > 1e-12 && n < 10)){std::cout << "we stop" << std::endl;}
		std::cout << "" << std::endl;
		
	}
	
	std::cout << "n = " << n << std::endl;
	
	for(int i=0; i<N; i++){
		a_final[i] = a_n[i];
	}
	

}

double f(double xVec[2]){
	double result = pow(xVec[0]-1,2) + pow(xVec[1]-10,2);
	//~ double result = pow(xVec[0]-1,2);
	return result;
}

void dfdx(const double xVec[], const int N, double dfdx[]){
	double eps = 100000*std::numeric_limits<double>::epsilon();
	eps = 0.1;
	double x1[N];
	double x2[N];
	
	
	
	const double R1 = 0.04125;			// inner radius
	const double R2 = 0.0463701;		// outer radius
	const double L = 0.0346811; 		// Length
	const double current = 600; 		// current per loop/current in the wire
	const int N_z = 30; 				// number of wire turns/loops per layer/shell
	const int N_rho = 4; 				// number of layers/shell
	const int N_wires = N_z*N_rho;		// total number of wires
	//~ const double I = N_z*N_rho*current;	// "total" amount of current
	const double x = 0.;				// x coordinate of magnet centre
	const double y = 0.;				// y coordinate of magnet centre
	const double z = 0.;				// z coordinate of magnet centre			
	const int McDOrder = 5;				// number of terms to use in the McDonald model
	const int N_BS = 10000;				// number of segments to be used in the Biot-Savart model

	Helix helix = Helix(N_z,N_rho,N_BS,R1,R2,N_wires,current,L,x,y,z); // IN THEORY THIS SHOULD ONLY BE INSTANTIATED ONCE. MAKE IS AS A CLASS
	double BCarBS[3] = {0.,0.,0.}; 		// placeholder for the field in cartesian coordinates	
		
	// loop over all parameters
	for(int i=0; i<N; i++){
		// copy the default function parameters
		for(int j=0; j<N; j++){
			x1[j] = xVec[j];
			x2[j] = xVec[j];		
		}
		
		x1[i] -= eps;
		x2[i] += eps;
		printArr(x1,N,"x1");
		printArr(x2,N,"x2");
		
		//~ double f1 = f(x1);
		//~ double f2 = f(x2);
		
		
		
		McD_Tube mcD_Tube1 = McD_Tube(McDOrder,R1+x1[0],R2+x1[0],N_wires,current,L+x1[1],x,y,z);
		McD_Tube mcD_Tube2 = McD_Tube(McDOrder,R1+x2[0],R2+x2[0],N_wires,current,L+x2[1],x,y,z);
		double BCarMcD1[3] = {0.,0.,0.};
		double BCarMcD2[3] = {0.,0.,0.};
		double BCylMcD1[3] = {0.,0.,0.};
		double BCylMcD2[3] = {0.,0.,0.};
		
		double f1 = 0.0;
		double f2 = 0.0;
		
		const int Nk = 100;
		for(int k=0; k<Nk; k++){
			double carP[3] = {0.,0.,-5*R1 + (double)k/(Nk-1.0)*10*R1}; 	// point to calculate field at in cartesian coordinates
			double cylP[3];										// point to calculate field at in cylindrical coordinates
			carPToCylP(carP,cylP);								// converting the point in cartesian coor to cylindrical coordinates	
			printVec(cylP,"cylP_k");
			
			
			mcD_Tube1.getB(cylP,BCylMcD1);
			mcD_Tube2.getB(cylP,BCylMcD2);
			printVec(BCylMcD1,"BMcD1");
			printVec(BCylMcD2,"BMcD2");
			cylVecToCarVec(BCylMcD1,cylP,BCarMcD1);
			cylVecToCarVec(BCylMcD2,cylP,BCarMcD2);
			helix.getB(carP,BCarBS);		
			f1 += vecNorm(BCarBS) - vecNorm(BCarMcD1);
			f2 += vecNorm(BCarBS) - vecNorm(BCarMcD2);
		}
		
		std::cout << "f1 " << f1 << std::endl;
		std::cout << "f2 " << f2 << std::endl;		
		
		dfdx[i] = (f2-f1)/2.0*eps;
		
		std::cout << "dfdx_" << i << " = " << dfdx[i] << std::endl;
		
	}
}


int main(){
	std::cout.precision(15);	// sets the number of significant digits of cout
	std::cout << "Hello world\n";
	
	const int N = 2;
	double a_0[N] = {0.02,0.02};
	double a_min[N] = {0.,0.};
	//~ int N = 1;
	//~ double a_0[N] = {9};
	//~ double a_min[N] = {0};
	
	
	gradDesc(a_0, N, a_min);
	
	printArr(a_min,N,"a_min");
	
	return 0;
}
