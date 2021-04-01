#ifndef BFIELDMODELS_H
#define BFIELDMODELS_H

#include <iostream>
#include <cmath>
#include <math.h>
#include <tr1/cmath>
#include <tgmath.h>
#include <chrono>

#include "Loop.h"
#include "Shell.h"
#include "Tube.h"
#include "PhysicsConstants.hpp"
#include "mcDTubeSupportClass.h"
#include "mathFuncs.h"

// Basic

void biotSavart(const double s0[3], const double s1[3], const double r[3], const double I, double B_vec[3]);

// Current Loop
void loopExactSAM(const Loop& loop, const double cylP[3], double BCylVec[3]);

void loopBiotSavart(const Loop& loop, const int NSegments, const double carP[3], double BCarVec[3]);

void mcDonald(const Loop& loop, const int nmax, const double cylP[3], double BCylVec[3]);
void mcDonaldLoopSupFunc(const int n, const double z, const double R, double an[]);

// Shell

void Conway1D(const Shell& shell, const double cylP[3], double BCylVec[3]);

//~ void mcDonald(const Shell& shell, const int nmax, const double cylP[3], double BCylVec[3]);
//~ void mcDonaldShellSupFunc(const int n, const double z, const double R, double an[]);

void NWire(const Shell& shell, const int N_z, const double cylP[3], double BCylVec[3]);

void GaussianQuadratureLoop(const Shell& shell, const int N_z, const int NG_z, const double cylP[3], double BCylVec[3]);

// Tube

void Helix(const Tube& tube, const int N_rho, const int N_z, const int N_BS, const double carP[3], double BCarVec[3]);
void mcDonald(const Tube& tube, const double cylP[3], double BCylVec[3], const McD_Tube_Support& McDSupport);

void loopApproxFrancis(const Loop& loop, const double lam, const double sphP[3], double BSphVec[3]);

void NWire(const Tube& tube, const int N_rho, const int N_z, const double cylP[3], double BCylVec[3]);

void GaussianQuadratureLoop(const Tube& tube, const int N_rho, const int N_z, const int NG_rho, const int NG_z, const double cylP[3], double BCylVec[3]);
void GaussianQuadratureShell(const Tube& tube, const int N_rho, const int NG_rho, const double cylP[3], double BCylVec[3]);

// Utils

void getGaussianQuadratureParams(const int N, double points[], double weights[], const double dimLength, const double NWires);

// Classes

class McDShell : public Shell {
	private:
	const int McDOrder;
	
	public:
	McDShell(const int McDOrder, const double R, const int N, const double i, const double L, const double x, const double y, const double z):McDOrder(McDOrder), Shell(R,N,i,L,x,y,z){
	} // end of constructor
	
	~McDShell(){
	} // end of destructor
	
	void getB(const double cylP[3], double BCylVec[3]); // what does adding const here do?
	void mcDonaldShellSupFunc(const int n, const double z, const double R, double an[]);
};


#endif
