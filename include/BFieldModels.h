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

void biotSavart(double s0[3], double s1[3], double r[3], const double I, double B_vec[3]);

// Current Loop
void loopExactSAM(Loop loop, const double cylP[3], double BCylVec[3]);

void loopBiotSavart(Loop loop, const int NSegments, double carP[3], double BCarVec[3]);

void mcDonald(Loop loop, const int nmax, const double cylP[3], double BCylVec[3]);
void mcDonaldLoopSupFunc(const int n, const double z, const double R, double an[]);

// Shell

void Conway1D(Shell shell, const double cylP[3], double BCylVec[3]);

void mcDonald(Shell shell, const int nmax, const double cylP[3], double BCylVec[3]);
void mcDonaldShellSupFunc(const int n, const double z, const double R, double an[]);

// Tube

void Helix(Tube tube, const int N_rho, const int N_z, const int N_BS, double carP[3], double BCarVec[3]);

void mcDonald(Tube tube, const int nmax, const double cylP[3], double BCylVec[3], McD_Tube_Support& McDSupport);
void mcDonaldTubeSupFunc(const int n, const double z, const double R, double an[]);

void loopApproxFrancis(Loop loop, const double lam, const double sphP[3], double BSphVec[3]);

// Shell and Tube

void NWire(Tube tube, const int N_rho, const int N_z, const double cylP[3], double BCylVec[3]);

void GaussianQuadratureLoop(Tube tube, const int N_rho, const int N_z, const int NG_rho, const int NG_z, const double cylP[3], double BCylVec[3]);
void GaussianQuadratureShell(Tube tube, const int N_rho, const int N_z,  const int NG_rho, const double cylP[3], double BCylVec[3]);
void getGaussianQuadratureParams(const int N, double points[], double weights[], const double dimLength, const double NWires);

#endif
