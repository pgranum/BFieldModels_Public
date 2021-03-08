#ifndef BFIELDMODELS_H
#define BFIELDMODELS_H

#include <iostream>
#include <cmath>
#include <math.h>
#include "mathFuncs.h"
#include "Loop.h"
#include "Shell.h"
#include "Tube.h"
#include "Generic3Vector.hpp"
#include <tr1/cmath>
#include <tgmath.h>
#include <chrono>
#include "coorTransf.h"
#include "writeToFile.h"
#include "ALPHAPhysicalConstants.hpp"
#include "printFuncs.h"
#include "mcDTubeSupportClass.h"

// functions should later be added to the magnet parent class.
// keep here for now to avoid merge conflicts during development

// Biot-Savart law
void biotSavart(double s0[3], double s1[3], double r[3], const double I, double B_vec[3]);
void biotSavart_vecClass(const CartesianVector& s0, const CartesianVector& s1, const CartesianVector& r, double I, CartesianVector& B_vec);

// Misc Models
void currentLoopFieldOnAxisAnalytic(double R, double I, double z, double B_vec[3]);
void currentLoopFieldOnAxisAnalytic(const double R, const double I, const double z, CartesianVector& B_vec);
void lineSegmentsOct(double R, double L, double z, double dPhi, const int N, double s_arr[][3]);
void mirrorBFieldLine_BS(Loop loop, double r0[3], double r1[3]);

// Current Loop
void loopExactSAM(Loop loop, const double cylP[3], double BCylVec[3]);
void loopExactSAM_vecClass(Loop loop, CylindricalVector& cylP, CylindricalVector& BCylVec);
const CartesianVector SAM(Loop loop, CartesianVector& carP);

void loopApproxMcDonald_vecClass(double z, const double I, const double R, const int nmax, const CylindricalVector& cylP, CylindricalVector& BCylVec);
void mcDonald(Loop loop, const int nmax, const double cylP[3], double BCylVec[3]);
void mcDonaldLoopSupFunc(const int n, const double z, const double R, double an[]);

void mcDonald(Shell shell, const int nmax, const double cylP[3], double BCylVec[3]);
void mcDonaldShellSupFunc(const int n, const double z, const double R, double an[]);

void mcDonald(Tube tube, const int nmax, const double cylP[3], double BCylVec[3]);
void mcDonald(Tube tube, const int nmax, const double cylP[3], double BCylVec[3], McD_Tube_Support& McDSupport);
void mcDonaldTubeSupFunc(const int n, const double z, const double R, double an[]);

void loopBiotSavart(Loop loop, const int NSegements, double carP[3], double BCarVec[3]);
void loopExactJack510(Loop loop, double cylP[3], double BCylVec[3]);
void loopApproxFrancis(Loop loop, const double lam, const double sphP[3], double BSphVec[3]);
void loopApproxJackson(Loop loop, double cylP[3], double BCylVec[3]);
void Conway1D(Shell shell, const double cylP[3], double BCylVec[3]);

#endif
