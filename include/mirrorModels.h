#ifndef MIRRORMODELS_H
#define MIRRORMODELS_H

#include <iostream>
#include "Loop.h"
#include "Shell.h"
#include "BFieldModels.h"
#include "calcWrappers.h"
#include "writeToFile.h"
#include "mirrorMcD_v1.h"
#include "printFuncs.h"
#include "ALPHAPhysicalConstants.hpp"
#include "coorTransf.h"
#include <cassert>
#include "easy_eos.h"

//NWire models
void exactMirrorModel(std::string path, double carP1[3], double carP2[3]);
void exactMirrorModel_vecClass();

// Gaussian Quadrature models
void exactMirrorModel_GaussQuad_Loop_Line(std::string path, double carP1[3], double carP2[3]);
void exactMirrorModel_GaussQuad_Shell(std::string path, double carP1[3], double carP2[3]);
void getGaussianQuadratureParams(const int N, double points[], double weights[], const double dimLength, const double NWires);

// A detailed Biot-Savart model of a solenoid. Wire wound in helix shape
void exactMirrorModel_BiotSavart(std::string path, double carP1[3], double carP2[3]);
void exactMirrorModel_BiotSavart(double carP[3], double BCarVec[3]);
void testExactMirror();



#endif
