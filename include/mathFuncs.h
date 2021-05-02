#ifndef MATH_FUNCS_H
#define MATH_FUNCS_H

#include <iostream>
#include <cmath>
#include <tr1/cmath> // adds the special elliptic integrals

#include "PhysicsConstants.hpp"

// vector manipulation

//~ void vecAdd(const double a[3], const double b[3], double c[3]);
//~ void vecAddOvrwrt(double a[3], const double b[3]);
void vecSub(const double a[3], const double b[3], double c[3]);
//~ void vecMult(const double a[3], const double b[3], double c[3]);
double vecDotP(const double a[3], const double b[3]);
//~ void vecCrsP(const double  a[3], const double  b[3], double  c[3]);
double vecNorm(const double a[3]);
//~ void vecMultScal(const double a[3], const double k, double c[3]);
//~ void vecMultScalOvrwrt(double a[3], const double k);

// Conversion of points

void carPToCylP(const double carP[3], double cylP[3]);
void carPToSphP(const double carP[3], double sphP[3]);
void cylPToCarP(const double cylP[3], double carP[3]);
void sphPToCarP(const double sphP[3], double carP[3]);

// Conversion of vectors

void carVecToCylVec(const double carVec[3], const double carP[3], double cylVec[3]);
void carVecToSphVec(const double carVec[3], const double carP[3], double sphVec[3]);
void cylVecToCarVec(const double cylVec[3], const double cylP[3], double carVec[3]);
void sphVecToCarVec(const double sphVec[3], const double sphP[3], double carVec[3]);

//~ // interpolation

//~ void calcPointsLine3D(double r1[3], double r2[3], int N, double r_arr[][3]);
//~ void calcPointsLine3D(CartesianVector r1, CartesianVector r2, int N, std::vector<CartesianVector> *points);
//~ void calcPointsLine1D(double r1, double r2, int N, double r_arr[]);
//~ void calcPointsRhoZCylinder(  double phi, double rho_min, double rho_max, double z_min, double z_max, int NRho, int NZ, double arr[][3]);
//~ void calcPointsPhiZCylinder(  double rho, double phi_min, double phi_max, double z_min, double z_max, int NPhi, int NZ, double arr[][3]);
//~ void calcPointsRhoPhiCylinder(double z,   double rho_min, double rho_max, double phi_min, double phi_max, int NRho, int NPhi, double arr[][3]);

//~ // linspace, trapz and diff
//~ void linspace(double min, double  max, int N, double result[]);
//~ double trapz(int N, double xarr[], double yarr[]);
//~ double trapz(func f, double a, double b, int N);
//~ void cumtrapz(int N, double xarr[], double yarr[], double areaArr[]);
//~ void gradient(int N, double xarr[], double yarr[], double gradArr[]);

//~ // adaptive integration methods
//~ double right_riemann(func f, double a, double b, double tol);
//~ double left_riemann(func f, double a, double b, double tol);
//~ double recur_left_riemann(func f, double a, double fa, double b, double tol);
//~ double trapz(func f, double a, double b, double tol);
//~ double recur_trapz(func f, double a, double fa, double b, double fb, double tol);
//~ double simpson(func f, double a , double b, double tol);
//~ double recur_simpson(func f, double a, double fa, double m, double fm, double b, double fb, double tol);


//~ // vector array manipulation
//~ void readArray(int i, double vecArr[][3], double vec[3]);
//~ void assignArray(int i, double vec[3], double vecArr[][3]);

//~ double mean(int N, double arr[]);
//~ double standdev(int N, double arr[]);

// integrals
//~ double I_010(const double R, const double zDiff, const double cylP[3]);
//~ double I_011(const double R, const double zDiff, const double cylP[3]);
//~ double HeumansLambda(const double beta, const double k);

//~ // elliptical integrals
//~ double complEllipK1(double x2);
//~ double complEllipK2(double x2);
//~ double complEllipK3(double x2);
//~ double complEllipK4(double x2);
//~ double complEllipE1(double x2);
//~ double complEllipE2(double x2);
//~ double complEllipE3(double x2);



#endif
