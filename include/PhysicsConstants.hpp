#ifndef PHYSICSCONSTANTS_H
#define PHYSICSCONSTANTS_H
#include <math.h>

namespace PhysicsConstants {
  static constexpr double pi = M_PI;// double(2.0)*asin(double(1.0)); 
  //~ static constexpr double pi_2 = pi*pi;
  //~ static constexpr double kboltz = 1.3806488e-23;
  //~ static constexpr double mprot = 1.672621777e-27;
  //~ static constexpr double kboltz_2_mprot = 2.0 * kboltz / mprot;
  //~ static constexpr double magmom0 = 9.27400968e-24;
  static constexpr double mu0 = 4*pi*1e-7;
};

#endif
