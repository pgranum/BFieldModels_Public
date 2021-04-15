#include "Shell.h"

// Constructor
Shell::Shell(const double R, const int N, const double i, const double L, const double x, const double y, const double z) : R(R), N(N),i(i),L(L),x(x),y(y),z(z),Z1(z-L*0.5),Z2(z+L*0.5) {}
Shell::Shell(const double R, const double I, const double L, const double x, const double y, const double z) : R(R), N(1),i(I),L(L),x(x),y(y),z(z),Z1(z-L*0.5),Z2(z+L*0.5) {}
Shell::Shell() : R((0.0463701-0.04125)/2.0+0.04125),N(120),i(600),L(0.0346811),x(0),y(0),z(0),Z1(z-L*0.5),Z2(z+L*0.5){} // Specifications for A2 where the radius has been set to the centre of the mirror

// Get Methods
double Shell::getR() const {
		return R;
}
double Shell::getRR() const {
		return R*R;
}
int    Shell::getN() const {
		return N;
}
double Shell::geti() const {
		return i;
}
double Shell::getI() const {
		return i*N;
}
double Shell::getL() const {
		return L;
}
double Shell::getx() const {
		return x;
}
double Shell::gety() const {
		return y;
}
double Shell::getz() const {
		return z;
}
void Shell::print(){
	std::cout << "Tube: " << " R = " << R << "; N = " << N << "; i = " << i << "; x = " << x << "; y = " << y << "; z = " << z << std::endl;
}
