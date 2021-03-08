#include "Loop.h"

// Constructor
Loop::Loop(const double R, const double I, const double x, const double y, const double z): R(R), I(I), x(x), y(y), z(z){}
Loop::Loop():R((0.0463701-0.04125)/2.0+0.04125),I(120*600),x(0),y(0),z(0){}
// Get Methods

double Loop::getR() const {
		return R;
}

double Loop::getI() const {
		return I;
}

double Loop::getx() const {
		return x;
}

double Loop::gety() const {
		return y;
}

double Loop::getz() const {
		return z;
}

void Loop::print(){
	std::cout << "Loop: " << "R = " << R << "; I = " << I << "; x = " << x << "; y = " << y << "; z = " << z << std::endl;
}
