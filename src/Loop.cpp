#include "Loop.h"

// Constructor
Loop::Loop(const double R, const double I, const double x, const double y, const double z): R(R), I(I), x(x), y(y), z(z),NSegments(0){}
Loop::Loop(const double R, const double I, const double x, const double y, const double z, const int NSegments): R(R), I(I), x(x), y(y), z(z), NSegments(NSegments){	
	// setting up the straght line segments used by the Biot-Savart method
	const double dPhi = 2*PhysicsConstants::pi/(NSegments-1);
	for(int i=0; i<NSegments; i++){
		const double phi = i*dPhi;
		xs.push_back(cos(phi)*R);
		ys.push_back(sin(phi)*R);
	}
}
Loop::Loop():R((0.0463701-0.04125)/2.0+0.04125),I(120*600),x(0),y(0),z(0),NSegments(0){}
// Get Methods

double Loop::getR() const {
		return R;
}

double Loop::getRR() const {
		return R*R;
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

double Loop::getxi(const int i) const {
		return xs[i];
}

double Loop::getyi(const int i) const {
		return ys[i];
}

void Loop::print(){
	std::cout << "Loop: " << "R = " << R << "; I = " << I << "; x = " << x << "; y = " << y << "; z = " << z << std::endl;
}
