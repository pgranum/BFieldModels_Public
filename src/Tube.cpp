#include "Tube.h"

// Constructor
Tube::Tube(double R1, double R2, int N, double i, double L, double x, double y, double z): R1(R1), R2(R2),N(N), i(i), L(L), x(x), y(y), z(z){}
Tube::Tube(double R1, double R2, double I, double L, double x, double y, double z): R1(R1), R2(R2), N(1),i(I), L(L), x(x), y(y), z(z){}
Tube::Tube():R1(0.04125),R2(0.04637),N(4*30),i(600),L(0.03468),x(0.),y(0.),z(0.){}

// Get Methods
double Tube::getR1() const{
		return R1;
}
double Tube::getR2() const{
		return R2;
}
int Tube::getN() const{
		return N;
}
double Tube::geti() const{
		return i;
}
double Tube::getI() const{
		return N*i;
}
double Tube::getL() const{
		return L;
}
double Tube::getx() const{
		return x;
}
double Tube::gety() const{
		return y;
}
double Tube::getz() const{
		return z;
}
void Tube::print() {
	std::cout << "Tube: " << "R1 = " << R1 << "; R2 = " << R2 << "; N = " << N << "; i = " << i << "; x = " << x << "; y = " << y << "; z = " << z <<std:: endl;
}


