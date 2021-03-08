#ifndef TUBE_H
#define TUBE_H

#include <iostream>

/** 
 * Implementation of a finite solenoid
 *
 * The finite solenoid contains information about its two radii, length, number of wires, the current in the wires, and the posistion of the center of the tube 
 */
class Tube {
	
private: 
// Class Variables
const double R1; // The radius of the inside of the tube
const double R2; // The radius of the outside of the tube
const int	 N; // The number of wires in the tube
const double i; // The current in each wire
const double L; // The length of the tube
const double x; // The the x position of the center of the tube
const double y; // The the y position of the center of the tube
const double z; // The the z position of the center of the tube
	
public:
// Class Constructor
Tube(const double R1, const double R2, const int N,const double i, const double L, const double x, const double y, const double z);
Tube(const double R1, const double R2, const double I, const double L, const double x, const double y, const double z);
Tube();
	
		// Get Methods
	
// Get Methods
/** This function return the inner radius
* 
* @return the inner radius
*/
double getR1() const;
/** This function return the outer radius
* 
* @return the outer radius
*/
double getR2() const;
/** This function return the number of wires
* 
* @return the number of wires
*/
int getN() const;
/** This function return the total current
* 
* @return the total current
*/
double getI() const;
/** This function return the total current
* 
* @return the total current
*/
double geti() const; 
/** This function return the length
* 
* @return the length
*/
double getL() const;
/** This function return the x coordinate the center
* 
* @return the current of a Loop object
*/
double getx() const;
/** This function return the y coordinate of the center
* 
* @return the y coordinate of the center
*/
double gety() const;
/** This function return the z coordinate of the center
* 
* @return the z coordinate of the center
*/
double getz() const;

void print();

};


#endif // TUBE_H
