#ifndef SHELL_H
#define SHELL_H

#include <iostream>

/** 
 * Implementation of a thin shell solenoid
 *
 * The thin shell solenoid contains information about its radius, length, number of wires, the current in the wires, and the posistion of the center of the tube 
 */
class Shell {
	
protected: 
// Class Variables
const double R; // The radius of the inside of the shell
const int	 N; // The number of wires in the shell
const double i; // The current in each wire
const double L; // The length of the shell
const double x; // The the x position of the center of the shell
const double y; // The the y position of the center of the shell
const double z; // The the z position of the center of the shell
	
public:
// Class Constructor
Shell(const double R, const int N,const double i, const double L, const double x, const double y, const double z);
Shell(const double R, const double I, const double L, const double x, const double y, const double z);
Shell();
	
// Get Methods
/** This function return the radius
* 
* @return the radius
*/
double getR() const;
/** This function return the radius
* 
* @return the radius
*/
double getRR() const;
/** This function return the radius
* 
* @return the radius squared
*/
int getN() const;
/** This function return the  number of wires
* 
* @return the number of wires
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


#endif
