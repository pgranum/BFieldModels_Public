#ifndef LOOP_H
#define LOOP_H

#include <iostream>
#include <vector>

#include "PhysicsConstants.hpp"

/** 
 * Implementation of a current loop
 *
 * The current loop contains information about its radius current and the posistion of the center of the loop 
 */
class Loop{
	
private: 
// Class Variables
const double R; // The radius of the loop
const double I; // The current in the loop
const double x; // The the x position of the center of the loop
const double y; // The the y position of the center of the loop
const double z; // The the z position of the center of the loop
	
public:
// Class Constructor

Loop(const double R, const double I, const double x, const double y, const double z);
Loop(const double R, const double I, const double x, const double y, const double z, const int NSegments);
Loop();


void print();

// Get Methods
/** This function return the radius of a Loop object
* 
* @return the radius of a Loop object
*/
double getR() const;
/** This function return the radius of a Loop object
* 
* @return the radius of a Loop object
*/
double getRR() const;
/** This function return the radius of a Loop object
* 
* @return the radius of a Loop object squared
*/
double getI() const;
/** This function return the current of a Loop object
* 
* @return the current of a Loop object
*/
double geti() const; 
/** This function return the x coordinate the center of a Loop object
* 
* @return the current of a Loop object
*/
double getx() const;
/** This function return the y coordinate of the center of a Loop object
* 
* @return the y coordinate of the center of a Loop object
*/
double gety() const;
/** This function return the z coordinate of the center of a Loop object
* 
* @return the z coordinate of the center of a Loop object
*/
double getz() const;

double getxi(const int i) const;

double getyi(const int i) const;

};

#endif // LOOP_H
