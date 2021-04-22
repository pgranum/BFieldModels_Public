#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <vector>

//~ #include <cmath>
//~ #include <math.h>
//~ #include "mathFuncs.h"
//~ #include "Loop.h"
//~ #include "Generic3Vector.hpp"
//~ #include <tr1/cmath>
//~ #include <tgmath.h>
//~ #include <chrono>
//~ #include "coorTransf.h"
//~ #include "writeToFile.h"
//~ #include "ALPHAPhysicalConstants.hpp"

// print functions
//~ void print(double scal, std::string scalarName);
//~ void print(int scal, std::string scalarName);
//~ void printVec(double vec[3]);
//~ void printVec(long double vec[3]);

void printVec(const double vec[3], const std::string vecName);

//~ void printVec(long double vec[3], std::string vecName);
//~ void printArr(double arr[], int N);
//~ void printArr(long double arr[], int N);
//~ void printArr(int arr[], int N, std::string vecArrName);
//~ void printArr(double arr[], int N, std::string vecArrName);
//~ void printArr(long double arr[], int N, std::string vecName);
//~ void printStdVec(std::vector<int> vec, std::string vecName);
//~ void printStdVec(std::vector<double> vec, std::string vecName);
//~ void printStdVec(std::vector<long double> vec, std::string vecName);
//~ void printVecArr(double vecArr[][3], int N);
//~ void printVecArr(double vecArr[][3], int N, std::string vecArrName);

class writeToFile{
private:
	std::ofstream myFile;
public:	

	
protected:
public:
	writeToFile(std::string basePath, std::string fileName){
		basePath.append("/");
		basePath.append(fileName);
		basePath.append(".bin");
		myFile.open(basePath, std::ios::binary | std::ios::out | std::ios::app); // appends to binary file
	};
	
	~writeToFile(){
		myFile.close(); 
	};
	
	template <typename T>
	int write(T data){
		myFile.write(reinterpret_cast <const char*> (&data) , sizeof(data));
		return 0;
	}
	
	template<typename T>
	int write(T arr[], int n){
		myFile.write(reinterpret_cast <const char*> (arr) , sizeof(T)*n);
		return 0;
	}
	
	template<typename T>
	int write(std::vector<T> arr){
		myFile.write(reinterpret_cast <const char*> (arr.data()), sizeof(T)*arr.size());
		return 0;
	}

	
};

//~ class writeVecArrToFile{
//~ private:
	//~ writeToFile* x;
	//~ writeToFile* y;
	//~ writeToFile* z;

//~ public:
	//~ writeVecArrToFile(std::string basePath, std::string fileName){
		//~ x = new writeToFile(basePath,fileName+"_x");
		//~ y = new writeToFile(basePath,fileName+"_y");
		//~ z = new writeToFile(basePath,fileName+"_z");
	//~ };
	
	//~ ~writeVecArrToFile(){
		//~ delete x;
		//~ delete y;
		//~ delete z;
	//~ }
	
	//~ void write(double r[][3], int N){
		//~ for(int i=0;i<N;i++){
			//~ x -> write(r[i][0]);
			//~ y -> write(r[i][1]);
			//~ z -> write(r[i][2]);
		//~ }
	//~ }
	
	//~ void write(std::vector<CartesianVector> r){
		//~ int N = r.size();
		//~ for(int i=0;i<N;i++){
			//~ x -> write(r.at(i).GetX());
			//~ y -> write(r.at(i).GetY());
			//~ z -> write(r.at(i).GetZ());
		//~ }
	//~ }
	
//~ };

class writeBFieldToFile{
private:
	writeToFile* x;
	writeToFile* y;
	writeToFile* z;
	writeToFile* Bx;
	writeToFile* By;
	writeToFile* Bz;

public:
	writeBFieldToFile(std::string basePath, std::string fileName){
		x = new writeToFile(basePath,fileName+"_x");
		y = new writeToFile(basePath,fileName+"_y");
		z = new writeToFile(basePath,fileName+"_z");
		Bx = new writeToFile(basePath,fileName+"_Bx");
		By = new writeToFile(basePath,fileName+"_By");
		Bz = new writeToFile(basePath,fileName+"_Bz");		
	};
	
	~writeBFieldToFile(){
		delete x;
		delete y;
		delete z;
		delete Bx;
		delete By;
		delete Bz;
	}
	
	void write(double r[][3], double B[][3], int N){
		for(int i=0;i<N;i++){
			x -> write(r[i][0]);
			y -> write(r[i][1]);
			z -> write(r[i][2]);
			Bx -> write(B[i][0]);
			By -> write(B[i][1]);
			Bz -> write(B[i][2]);
		}
	}	
};


#endif
