#include "utils.h"

// print functions

//~ void printVec(double vec[3]){
//~ /* This function print the entries of a vector
//~ * 
//~ * @param vec 		the vector being printed
//~ */

	//~ std::cout << " [" << vec[0]
		//~ << ", " << vec[1]
		//~ << ", " << vec[2]
		//~ << "] " << std::endl;
//~ }
//~ void printVec(long double vec[3]){
//~ /* This function print the entries of a vector
//~ * 
//~ * @param vec 		the vector being printed
//~ */

	//~ std::cout << " [" << vec[0]
		//~ << ", " << vec[1]
		//~ << ", " << vec[2]
		//~ << "] " << std::endl;
//~ }

void printVec(const double vec[3], const std::string vecName){
/* This function print the entries of a vector with its name
* 
* @param vec 		the vector being printed
* @param vecName 	the name of the vector
*/

	std::cout << vecName << " = [" << vec[0]
		<< ", " << vec[1]
		<< ", " << vec[2]
		<< "] " << std::endl;
}

//~ void printVec(long double vec[3], std::string vecName){
//~ /* This function print the entries of a vector with its name
//~ * 
//~ * @param vec 		the vector being printed
//~ * @param vecName 	the name of the vector
//~ */

	//~ std::cout << vecName << " = [" << vec[0]
		//~ << ", " << vec[1]
		//~ << ", " << vec[2]
		//~ << "] " << std::endl;
//~ }

//~ void printArr(double arr[], int N){
//~ /* This function print the N entries of an array
//~ * 
//~ * @param arr		the array being printed
//~ * @param N 			the length of the array
//~ */
	//~ std::cout << "[";
	//~ for(int i = 0; i<N; i++){
		//~ std::cout << arr[i];
		//~ if(i < N-1){
			//~ std::cout << ", ";
		//~ }
	//~ }	
	//~ std::cout << "]" << std::endl;
//~ }
//~ void printArr(long double arr[], int N){
//~ /* This function print the N entries of an array
//~ * 
//~ * @param arr		the array being printed
//~ * @param N 			the length of the array
//~ */
	//~ std::cout << "[";
	//~ for(int i = 0; i<N; i++){
		//~ std::cout << arr[i];
		//~ if(i < N-1){
			//~ std::cout << ", ";
		//~ }
	//~ }	
	//~ std::cout << "]" << std::endl;
//~ }
//~ void printArr(int arr[], int N, std::string vecName){
//~ /* This function print the N entries of an array with its name in front
//~ * 
//~ * @param arr		the array being printed
//~ * @param N 			the length of the array
//~ * @param vecName 	the name of the vector
//~ */
	//~ std::cout << vecName << " = [";
	//~ for(int i = 0; i<N; i++){
		//~ std::cout << arr[i];
		//~ if(i < N-1){
			//~ std::cout << ", ";
		//~ }
	//~ }	
	//~ std::cout << "]" << std::endl;
//~ }
void printArr(const double arr[], const int N, std::string vecName){
/* This function print the N entries of an array with its name in front
* 
* @param arr		the array being printed
* @param N 			the length of the array
* @param vecName 	the name of the vector
*/
	std::cout << vecName << " = [";
	for(int i = 0; i<N; i++){
		std::cout << arr[i];
		if(i < N-1){
			std::cout << ", ";
		}
	}	
	std::cout << "]" << std::endl;
}
//~ void printArr(long double arr[], int N, std::string vecName){
//~ /* This function print the N entries of an array with its name in front
//~ * 
//~ * @param arr		the array being printed
//~ * @param N 			the length of the array
//~ * @param vecName 	the name of the vector
//~ */
	//~ std::cout << vecName << " = [";
	//~ for(int i = 0; i<N; i++){
		//~ std::cout << arr[i];
		//~ if(i < N-1){
			//~ std::cout << ", ";
		//~ }
	//~ }	
	//~ std::cout << "]" << std::endl;
//~ }
//~ void printStdVec(std::vector<int> vec, std::string vecName){
//~ /* This function print the N entries of a std:vector with its name in front
//~ * 
//~ * @param vec		the vector being printed
//~ * @param vecName 	the name of the vector
//~ */
	//~ std::cout << vecName << " = [";
	//~ for(unsigned int i = 0; i<vec.size(); i++){
		//~ std::cout << vec[i];
		//~ if(i < vec.size()-1){
			//~ std::cout << ", ";
		//~ }
	//~ }	
	//~ std::cout << "]" << std::endl;
//~ }
//~ void printStdVec(std::vector<double> vec, std::string vecName){
//~ /* This function print the N entries of a std:vector with its name in front
//~ * 
//~ * @param vec		the vector being printed
//~ * @param vecName 	the name of the vector
//~ */
	//~ std::cout << vecName << " = [";
	//~ for(unsigned int i = 0; i<vec.size(); i++){
		//~ std::cout << vec[i];
		//~ if(i < vec.size()-1){
			//~ std::cout << ", ";
		//~ }
	//~ }	
	//~ std::cout << "]" << std::endl;
//~ }

//~ void printStdVec(std::vector<long double> vec, std::string vecName){
//~ /* This function print the N entries of a std:vector with its name in front
//~ * 
//~ * @param vec		the vector being printed
//~ * @param vecName 	the name of the vector
//~ */
	//~ std::cout << vecName << " = [";
	//~ for(unsigned int i = 0; i<vec.size(); i++){
		//~ std::cout << vec[i];
		//~ if(i < vec.size()-1){
			//~ std::cout << ", ";
		//~ }
	//~ }	
	//~ std::cout << "]" << std::endl;
//~ }

//~ void printVecArr(double vecArr[][3], int N){
//~ /* This function print the N vectors of a vector array with its name in front
//~ * 
//~ * @param vecArr 		the vector array being printed
//~ * @param N 				the length of the vector array
//~ */
		//~ std::cout << "Vector Array size "<< N << "x"
			  //~ << 3
			  //~ << " = "
			  //~ << std::endl;
			  
	//~ for(int i=0; i<N; i++){
		//~ std::cout << vecArr[i][0]
				  //~ << ", " << vecArr[i][1]
			      //~ << ", " << vecArr[i][2]
			      //~ <<  std::endl;
	//~ }

//~ }
//~ void printVecArr(double vecArr[][3], int N, std::string vecArrName){
/* This function print the N vectors of a vector array
* 
* @param vecArr 		the vector array being printed
* @param N 				the length of the vector array
* @param vecArrName 	the name of the vector array
*/
		//~ std::cout << vecArrName << " "
			  //~ << N << "x"
			  //~ << 3
			  //~ << " = "
			  //~ << std::endl;
			  
	//~ for(int i=0; i<N; i++){
		//~ std::cout << vecArr[i][0]
				  //~ << ", " << vecArr[i][1]
			      //~ << ", " << vecArr[i][2]
			      //~ <<  std::endl;
	//~ }

//~ }

