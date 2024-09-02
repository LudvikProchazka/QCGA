#ifndef __QCGA_H__
#define __QCGA_H__

#include <string>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <map>



#define TOTAL_DIMENSION 32768 //total dimension of algebra (2^5)
#define GENERATING_BASIS_DIMENSION 15 //dimension of generating space R^5
#define ALGEBRA_P 9 //number of positive squared vectors
#define ALGEBRA_Q 6 //number of negative squared vectors
#define PRECISION 100000 //for rounding

#define zero (QCGA())
#define one (QCGA::generatingBlades[0])

#define	e1 (QCGA::generatingBlades[1])
#define	e2 (QCGA::generatingBlades[2])
#define	e3 (QCGA::generatingBlades[3])

#define	e4 (QCGA::generatingBlades[4])//e+
#define	e5 (QCGA::generatingBlades[5])
#define	e6 (QCGA::generatingBlades[6])
#define	e7 (QCGA::generatingBlades[7])
#define	e8 (QCGA::generatingBlades[8])
#define	e9 (QCGA::generatingBlades[9])

#define	e10 (QCGA::generatingBlades[10])//e-
#define	e11 (QCGA::generatingBlades[11])
#define	e12 (QCGA::generatingBlades[12])
#define	e13 (QCGA::generatingBlades[13])
#define	e14 (QCGA::generatingBlades[14])
#define	e15 (QCGA::generatingBlades[15])

#define	eo6 (0.5*(QCGA::generatingBlades[15]-QCGA::generatingBlades[9])) 
#define	eo5 (0.5*(QCGA::generatingBlades[14]-QCGA::generatingBlades[8])) 
#define	eo4 (0.5*(QCGA::generatingBlades[13]-QCGA::generatingBlades[7])) 
#define	eo3 (0.5*(QCGA::generatingBlades[12]-QCGA::generatingBlades[6])) 
#define	eo2 (0.5*(QCGA::generatingBlades[11]-QCGA::generatingBlades[5])) 
#define	eo1 (0.5*(QCGA::generatingBlades[10]-QCGA::generatingBlades[4])) 
#define	ei1 (QCGA::generatingBlades[10]+QCGA::generatingBlades[4]) 
#define	ei2 (QCGA::generatingBlades[11]+QCGA::generatingBlades[5]) 
#define	ei3 (QCGA::generatingBlades[12]+QCGA::generatingBlades[6]) 
#define	ei4 (QCGA::generatingBlades[13]+QCGA::generatingBlades[7]) 
#define	ei5 (QCGA::generatingBlades[14]+QCGA::generatingBlades[8]) 
#define	ei6 (QCGA::generatingBlades[15]+QCGA::generatingBlades[9]) 


#define I Blade(e1*e2*e3*e4*e5*e6*e7*e8*e9*e10*e11*e12*e13*e14*e15) //Pseaudoscalar of an algebra

class QCGA
{
public:
	static QCGA generatingBlades[]; //stores 1,e1,e2,...,en
	static void generateGeneratingBlades(); //generates generatingBlades

	explicit QCGA(std::string input); //constructor used for construct generatingBlades
	explicit QCGA(std::map<std::string, double> map);
	QCGA(const QCGA& Multivector); //instanciate CGA object from given objects
	QCGA(); //default constructor, calls CGA("0");

	std::map<std::string, double> getSTDmapLabelToCoefficient() const; //returns map (=representation of multivector)
	double toNumeric(); //returs coefficient at basis blade "1"


	//**********************************OPERATORS**********************************\\

	bool operator==(const QCGA& other) const; //equals operator
	bool operator!=(const QCGA& other) const; //not equals operator
	QCGA operator[](const int grade) const; //grade projection
	QCGA operator*(const QCGA& other) const; //geometric product operator
	QCGA operator*(const double scalar) const; //multiplying by scalar from the right operator
	QCGA operator~() const; //reverse operator
	QCGA operator+(const QCGA& other) const; //multivector addition operator
	QCGA operator-(const QCGA& other) const; //multivector substraction operator
	QCGA operator|(const QCGA& other) const; //inner product operator
	QCGA operator^(const QCGA& other) const; //outer product operator
	QCGA operator^(const int exponent) const; //exponent operator
	QCGA operator/(const double divider) const; //dividing operator

	int grade(std::string label) const; //returns grade of basis blade (if we give it appropriate label...)
	std::string log() const; //returns multivector, used in << operator
protected:
	//**********************************STATIC_SUPPORT_FUNCTIONS**********************************\\

	static int calculateSign(const std::vector<int>& permutation); //Helps in validating basis Element, calculates sign of permutation
	static bool searchString(std::string string[], int size, std::string target); //returns bool depending if target string is present in string array
	static void simplifyBasisBlade(std::string& label, int& sign); //simplifies label in a form of for example  e1e2e3e2e3 into e1
	static void processVector(std::vector<int>& vec, int& sign); //used when simplifying results of geometric product: e1e2e5e2e3e4e5 -> e1e5e3e4e5 -> e1e3e4 represented byjust numbers (1252345 -> 15345 ...)
	static std::vector<int> extractIntegersFromBasisBlades(std::string label); //from a given label, for example e1*e2*e3, returns vector {1,2,3}

	//**********************************STATIC_SUPPORT_OPERATORS**********************************\\

	QCGA operator||(const QCGA& other) const; //inner product of two basis blades operator
	QCGA operator &&(const QCGA& other) const; //outer product of two basis blades operator
	QCGA operator ()(const int grade) const; //grade projection of basis blade operator

	//**********************************ACTUAL_ATRIBUTES**********************************\\

	std::map<std::string, double> STDmapLabelToCoefficient; //representation of a general multivector

	//**********************************SUPPORT_FUNCTIONS**********************************\\

	void deleteZeroFromVector(); //if multivector is of a form 0*1 + c1e1+ c2e1*e2 +... it removes 0*1
};
//**********************************NON-MEMBER_OPERATORS**********************************\\

QCGA operator*(const double scalar, const QCGA& onther); //multiplying by scalar from the left
std::ostream& operator<<(std::ostream& stream, const QCGA& vector); //operator for printing
std::vector<QCGA> makeCGAFromBasisBlades(const QCGA& multivector); //returns vector of basis blades in linear combination of general multivector
void removeOccurences(std::string& str, const std::string substr); //removes occurences of substring in string


#endif