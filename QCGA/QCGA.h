#pragma once
#include <string>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <map>
#include <iomanip>
#include <utility>
#include <string_view>

constexpr size_t TOTAL_DIMENSION{32'768};  //total dimension of an algebra (2^5)
constexpr size_t GENERATING_BASIS_SIZE{15};  //dimension of generating space R^(9,6)
constexpr size_t ALGEBRA_P{9};  //number of positive squared vectors
constexpr size_t ALGEBRA_Q{6};  //number of negative squared vectors
constexpr long long PRECISION{1'000'000'000'000};  //constant for rounding

#define zero_vector (QCGA()) //zero vecor
#define one (QCGA::generatingBlades[0]) //scalar

#define	e1 (QCGA::generatingBlades[1]) //euclidean vectors
#define	e2 (QCGA::generatingBlades[2])
#define	e3 (QCGA::generatingBlades[3])

#define	e4 (QCGA::generatingBlades[4])//e+: positive squared vectors
#define	e5 (QCGA::generatingBlades[5])
#define	e6 (QCGA::generatingBlades[6])
#define	e7 (QCGA::generatingBlades[7])
#define	e8 (QCGA::generatingBlades[8])
#define	e9 (QCGA::generatingBlades[9])

#define	e10 (QCGA::generatingBlades[10])//e-: negative squared vectors
#define	e11 (QCGA::generatingBlades[11])
#define	e12 (QCGA::generatingBlades[12])
#define	e13 (QCGA::generatingBlades[13])
#define	e14 (QCGA::generatingBlades[14])
#define	e15 (QCGA::generatingBlades[15])

#define	eo6 (0.5*(QCGA::generatingBlades[15]-QCGA::generatingBlades[9]))//change of basis
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

#define rxy ((e1 ^ e2) + (eo6 ^ ei5) + (ei6 ^ eo5) + (2 * (eo4 ^ ei2)) + (2 * (ei4 ^ eo2)) + (eo4 ^ ei3)) //generator for rotation in the xy-plane
#define rxz ((-1 * (e3 ^ e1)) + (-1 * (ei4 ^ eo6)) + (-2 * (ei3 ^ eo5)) + (-1 * (ei2 ^ eo5)) + (-1 * (eo4 ^ ei6)) + (-2 * (eo3 ^ ei5))) //generator for rotation in the xz-plane
#define ryz ((e2 ^ e3) + (eo6 ^ ei3) + (ei5 ^ eo4) + (2 * (ei6 ^ eo3)) + (2 * (eo2 ^ ei6)) + (ei2 ^ eo6) + (eo5 ^ ei4)) //generator for rotation in the yz-plane

#define T1x (one - 0.5 * distance * (e1 ^ ei1))
#define T2x (one - 0.5 * distance * (e1 ^ ei2) + 0.25 * pow(distance, 2) * (ei1 ^ ei2))
#define T3x (one - 0.5 * distance * (e1 ^ ei3) + 0.25 * pow(distance, 2) * (ei1 ^ ei3) + 0.25 * pow(distance, 2) * (ei2 ^ ei3))
#define T4x (one - 0.5 * distance * (e2 ^ ei4))
#define T5x (one - 0.5 * distance * (e3 ^ ei5))
#define Tx (T1x*T2x*T3x*T4x*T5x) //Translator in x

#define T1y (one - 0.5 * distance * (e2 ^ ei1))
#define T2y (one + 0.5 * distance * (e2 ^ ei2) - 0.25 * pow(distance, 2) * (ei1 ^ ei2))
#define T3y (one - 0.5 * distance * (e1 ^ ei4))
#define T4y (one - 0.5 * distance * (e3 ^ ei6))
#define Ty (T1y*T2y*T3y*T4y) //Translator in y

#define T1z (one - 0.5 * distance * (e3 ^ ei1))
#define T2z	(one + 0.5 * distance * (e3 ^ ei3) - 0.25 * pow(distance, 2) * (ei1 ^ ei3))
#define T3z	(one - 0.5 * distance * (e1 ^ ei5))
#define T4z	(one - 0.5 * distance * (e2 ^ ei6))
#define Tz (T1z*T2z*T3z*T4z) //Translator in z

#define I Blade(e1*e2*e3*e4*e5*e6*e7*e8*e9*e10*e11*e12*e13*e14*e15) //Pseaudoscalar

enum rotation_planes
{
	xy, xz, yz
};
enum translation_directions
{
	x, y, z
};

class QCGA
{
public:
	static QCGA generatingBlades[]; //stores 1,e1,e2,...,en. Will be made protected after all is done. I want defined vectors to work in Blade now
protected:
	std::map<std::string, long double> m_mapLabelToCoefficient; //representation of a general multivector
	// 3 + 2e1 - e1*e2*e3
	// ==================
	// 1 |  e1 | e1*e2*e3 
	// 3 |   2 |       -1
public:
	static void generateGeneratingBlades();
	QCGA(); //creates zero vector;
	QCGA(const std::string& input); 
	QCGA(std::string&& input) noexcept;
	QCGA(const std::map<std::string, long double>& map);
	QCGA(std::map<std::string, long double>&& map);
	QCGA(const std::pair<std::string,long double>& basis_blade); 
	QCGA(std::pair<std::string,long double>&& basis_blade);

	QCGA(const QCGA& instance); 
	QCGA(QCGA&& instance) noexcept; 
	virtual ~QCGA() = default;

	long double toNumeric(); //returs coefficient at basis blade "1"

	const std::map<std::string, long double>& getSTDmapLabelToCoefficient() const; //returns map (=representation of multivector)

	//**********************************OPERATORS**********************************\\
	
	QCGA rotorExponential(unsigned int degree, long double phi) const;
	QCGA translatorExponential(unsigned int degree, long double distance) const;
	QCGA bivectorExponential(unsigned int degree, long double parameter) const;
	QCGA& operator=(const QCGA& other); 
	QCGA& operator=(QCGA&& other) noexcept; 
	bool operator==(const QCGA& other) const; //equals operator
	bool operator!=(const QCGA& other) const; //not equals operator
	QCGA operator[](int grade) const; //grade projection
	QCGA operator[](const QCGA& other) const; //basis blade selection
	QCGA operator*(const QCGA& other) const; //geometric product operator
	QCGA operator*(QCGA&& other) const; 
	QCGA operator*(long double scalar) const; //multiplying by scalar from the right operator
	QCGA operator~() const; //reverse operator
	QCGA operator+(const QCGA& other) const; //multivector addition operator
	QCGA operator-(const QCGA& other) const; //multivector substraction operator
	QCGA operator|(const QCGA& other) const; //inner product operator
	QCGA operator^(const QCGA& other) const; //outer product operator
	QCGA operator^(int exponent) const; //exponent operator
	QCGA operator/(long double divider) const; //dividing operator
	QCGA scalarProduct(const QCGA& b) const;

	static QCGA rotate(const QCGA& point, rotation_planes plane, long double angle);
	static QCGA translate(const QCGA& point, translation_directions plane, long double angle);

	int grade(std::string_view label) const; //returns grade of basis blade (if we give it appropriate label...)
	std::string log() const; //returns multivector, used in << operator
protected:
	//**********************************STATIC_SUPPORT_FUNCTIONS**********************************\\

	static int calculateSign(const std::vector<int>& permutation); //Helps in validating basis Element, calculates sign of permutation
	static void simplifyBasisBlade(std::string& label, int& sign); //simplifies label in a form of for example  e1e2e3e2e3 into e1
	static void processVector(std::vector<int>& vec, int& sign); //used when simplifying results of geometric product: e1e2e5e2e3e4e5 -> e1e5e3e4e5 -> e1e3e4 represented byjust numbers (1252345 -> 15345 ...)
	static std::vector<int> extractIntegersFromBasisBlades(std::string_view label); //from a given label, for example e1*e2*e3, returns vector {1,2,3}

	//**********************************STATIC_SUPPORT_OPERATORS**********************************\\

	QCGA operator||(const QCGA& other) const; //inner product of two basis blades operator
	QCGA operator &&(const QCGA& other) const; //outer product of two basis blades operator
	QCGA operator ()(int grade) const; //grade projection of basis blade operator

	//**********************************ACTUAL_ATRIBUTES**********************************\\

	

	//**********************************SUPPORT_FUNCTIONS**********************************\\

	void deleteZeroFromVector(); //if multivector is of a form 0*1 + c1e1+ c2e1*e2 +... it removes 0*1
};
//**********************************NON-MEMBER_OPERATORS**********************************\\

QCGA operator*(const long double& scalar, const QCGA& onther); //multiplying by scalar from the left
std::ostream& operator<<(std::ostream& stream, const QCGA& vector); //operator for printing
std::vector<QCGA> makeQCGAFromBasisBlades(const QCGA& multivector); //returns vector of basis blades in linear combination of general multivector
