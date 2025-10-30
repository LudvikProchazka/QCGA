#pragma once
#include <string>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <map>

constexpr size_t TOTAL_DIMENSION{32'768};				//total dimension of an algebra (2^5)
constexpr size_t GENERATING_BASIS_SIZE{15};				//dimension of generating space R^(9,6)
constexpr size_t ALGEBRA_P{9};							//number of positive squared vectors
constexpr size_t ALGEBRA_Q{6};							//number of negative squared vectors
constexpr long double PRECISION{1'000'000'000'000.0};	//constant for rounding

#define zero_vector (GAQ())			//zero vecor
#define one (GAQ::generatingBlades[0]) //scalar

#define	e1 (GAQ::generatingBlades[1]) //euclidean vectors
#define	e2 (GAQ::generatingBlades[2])
#define	e3 (GAQ::generatingBlades[3])

#define	e4 (GAQ::generatingBlades[4]) //e+: positive squared vectors
#define	e5 (GAQ::generatingBlades[5])
#define	e6 (GAQ::generatingBlades[6])
#define	e7 (GAQ::generatingBlades[7])
#define	e8 (GAQ::generatingBlades[8])
#define	e9 (GAQ::generatingBlades[9])

#define	e10 (GAQ::generatingBlades[10]) //e-: negative squared vectors
#define	e11 (GAQ::generatingBlades[11])
#define	e12 (GAQ::generatingBlades[12])
#define	e13 (GAQ::generatingBlades[13])
#define	e14 (GAQ::generatingBlades[14])
#define	e15 (GAQ::generatingBlades[15])

#define	eo6 (0.5*(GAQ::generatingBlades[15]-GAQ::generatingBlades[9])) //change of basis
#define	eo5 (0.5*(GAQ::generatingBlades[14]-GAQ::generatingBlades[8])) 
#define	eo4 (0.5*(GAQ::generatingBlades[13]-GAQ::generatingBlades[7])) 
#define	eo3 (0.5*(GAQ::generatingBlades[12]-GAQ::generatingBlades[6])) 
#define	eo2 (0.5*(GAQ::generatingBlades[11]-GAQ::generatingBlades[5])) 
#define	eo1 (0.5*(GAQ::generatingBlades[10]-GAQ::generatingBlades[4])) 
#define	ei1 (GAQ::generatingBlades[10]+GAQ::generatingBlades[4]) 
#define	ei2 (GAQ::generatingBlades[11]+GAQ::generatingBlades[5]) 
#define	ei3 (GAQ::generatingBlades[12]+GAQ::generatingBlades[6]) 
#define	ei4 (GAQ::generatingBlades[13]+GAQ::generatingBlades[7]) 
#define	ei5 (GAQ::generatingBlades[14]+GAQ::generatingBlades[8]) 
#define	ei6 (GAQ::generatingBlades[15]+GAQ::generatingBlades[9]) 

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

class GAQ
{
public:
	static void GenerateGeneratingBlades();
	GAQ(); //creates zero vector;
	GAQ(const std::string& input); 
	GAQ(std::string&& input) noexcept;
	GAQ(const std::map<std::string, long double>& map);
	GAQ(std::map<std::string, long double>&& map);
	GAQ(const std::pair<std::string,long double>& basisBlade); 
	GAQ(std::pair<std::string,long double>&& basisBlade);

	GAQ(const GAQ& instance); 
	GAQ(GAQ&& instance) noexcept; 
	virtual ~GAQ() = default;

	long double ToNumeric(); //returs coefficient at basis blade "1"
	bool IsEqual(const GAQ& second, double precision) const;

	const std::map<std::string, long double>& GetSTDmapLabelToCoefficient() const; //returns map (=representation of multivector)

	//**********************************OPERATORS**********************************\\
	
	GAQ RotorExponential(unsigned int degree, long double phi) const;			// Use carefully! Only works for specific elements, may crash otherwise
	GAQ TranslatorExponential(unsigned int degree, long double distance) const; // Use carefully! Only works for specific elements, may crash otherwise
	GAQ& operator=(const GAQ& other); 
	GAQ& operator=(GAQ&& other) noexcept; 
	bool operator==(const GAQ& other) const;	//equals operator
	bool operator!=(const GAQ& other) const;	//not equals operator
	GAQ operator[](int Grade) const;			//Grade projection
	GAQ operator[](const GAQ& other) const;	//basis blade selection
	GAQ operator*(const GAQ& other) const;	//geometric product operator
	GAQ operator*(GAQ&& other) const; 
	GAQ operator*(long double scalar) const;	//multiplying by scalar from the right operator
	GAQ operator~() const;						//reverse operator
	GAQ operator+(const GAQ& other) const;	//multivector addition operator
	GAQ operator-(const GAQ& other) const;	//multivector substraction operator
	GAQ operator|(const GAQ& other) const;	//inner product operator
	GAQ operator^(const GAQ& other) const;	//outer product operator
	GAQ operator^(int exponent) const;			//exponent operator
	GAQ operator/(long double divider) const;	//dividing operator
	GAQ ScalarProduct(const GAQ& b) const;

	static GAQ Rotate(const GAQ& point, rotation_planes plane, long double angle);
	static GAQ Translate(const GAQ& point, translation_directions plane, long double angle);

	int Grade(std::string_view label) const;	//returns Grade of basis blade (if we give it appropriate label...)
	std::string Log() const;					//returns multivector, used in << operator
	
	static GAQ generatingBlades[];				//stores 1,e1,e2,...,en.

protected:

	static int CalculateSign(const std::vector<int>& permutation); //Helps in validating basis Element, calculates sign of permutation
	static void SimplifyBasisBlade(std::string& label, int& sign); //simplifies label in a form of for example  e1e2e3e2e3 into e1
	static void processVector(std::vector<int>& vec, int& sign); //used when simplifying results of geometric product: e1e2e5e2e3e4e5 -> e1e5e3e4e5 -> e1e3e4 represented byjust numbers (1252345 -> 15345 ...)
	static std::vector<int> ExtractIntegersFromBasisBlades(std::string_view label); //from a given label, for example e1*e2*e3, returns vector {1,2,3}

	GAQ operator||(const GAQ& other) const; //inner product of two basis blades
	GAQ operator &&(const GAQ& other) const; //outer product of two basis blades
	GAQ operator ()(int Grade) const; //Grade projection of basis blade
	
	std::map<std::string, long double> m_mapLabelToCoefficient; //representation of a general multivector
	// 3 + 2e1 - e1*e2*e3
	// ==================
	// 1 |  e1 | e1*e2*e3 
	// 3 |   2 |       -1

private:
	void DeleteZeroFromVector(); //if multivector is of a form 0*1 + c1e1+ c2e1*e2 +... it removes 0*1
};
//**********************************NON-MEMBER_OPERATORS**********************************\\

GAQ operator*(long double scalar, const GAQ& onther); //multiplying by scalar from the left
std::ostream& operator<<(std::ostream& stream, const GAQ& vector); //operator for printing
std::vector<GAQ> MakeQCGAFromBasisBlades(const GAQ& multivector); //returns vector of basis blades in linear combination of general multivector
