#include "QCGA.h"
#include <cmath>
#include <numeric>

QCGA QCGA::generatingBlades[GENERATING_BASIS_DIMENSION + 1];
//generates generatingBlades
void QCGA::generateGeneratingBlades()
{
	QCGA::generatingBlades[0] = QCGA("1");
	for (int i = 1; i < GENERATING_BASIS_DIMENSION + 1; i++)
	{
		QCGA::generatingBlades[i] = QCGA("e" + std::to_string(i));
	}
}

//Constructor which takes basis blade
QCGA::QCGA(const std::string& input)
{
	if (input == "0")
	{
		STDmapLabelToCoefficient["1"] = 0;
	}
	else
	{
		STDmapLabelToCoefficient[input] = std::round(1 * PRECISION) / PRECISION;
	}
}

QCGA::QCGA(const std::map<std::string, long double>& map)
{
	std::map<std::string, long double> copyOfMap = map;
	for (const auto& [basisBlade, coef] : map)
	{
		//copyOfMap.at(basisBlade) = std::round(coef * PRECISION) / PRECISION; //tried to fix rounding error in another way
		if (abs(coef) <= long double(1) / PRECISION)
		{
			copyOfMap.erase(basisBlade);
		}
	}
	if (copyOfMap.empty())
	{
		copyOfMap["1"] = 0;
	}
	this->STDmapLabelToCoefficient = copyOfMap;
}

QCGA::QCGA(const std::pair<std::string, long double>& basis_blade)
{
	this->STDmapLabelToCoefficient.emplace(basis_blade);
}

//default constructor calls CGA("0")
QCGA::QCGA() { STDmapLabelToCoefficient["1"] = 0; }

const std::map<std::string, long double>& QCGA::getSTDmapLabelToCoefficient() const
{
	return this->STDmapLabelToCoefficient;
}

long double QCGA::toNumeric() // returns coefficient at basis blade 1
{
	if (this->STDmapLabelToCoefficient.find("1") != this->STDmapLabelToCoefficient.end())
	{
		return this->STDmapLabelToCoefficient.at("1");
	}
	else
	{
		std::cout << "WARNING! toNumeric found no zero-degree basis blade, returned number 1";
		return 1;
	}
}

//************************************MEMBER_OPERATOR************************************\\

QCGA QCGA::rotorExponential(unsigned int degree, long double phi)
{
	QCGA res = one;
	for (int i = 1; i < degree+1; i++)
	{
		long double factorial = i;
		for (int j = i; j > 1; j--)
			factorial *= (j - 1);
		res = res + (((phi / 2) * (*this)) ^ i)* (long double(1) / factorial);
	}
	return res;
}

QCGA QCGA::translatorExponential(unsigned int degree, long double distance)
{
	QCGA res = one;
	for (int i = 1; i < degree + 1; i++)
	{
		long long unsigned int factorial = i;
		for (int j = i; j > 1; j--)
			factorial *= (j - 1);
		long double factor = (pow(distance/2, i)) * (long double(1) / factorial);
		res = res + (factor * ((*this) ^ i));
	}
	return res;
}


//equal operator
bool QCGA::operator==(const QCGA& other) const
{
	bool equal = true;
	for (const auto& [basisBlade, value] : this->STDmapLabelToCoefficient)
	{
		if (other.STDmapLabelToCoefficient.find(basisBlade) == other.STDmapLabelToCoefficient.end()) //if there is on the right not the same basis blade as on the left
		{
			equal = false;
			break;
		}
		//if there is, check for coefs if they are the same
		if (abs(other.STDmapLabelToCoefficient.at(basisBlade) - value) > long double(10000) / (PRECISION)) // if coefs are different, 10000 is for not being so strict due to rounding errors
		{
			equal = false;
			break;
		}
	}
	for (const auto& [basisBlade, value] : other.STDmapLabelToCoefficient) //now vice versa
	{
		if (this->STDmapLabelToCoefficient.find(basisBlade) == this->STDmapLabelToCoefficient.end()) //if there is on the right not the same basis blade as on the left
		{
			equal = false;
			break;
		}
		//if there is, check for coefs if they are the same
		if (abs(other.STDmapLabelToCoefficient.at(basisBlade) - value) > long double(10000) / (PRECISION)) // if coefs are different, 10000 is for not being so strict due to rounding errors
		{
			equal = false;
			break;
		}

	}
	return equal;
}

//not equal operator
bool QCGA::operator!=(const QCGA& other) const
{
	return !(*this == other);
}

//grade projection operator
QCGA QCGA::operator[](const int _grade) const
{
	std::vector<QCGA> left = makeCGAFromBasisBlades(*this);

	QCGA res = zero_vector;
	for (int i = 0; i < this->STDmapLabelToCoefficient.size(); i++)
	{
		res = res + (left[i])(_grade); //equivalent to standard formula 
	}
	res.deleteZeroFromVector();
	return res;
}

//geometric product operator
QCGA QCGA::operator*(const QCGA& other) const
{
	//these vectors will be used for constructor
	std::map<std::string, long double> map;

	for (const auto& [thisBasisBlade, thisCoef] : this->STDmapLabelToCoefficient)
		for (const auto& [rightBasisBlade, rightCoef] : other.STDmapLabelToCoefficient)
			map[thisBasisBlade + "*" + rightBasisBlade] = thisCoef * rightCoef;

	QCGA res = zero_vector;
	//now, try to simplify individual label
	for (auto& [basisBlade, coef] : map) //each label in newLabel will be modified
	{ 
		long double oldCoef = coef;
		std::string copyOfBasisBlade = basisBlade;
		int sign = 1; //sign for controling sing when swaps happen

		if (basisBlade.find("e") == std::string::npos)//if there are just greade 0 elements, label is 1
			copyOfBasisBlade = "1";
		
		else if ((basisBlade.find_last_of("*") + 1) != basisBlade.find_last_of("e")) //if there is ei*1, label will be ei
			copyOfBasisBlade = basisBlade.substr(0, basisBlade.find_last_of("*"));

		else if (basisBlade.find("e") > 0) //if there is 1*ei, label will be ei
			copyOfBasisBlade = basisBlade.substr(2, basisBlade.size());

		if (copyOfBasisBlade != "1")
			simplifyBasisBlade(copyOfBasisBlade, sign);//in case there is ei, for example e1e2e3e2e3 needs to be simplified in e1

		res = res + QCGA(std::make_pair(copyOfBasisBlade, oldCoef * sign));
	}
	res.deleteZeroFromVector();
	return res;
}

//multiplying by scalar from the right
QCGA QCGA::operator*(const long double scalar) const
{
	//std::map<std::string, long double> map;
	std::map<std::string, long double> map = this->STDmapLabelToCoefficient;

	for (const auto& [basisBlade, coef] : STDmapLabelToCoefficient)
	{
		//map.emplace(std::make_pair(basisBlade, coef * scalar));// *= scalar;
		map.at(basisBlade) *= scalar;
	}
	return QCGA(map);
}

//reverse operator
QCGA QCGA::operator~() const
{
	std::map<std::string, long double> map;

	for (const auto& [basisBlade, coef] : STDmapLabelToCoefficient) //it is sufficient to check sign of reversed permutation of basis blades
	{
		if (basisBlade == "1")
			map["1"] = coef;//coefficients.push_back(this->coefficients[0]);
		else
		{
			std::vector<int> permutation = extractIntegersFromBasisBlades(basisBlade);
			std::reverse(permutation.begin(), permutation.end());
			map[basisBlade] = coef * calculateSign(permutation);//coefficients.push_back(this->coefficients[iter] * calculateSign(permutation));
		}
	}
	return QCGA(map);
}

//addition operator
QCGA QCGA::operator+(const QCGA& other) const
{
	std::map<std::string, long double> map = this->STDmapLabelToCoefficient;
	for (const auto& [basisBlade, coefficient] : other.STDmapLabelToCoefficient)
	{
		//if labels are the same
		if (map.find(basisBlade) != map.end())
		{
			long double value = map.at(basisBlade);
			if (value + coefficient == 0) // if e1 - e1 for example.
			{
				map.erase(basisBlade);
				if (map.empty()) //if map is null vector
				{
					map["1"] = 0;
				}
			}
			else
			{
				map.at(basisBlade) += coefficient;
			}
		}
		else
		{
			map[basisBlade] = coefficient;
		}
	}
	return QCGA(map);
}

//substraction operator
QCGA QCGA::operator-(const QCGA& other) const
{
	return *this + ((-1) * other);
}

//inner product operator
QCGA QCGA::operator|(const QCGA& other) const
{
	//Make vectors from left and right operadns, they will store basis blades as CGA object and they will participate in product
	std::vector<QCGA> left = makeCGAFromBasisBlades(*this);
	std::vector<QCGA> right = makeCGAFromBasisBlades(other);

	QCGA res = zero_vector;
	for (int i = 0; i < this->STDmapLabelToCoefficient.size(); i++)
	{
		for (int j = 0; j < other.STDmapLabelToCoefficient.size(); j++)
		{
			res = res + (left[i] || right[j]); //equivalent to standard formula 
		}
	}
	res.deleteZeroFromVector();
	return res;
}

//outer product operator
QCGA QCGA::operator^(const QCGA& other) const
{
	//Make vectors from left and right operadns, they will store basis blades as CGA object and they will participate in product
	std::vector<QCGA> left = makeCGAFromBasisBlades(*this);
	std::vector<QCGA> right = makeCGAFromBasisBlades(other);

	QCGA res = zero_vector;
	for (int i = 0; i < this->STDmapLabelToCoefficient.size(); i++)
	{
		for (int j = 0; j < other.STDmapLabelToCoefficient.size(); j++)
		{
			res = res + (left[i] && right[j]); //equivalent to standard formula 
		}
	}
	res.deleteZeroFromVector();
	return res;
}

QCGA QCGA::operator^(const int exponent) const
{
	if (exponent < 0)
	{
		std::cout << "Warning, calling an inverse of QCGA, not of a blade, this might fail!\n";
		QCGA res = (~*this) / ((*this * ~(*this)).toNumeric());
		return res;
	}
	else
	{
		QCGA res = *this;
		for (int i = 0; i < exponent - 1; i++)
		{
			res = res * *this;
		}
		return res;
	}
}

QCGA QCGA::operator/(const long double divider) const
{
	return *this * (1 / divider);
}

QCGA QCGA::scalarProduct(const QCGA& b)
{
	return (*this * b)[0];
}


QCGA QCGA::rotate(const QCGA& point, int plane, long double angle)
{
	QCGA rotor = zero_vector;
	switch (plane)
	{
	case 1:
		rotor = rxy.rotorExponential(20, angle);
		break;
	case 2:
		rotor = rxz.rotorExponential(20, angle);
		break;
	case 3:
		rotor = ryz.rotorExponential(20, angle);
		break;
	default:
		std::cout << "Wrong rotation plane" << std::endl;
		return point;
		break;
	}
	return (rotor * point * ~rotor)[1];
}

QCGA QCGA::translate(const QCGA& point, int direction, long double distance)
{
	switch (direction)
	{
	case 1:
		return Tx * point * ~Tx;
		break;
	case 2:
		return Ty * point * ~Ty;
		break;
	case 3:
		return Tz * point * ~Tz;
		break;
	default:
		std::cout << "Wrong direction for translating" << std::endl;
		return point;
		break;
	}
}

//returns grade of basis blade (if we give it appropriate label...)
int QCGA::grade(const std::string& label) const
{
	int grade = 0;
	for (const char c : label) {
		if (c == 'e') {
			grade++;
		}
	}
	return grade;
}

//returns multivector in e_inf, e_o basis, used in << operator
std::string QCGA::log() const
{
	std::string s = "";
	if (this->STDmapLabelToCoefficient.find("1") != this->STDmapLabelToCoefficient.end() && this->STDmapLabelToCoefficient.at("1") == 0)
	{
		s = "0";
	}
	else
	{
		std::string coef;
		for (auto& a : STDmapLabelToCoefficient)
		{
			if (abs(a.second) < (double(1000000)/PRECISION))
			{
				continue;
			}
			coef = std::to_string(a.second);
			int j = coef.length() - 1;
			while (coef[j] == '0' && coef[j] != '.')
			{
				coef.erase(coef.end() - 1);
				j--;
			}
			if (coef[j] == '.')
				coef.erase(coef.end() - 1);
			s += coef + "*" + a.first + " + ";
		}
	}
	s = s.substr(0, s.size() - 3); //deletes " + " at the end
	return s;
}

//************************************PRROTECTED************************************\\

//calculates sign of permutation
int QCGA::calculateSign(const std::vector<int>& permutation) {
	int sign = 1; // Initialize the sign to positive

	for (size_t i = 0; i < permutation.size(); i++) {
		for (size_t j = i + 1; j < permutation.size(); j++) {
			if (permutation[i] > permutation[j]) {
				sign = -sign; // Switch the sign for each inversion
			}
		}
	}
	return sign;
}

//simplifies label in a form of for example  e1e2e3e2e3 into e1
void QCGA::simplifyBasisBlade(std::string& label, int& sign)
{
	//std::cout << label << std::endl;
	std::vector<int> permutations; //keeps numbers next to individual e's
	permutations.reserve(30);
	std::string s = "";
	while (label.find("*") != std::string::npos)
	{
		permutations.emplace_back(std::stoi(label.substr(label.find("e") + 1, label.find("*") - 1)));
		label = label.substr(label.find("*") + 1, label.size());
	}
	permutations.emplace_back(std::stoi(label.substr(label.find("e") + 1, label.size())));
	processVector(permutations, sign); //e1e2e5e2e3e4e5 -> e1e5e3e4e5 -> e1e3e4 represented byjust numbers (1252345 -> 15345 ...)
	for (int i = 0; i < permutations.size(); i++) //creates new proper label
	{
		s += "e";
		s += std::to_string(permutations[i]);
		s += "*";
	}
	if (s.size() == 0)
	{
		s += "1";
	}
	else
		s = s.substr(0, s.length() - 1);
	label = s;
}

//used when simplifying results of geometric product: e1e2e5e2e3e4e5 -> e1e5e3e4e5 -> e1e3e4 represented byjust numbers (1252345 -> 15345 ...)
void QCGA::processVector(std::vector<int>& vec, int& sign) {
	for (int i = 0; i < vec.size() - 1; i++)
	{
		for (int j = i + 1; j < vec.size(); j++) {
			if (vec[i] == vec[j])
			{
				if (vec[i] > ALGEBRA_P)
					sign *= -1;
				for (int k = 0; k < j - i - 1; k++)
				{
					std::swap(vec[j - k - 1], vec[j - k]);
					sign *= -1;
				}
				vec.erase(vec.begin() + i, vec.begin() + i + 2);
				if (vec.size() == 0)
					break;
				processVector(vec, sign);
			}
		}
		if (vec.size() == 0)
			break;
	}
	//bubble sort, easy to track swaps... e3e2e1 -> e1e2e3
	if (vec.size() != 0)
	{
		int i, j;
		for (i = 0; i < vec.size() - 1; i++)
			for (j = 0; j < vec.size() - i - 1; j++)
				if (vec[j] > vec[j + 1])
				{
					std::swap(vec[j], vec[j + 1]);
					sign *= -1;
				}
	}
}

//from a given label, for example e1*e2*e3, returns vector {1,2,3}
std::vector<int> QCGA::extractIntegersFromBasisBlades(const std::string& label)
{
	std::vector<int> permutation;
	std::string vec = label;
	while (vec.size() != 0) //checks, if idividual ei are presented in algebra. Otherwise, it crashes
	{
		std::string test;
		if (vec.substr(1, vec.length()).find("e") == std::string::npos)
		{
			test = vec.substr(1, vec.length());
			vec = "";
		}
		else
		{
			size_t positionOfNext = vec.substr(1, vec.length()).find("e") + 1;
			test = vec.substr(1, positionOfNext - 2);
			vec = vec.substr(positionOfNext, vec.length());
		}
		permutation.emplace_back(std::stoi(test));
	}
	return permutation;
}

//**********************************STATIC_SUPPORT_OPERATORS**********************************\\

//inner product of two basis blades operator. Usefull in inner product of two general multivectors
QCGA QCGA::operator||(const QCGA& other) const
{
	return QCGA((*this * other)[abs(this->grade(this->STDmapLabelToCoefficient.begin()->first) - other.grade(other.STDmapLabelToCoefficient.begin()->first))]);
}

//outer product of two basis blades operator
QCGA QCGA::operator&&(const QCGA& other) const
{
	return QCGA((*this * other)[abs(this->grade(this->STDmapLabelToCoefficient.begin()->first) + other.grade(other.STDmapLabelToCoefficient.begin()->first))]);
}

//grade projection of basis blade
QCGA QCGA::operator()(const int grade) const
{
	return (this->grade(this->STDmapLabelToCoefficient.begin()->first) == grade) ? *this : QCGA("0");
}

void QCGA::deleteZeroFromVector()
{
	if (this->STDmapLabelToCoefficient.find("1") != this->STDmapLabelToCoefficient.end())
	{
		if (this->STDmapLabelToCoefficient.at("1") == 0 && this->STDmapLabelToCoefficient.size() > 1)
		{
			STDmapLabelToCoefficient.erase("1");
		}
	}
}


//**********************************NON-MEMBER_OPERATORS**********************************\\

//multiplying by scalar from the left
QCGA operator*(const long double scalar, const QCGA& other)
{
	return other * scalar;
}

// Operator for printing
std::ostream& operator<<(std::ostream& stream, const QCGA& vector)
{
	stream << vector.log();
	return stream;
}

//returns vector of basis blades in linear combination of general multivector
std::vector<QCGA> makeCGAFromBasisBlades(const QCGA& multivector)
{
	std::vector<QCGA> basisBlades;
	basisBlades.reserve(20);
	for (const auto& [basisBlade, coef] : multivector.getSTDmapLabelToCoefficient())
	{
		basisBlades.emplace_back(std::make_pair(basisBlade, coef));
	}
	return basisBlades;
}

//removes occurences of substring in string
void removeOccurences(std::string& str, const std::string subStr)
{
	while (str.find(subStr) != std::string::npos)
	{
		size_t pos = str.find(subStr);
		if (pos != std::string::npos)
		{
			// Remove found substring from string
			str.erase(pos, subStr.length());
		}
	}

}