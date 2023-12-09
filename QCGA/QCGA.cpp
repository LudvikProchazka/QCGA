#include "QCGA.h"

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
QCGA::QCGA(std::string input)
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

QCGA::QCGA(std::map<std::string, double> map)
{
	std::map<std::string, double> copyOfMap = map;
	for (const auto& [basisBlade, coef] : map)
	{
		if (abs(coef) <= double(1) / PRECISION)
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

//instanciate CGA object from given objects
QCGA::QCGA(const QCGA& Multivector) { *this = Multivector; }

//default constructor calls CGA("0")
QCGA::QCGA() { *this = QCGA("0"); }

std::map<std::string, double> QCGA::getSTDmapLabelToCoefficient() const
{
	return this->STDmapLabelToCoefficient;
}

double QCGA::toNumeric() // returns coefficient at basis blade 1
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

//equal operator
bool QCGA::operator==(const QCGA& other) const
{
	bool equal = true;
	for (const auto& [basisBlade, value] : this->STDmapLabelToCoefficient)
	{
		auto otherMap = other.getSTDmapLabelToCoefficient();
		if (otherMap.find(basisBlade) == otherMap.end()) //if there is on the right not the same basis blade as in the left
		{
			equal = false;
			break;
		}
		else //if there is, check for coefs if they are the same
		{
			if (otherMap.at(basisBlade) != value) // if coefs are different
			{
				equal = false;
				break;
			}
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

	QCGA res("0");
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
	std::map<std::string, double> map;

	for (const auto& [thisBasisBlade, thisCoef] : this->STDmapLabelToCoefficient)
		for (const auto& [rightBasisBlade, rightCoef] : other.getSTDmapLabelToCoefficient())
			map[thisBasisBlade + "*" + rightBasisBlade] = thisCoef * rightCoef;

	std::map<std::string, double> singleMap;
	std::vector<QCGA*> basisBlades;
	basisBlades.reserve(map.size());
	QCGA res("0");
	//now, try to simplify individual label
	int iter = 0;
	std::map<std::string, double>mapCopy = map;
	for (auto& [basisBlade, coef] : map) //each label in newLabel will be modified
	{ //Ted je tady problem, ze to map primo menim. ez fix: meni jinou map do ktere sypat vysledky
		double oldCoef = coef;
		std::string copyOfBasisBlade = basisBlade;
		int sign = 1; //sign for controling sing when swaps happen

		if (basisBlade.find("e") == std::string::npos)//if there are just greade 0 elements, label is 1
		{
			mapCopy.erase(basisBlade);
			mapCopy["1"] = oldCoef;
		}
		else if ((basisBlade.find_last_of("*") + 1) != basisBlade.find_last_of("e")) //if there is ei*1, label will be ei
		{
			mapCopy.erase(basisBlade);
			copyOfBasisBlade = basisBlade.substr(0, basisBlade.find_last_of("*"));
			mapCopy[copyOfBasisBlade] = oldCoef;

		}
		else if (basisBlade.find("e") > 0) //if there is 1*ei, label will be ei
		{

			mapCopy.erase(basisBlade);
			copyOfBasisBlade = basisBlade.substr(2, basisBlade.size());
			mapCopy[copyOfBasisBlade] = oldCoef;
		}

		if (copyOfBasisBlade != "1")
		{
			simplifyBasisBlade(copyOfBasisBlade, sign);//in case there is ei, for example e1e2e3e2e3 needs to be simplified in e1
			mapCopy.erase(basisBlade);
			mapCopy[copyOfBasisBlade] = oldCoef * sign;
		}
		//coefficients[iter] *= sign; //keep track of sign base on swaps in siplifying process

		singleMap[copyOfBasisBlade] = oldCoef * sign;
		basisBlades.emplace_back(new QCGA(singleMap));
		res = res + *basisBlades[0];
		for (QCGA* basisBlade : basisBlades)
		{
			delete basisBlade;
		}
		basisBlades.clear();
		singleMap.clear();
		iter++;
	}
	res.deleteZeroFromVector();
	return res;
}

//multiplying by scalar from the right
QCGA QCGA::operator*(const double scalar) const
{
	//std::vector<std::string> Label = this->Label;
	//std::vector<double> coefficients = this->coefficients;
	std::map<std::string, double> map = this->STDmapLabelToCoefficient;

	for (const auto& [basisBlade, coef] : STDmapLabelToCoefficient)
	{
		map.at(basisBlade) *= scalar;
	}
	return QCGA(map);
}

//reverse operator
QCGA QCGA::operator~() const
{
	//std::vector<double> coefficients; //these 2 vectors will be used for private constructor
	//std::vector<std::string> newLabel = this->getLabel();
	std::map<std::string, double> map;

	int iter = 0;
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
		iter++;
	}
	return QCGA(map);
}

//addition operator
QCGA QCGA::operator+(const QCGA& other) const
{
	std::map<std::string, double> map = this->STDmapLabelToCoefficient;
	for (const auto& [basisBlade, coefficient] : other.getSTDmapLabelToCoefficient())
	{
		//if labels are the same
		if (map.find(basisBlade) != map.end())
		{
			double value = map.at(basisBlade);
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

	QCGA res("0");
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

	QCGA res("0");
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
	QCGA res = *this;
	for (int i = 0; i < exponent - 1; i++)
	{
		res = res * *this;
	}
	return res;
}

QCGA QCGA::operator/(const double divider) const
{
	return *this * (1 / divider);
}

//returns grade of basis blade (if we give it appropriate label...)
int QCGA::grade(std::string label) const
{
	unsigned int grade = 0;
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
			coef = std::to_string(a.second);
			for (int j = coef.length() - 1; j >= 0; j--) //deletes zeros behing decimal point if possible
			{
				if (coef[j] == '0' || coef[j] == '.')
				{
					coef.erase(coef.end() - 1);
				}
				else break;
			}
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

//returns bool depending if target string is present in string array
bool QCGA::searchString(std::string string[], int size, std::string target)
{
	bool b = false;
	for (int j = 0; j < size; j++)
		if (string[j] == target)
		{
			b = true;
			break;
		}
	return b;
}

//simplifies label in a form of for example  e1e2e3e2e3 into e1
void QCGA::simplifyBasisBlade(std::string& label, int& sign)
{
	std::vector<int> permutations; //keeps numbers next to individual e's
	std::string s = "";
	while (label.find("*") != std::string::npos)
	{
		s = label.substr(label.find("e") + 1, label.find("*") - 1);
		label = label.substr(label.find("*") + 1, label.size());
		permutations.push_back(std::stoi(s));
	}
	s = label.substr(label.find("e") + 1, label.size());
	permutations.push_back(std::stoi(s));
	processVector(permutations, sign); //e1e2e5e2e3e4e5 -> e1e5e3e4e5 -> e1e3e4 represented byjust numbers (1252345 -> 15345 ...)
	s = ""; //used for new fixed label
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
std::vector<int> QCGA::extractIntegersFromBasisBlades(std::string label)
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
		permutation.push_back(std::stoi(test));
	}
	return permutation;
}

//**********************************STATIC_SUPPORT_OPERATORS**********************************\\

//inner product of two basis blades operator. Usefull in inner product of two general multivectors
QCGA QCGA::operator||(const QCGA& other) const
{
	return QCGA((*this * other)[abs(this->grade(this->STDmapLabelToCoefficient.begin()->first) - other.grade(other.getSTDmapLabelToCoefficient().begin()->first))]);
}

//outer product of two basis blades operator
QCGA QCGA::operator&&(const QCGA& other) const
{
	return QCGA((*this * other)[abs(this->grade(this->STDmapLabelToCoefficient.begin()->first) + other.grade(other.getSTDmapLabelToCoefficient().begin()->first))]);
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
QCGA operator*(const double scalar, const QCGA& other)
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
	std::map<std::string, double> map;

	std::vector<QCGA> basisBlades;
	int iter = 0;
	for (const auto& [basisBlade, coef] : multivector.getSTDmapLabelToCoefficient())
	{
		map[basisBlade] = coef;
		basisBlades.emplace_back(QCGA(map));
		map.clear();
		iter++;
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