#include "QCGA.h"
#include <cmath>
#include <numeric>

using namespace std::string_literals;

QCGA QCGA::generatingBlades[GENERATING_BASIS_SIZE + 1];

void QCGA::generateGeneratingBlades()
{
	QCGA::generatingBlades[0] = QCGA("1");
	for (char i = 1; i < GENERATING_BASIS_SIZE + 1; i++)
	{
		QCGA::generatingBlades[i] = QCGA("e" + std::to_string(i));
	}
}

QCGA::QCGA()
{
	m_mapLabelToCoefficient["1"] = 0;
}

QCGA::QCGA(const std::string& input)
{
	input == "0" ?
		  m_mapLabelToCoefficient["1"] = 0
		: m_mapLabelToCoefficient[input] = 1;
}

QCGA::QCGA(std::string&& input) noexcept
{
	input == "0" ?
		  m_mapLabelToCoefficient["1"] = 0
		: m_mapLabelToCoefficient[input] = 1;
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
	this->m_mapLabelToCoefficient = copyOfMap;
}

QCGA::QCGA(std::map<std::string, long double>&& map)
{
	this->m_mapLabelToCoefficient = std::move(map);

	auto ptr = &map;
	ptr = nullptr;
}

QCGA::QCGA(const std::pair<std::string, long double>& basis_blade)
{
	this->m_mapLabelToCoefficient.emplace(basis_blade);
}

QCGA::QCGA(std::pair<std::string, long double>&& basis_blade)
{
	this->m_mapLabelToCoefficient.emplace(std::move(basis_blade));
}


QCGA::QCGA(const QCGA& instance)
{
	this->m_mapLabelToCoefficient = instance.m_mapLabelToCoefficient;
}

QCGA::QCGA(QCGA&& instance) noexcept
{
	this->m_mapLabelToCoefficient = std::move(instance.m_mapLabelToCoefficient);
}


const std::map<std::string, long double>& QCGA::getSTDmapLabelToCoefficient() const
{
	return this->m_mapLabelToCoefficient;
}

long double QCGA::ToNumeric() // returns coefficient at basis blade 1
{
	if (this->m_mapLabelToCoefficient.find("1") != this->m_mapLabelToCoefficient.end())
	{
		return this->m_mapLabelToCoefficient.at("1");
	}
	std::cout << "WARNING! ToNumeric found no zero-degree basis blade, returned number 1";
	return 1.0;
}

//************************************MEMBER_OPERATOR************************************\\

QCGA QCGA::RotorExponential(unsigned int degree, long double phi) const
{
	QCGA res = one;
	for (unsigned int i = 1; i < degree + 1; i++)
	{
		long double factorial = i;
		for (int j = i; j > 1; j--)
		{
			factorial *= (j - 1);
		}
		res = res + (((phi / 2) * (*this)) ^ i) * (long double(1) / factorial);
	}
	return res;
}

QCGA QCGA::translatorExponential(unsigned int degree, long double distance) const
{
	QCGA res = one;
	for (unsigned int i = 1; i < degree + 1; i++)
	{
		long double factorial = i;
		for (int j = i; j > 1; j--)
		{
			factorial *= (j - 1);
		}
		long double factor = (pow(distance / 2, i)) * (long double(1) / factorial);
		res = res + (factor * ((*this) ^ i));
	}
	return res;
}

QCGA QCGA::BivectorExponential(unsigned int degree, long double parameter) const
{
	QCGA res = one;
	for (unsigned int i = 1; i < degree + 1; i++)
	{
		long double factorial = i;
		for (int j = i; j > 1; j--)
			factorial *= (j - 1);
		res = res + (((parameter / 2) * (*this)) ^ i) * (long double(1) / factorial);
	}
	return res;
}

QCGA& QCGA::operator=(const QCGA& other)
{
	this->m_mapLabelToCoefficient = other.m_mapLabelToCoefficient;
	return *this;
}

QCGA& QCGA::operator=(QCGA&& other) noexcept
{
	this->m_mapLabelToCoefficient = std::move(other.m_mapLabelToCoefficient);
	return *this;
}



//equal operator
bool QCGA::operator==(const QCGA& other) const
{
	bool equal = true;
	for (const auto& [basisBlade, value] : this->m_mapLabelToCoefficient)
	{
		if (other.m_mapLabelToCoefficient.find(basisBlade) == other.m_mapLabelToCoefficient.end()) //if there is on the right not the same basis blade as on the left
		{
			equal = false;
			break;
		}
		//if there is, check for coefs if they are the same
		if (abs(other.m_mapLabelToCoefficient.at(basisBlade) - value) > long double(10000) / (PRECISION)) // if coefs are different, 10000 is for not being so strict due to rounding errors
		{
			equal = false;
			break;
		}
	}
	for (const auto& [basisBlade, value] : other.m_mapLabelToCoefficient) //now vice versa
	{
		if (this->m_mapLabelToCoefficient.find(basisBlade) == this->m_mapLabelToCoefficient.end()) //if there is on the right not the same basis blade as on the left
		{
			equal = false;
			break;
		}
		//if there is, check for coefs if they are the same
		if (abs(other.m_mapLabelToCoefficient.at(basisBlade) - value) > long double(10000) / (PRECISION)) // if coefs are different, 10000 is for not being so strict due to rounding errors
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
QCGA QCGA::operator[](int _grade) const
{
	std::vector<QCGA> left = makeQCGAFromBasisBlades(*this);

	QCGA res = zero_vector;
	for (int i = 0; i < this->m_mapLabelToCoefficient.size(); i++)
	{
		res = res + (left[i])(_grade); //equivalent to standard formula 
	}
	res.deleteZeroFromVector();
	return res;
}

QCGA QCGA::operator[](const QCGA& other) const
{
	QCGA res;
	for (const auto& [key, coeff] : m_mapLabelToCoefficient) 
	{
		if (other.m_mapLabelToCoefficient.contains(key)) 
		{ 
			res = coeff * QCGA(key);
			break;
		}
	}
	return res;
}

//geometric product operator
QCGA QCGA::operator*(const QCGA& other) const
{
	std::map<std::string, long double> map;

	for (const auto& [thisBasisBlade, thisCoef] : this->m_mapLabelToCoefficient) //e1+e1*e2*e3
		for (const auto& [rightBasisBlade, rightCoef] : other.m_mapLabelToCoefficient) //2e2-e3*e4
			map[thisBasisBlade + "*" + rightBasisBlade] = thisCoef * rightCoef; // e1*e2 | e1*e3*e4 | e1*e2*e3*e2 | e1*e2*e3*e3*e4
																				//     2 |  	 -1 |			2 |				-1

	QCGA res = zero_vector;
	//now, try to simplify individual label
	for (const auto& [basisBlade, coef] : map) //each label in newLabel will be modified
	{ 
		//const long double oldCoef = coef;
		std::string copyOfBasisBlade = basisBlade;
		int sign = 1; //sign for controling sing when swaps happen

		if (!basisBlade.contains('e')) //if there are just greade 0 elements, label is 1
			copyOfBasisBlade = "1";
		
		else if ((basisBlade.find_last_of("*") + 1) != basisBlade.find_last_of("e")) //if there is ei*1, label will be ei
			copyOfBasisBlade = basisBlade.substr(0, basisBlade.find_last_of("*"));

		else if (basisBlade.find("e") > 0) //if there is 1*ei, label will be ei
			copyOfBasisBlade = basisBlade.substr(2, basisBlade.size());

		if (copyOfBasisBlade != "1")
			simplifyBasisBlade(copyOfBasisBlade, sign);//in case there is ei, for example e1e2e3e2e3 needs to be simplified in e1

		res = (res + std::move(QCGA(std::move(std::make_pair(copyOfBasisBlade, coef * sign)))));
	}
	res.deleteZeroFromVector();
	return res;
}

QCGA QCGA::operator*(QCGA&& other) const
{
	//these vectors will be used for constructor
	std::map<std::string, long double> map;

	for (const auto& [thisBasisBlade, thisCoef] : this->m_mapLabelToCoefficient)
		for (auto&& [rightBasisBlade, rightCoef] : other.m_mapLabelToCoefficient)
			map[std::move(thisBasisBlade + "*" + rightBasisBlade)] = std::move(thisCoef * rightCoef);

	QCGA res = zero_vector;
	//now, try to simplify individual label
	for (auto& [basisBlade, coef] : map) //each label in newLabel will be modified
	{
		const long double oldCoef = coef;
		std::string copyOfBasisBlade = basisBlade;
		int sign = 1; //sign for controling sing when swaps happen

		//if (basisBlade.find("e")==std::string::npos) //if there are just greade 0 elements, label is 1
		if (!basisBlade.contains('e')) //if there are just greade 0 elements, label is 1
			copyOfBasisBlade = "1";

		else if ((basisBlade.find_last_of("*") + 1) != basisBlade.find_last_of("e")) //if there is ei*1, label will be ei
			copyOfBasisBlade = basisBlade.substr(0, basisBlade.find_last_of("*"));

		else if (basisBlade.find("e") > 0) //if there is 1*ei, label will be ei
			copyOfBasisBlade = basisBlade.substr(2, basisBlade.size());

		if (copyOfBasisBlade != "1") //_LIKELY
			simplifyBasisBlade(copyOfBasisBlade, sign);//in case there is ei, for example e1e2e3e2e3 needs to be simplified in e1

		res = std::move((res + std::move(QCGA(std::move(std::make_pair(copyOfBasisBlade, oldCoef * sign))))));
	}
	res.deleteZeroFromVector();
	return res;
}

//multiplying by scalar from the right
QCGA QCGA::operator*(long double scalar) const
{
	std::map<std::string, long double> map = this->m_mapLabelToCoefficient;
	for (const auto& [basisBlade, coef] : m_mapLabelToCoefficient)
	{
		map.at(basisBlade) *= scalar;
	}
	return QCGA(std::move(map));
}

//reverse operator
QCGA QCGA::operator~() const
{
	std::map<std::string, long double> map;
	for (const auto& [basisBlade, coef] : m_mapLabelToCoefficient) //it is sufficient to check sign of reversed permutation of basis blades
	{
		if (basisBlade == "1")
			map["1"] = coef;
		else
		{
			std::vector<int> permutation = extractIntegersFromBasisBlades(basisBlade);
			std::reverse(permutation.begin(), permutation.end());
			map[basisBlade] = coef * calculateSign(permutation);
		}
	}
	return QCGA(std::move(map));
}

//addition operator
QCGA QCGA::operator+(const QCGA& other) const
{
	std::map<std::string, long double> map = this->m_mapLabelToCoefficient;
	for (const auto& [basisBlade, coefficient] : other.m_mapLabelToCoefficient)
	{
		//if labels are the same
		if (map.find(basisBlade) != map.end())
		{
			const long double value = map.at(basisBlade);
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
	std::vector<QCGA> left = makeQCGAFromBasisBlades(*this);
	std::vector<QCGA> right = makeQCGAFromBasisBlades(other);
	QCGA res = zero_vector;
	for (int i = 0; i < this->m_mapLabelToCoefficient.size(); i++)
	{
		for (int j = 0; j < other.m_mapLabelToCoefficient.size(); j++)
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
	std::vector<QCGA> left = makeQCGAFromBasisBlades(*this);
	std::vector<QCGA> right = makeQCGAFromBasisBlades(other);

	QCGA res = zero_vector;
	for (int i = 0; i < this->m_mapLabelToCoefficient.size(); i++)
	{
		for (int j = 0; j < other.m_mapLabelToCoefficient.size(); j++)
		{
			res = res + (left[i] && right[j]); //equivalent to standard formula 
		}
	}
	res.deleteZeroFromVector();
	return res;
}

QCGA QCGA::operator^(int exponent) const
{
	if (exponent < 0)
	{
		std::cout << "Warning, calling an inverse of QCGA, not of a blade, this might fail!\n";
		long double denominator = (((*this) * (~(*this))).ToNumeric());
		QCGA res = denominator * (~(*this));
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

QCGA QCGA::operator/(long double divider) const
{
	return *this * (1.0 / divider);
}

QCGA QCGA::scalarProduct(const QCGA& b) const
{
	return (*this * b)[0];
}


QCGA QCGA::rotate(const QCGA& point, rotation_planes plane, long double angle)
{
	QCGA rotor = zero_vector;
	switch (plane)
	{
	case xy:
		rotor = rxy.RotorExponential(20, angle);
		break;
	case xz:
		rotor = rxz.RotorExponential(20, angle);
		break;
	case yz:
		rotor = ryz.RotorExponential(20, angle);
		break;
	default:
		std::cout << "Wrong rotation plane" << std::endl;
		return point;
	}
	return (rotor * point * ~rotor)[1];
}

QCGA QCGA::translate(const QCGA& point, translation_directions direction, long double distance)
{
	switch (direction)
	{
	case x:
		return Tx * point * ~Tx;
	case y:
		return Ty * point * ~Ty;
	case z:
		return Tz * point * ~Tz;
	default:
		std::cout << "Wrong direction for translating" << std::endl;
		return point;
	}
}

//returns grade of basis blade (if we give it appropriate label...)
int QCGA::grade(std::string_view label) const
{
	int grade = 0;
	for (char c : label) 
	{
		if (c == 'e') 
		{
			grade++;
		}
	}
	return grade;
}

//returns multivector in e_inf, e_o basis, used in << operator
std::string QCGA::log() const
{
	std::string s;
	if (this->m_mapLabelToCoefficient.find("1") != this->m_mapLabelToCoefficient.end() && this->m_mapLabelToCoefficient.at("1") == 0)
	{
		s = "0";
	}
	else
	{
		std::string coef;
		for (auto& a : m_mapLabelToCoefficient)
		{
			if (abs(a.second) < (double(1000000)/PRECISION))
			{
				continue;
			}
			coef = std::to_string(a.second);
			size_t j{ coef.length() - 1 };
			while (coef[j] == '0' && coef[j] != '.')
			{
				coef.erase(coef.end() - 1);
				j--;
			}
			if (coef[j] == '.')
			{
				coef.erase(coef.end() - 1);
			}
			s += coef + "*" + std::string(a.first) + " + ";
		}
	}
	return s.substr(0, s.size() - 3);
}

//************************************PRROTECTED************************************\\

//calculates sign of permutation
int QCGA::calculateSign(const std::vector<int>& permutation) 
{
	int sign = 1;
	for (size_t i = 0; i < permutation.size(); i++) 
	{
		for (size_t j = i + 1; j < permutation.size(); j++) 
		{
			if (permutation[i] > permutation[j]) 
			{
				sign = -sign; // Switch the sign for each inversion
			}
		}
	}
	return sign;
}

//simplifies label in a form of for example  e1e2e3e2e3 into e1
void QCGA::simplifyBasisBlade(std::string& label, int& sign)
{
	std::vector<int> permutations; //keeps numbers next to individual e's
	permutations.reserve(30);
	std::string s;
	while (label.contains('*'))
	{
		permutations.emplace_back(std::stoi(label.substr(label.find("e") + 1, label.find("*") - 1)));
		label = label.substr(label.find("*") + 1, label.size());
	}
	permutations.emplace_back(std::stoi(label.substr(label.find("e") + 1, label.size())));
	processVector(permutations, sign); //e1e2e5e2e3e4e5 -> e1e5e3e4e5 -> e1e3e4 represented by numbers (1252345 -> 15345 ...)
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
	{
		s = s.substr(0, s.length() - 1);
	}
	label = s;
}

//used when simplifying results of geometric product: e1e2e5e2e3e4e5 -> e1e5e3e4e5 -> e1e3e4 represented byjust numbers (1252345 -> 15345 ...)
void QCGA::processVector(std::vector<int>& vec, int& sign) 
{
	for (int i = 0; i < vec.size() - 1; i++)
	{
		for (int j = i + 1; j < vec.size(); j++) 
		{
			if (vec[i] != vec[j])
			{
				continue;
			}
			if (vec[i] > ALGEBRA_P)
			{
				sign *= -1;
			}
			for (int k = 0; k < j - i - 1; k++)
			{
				std::swap(vec[j - k - 1], vec[j - k]);
				sign *= -1;
			}
			vec.erase(vec.begin() + i, vec.begin() + i + 2);
			if (vec.empty())
			{
				break;
			}
			processVector(vec, sign);
		}
		if (vec.empty())
		{
			break;
		}
	}
	if (vec.empty())
	{
		return;
	}
	//bubble sort, easy to track swaps... e3e2e1 -> e1e2e3
	int i, j;
	for (i = 0; i < vec.size() - 1; i++)
	{
		for (j = 0; j < vec.size() - i - 1; j++)
		{
			if (vec[j] > vec[j + 1])
			{
				std::swap(vec[j], vec[j + 1]);
				sign *= -1;
			}
		}
	}
}

//from a given label, for example e1*e2*e3, returns vector {1,2,3}
std::vector<int> QCGA::extractIntegersFromBasisBlades(std::string_view label)
{
	std::vector<int> permutation;
	permutation.reserve(15);
	int basis_vec_number;
	auto position = label.begin();
	while (position < label.end()) //checks, if idividual ei are presented in algebra. Otherwise, it crashes
	{
		auto [ptr, error] {std::from_chars(position._Unwrapped(), position._Unwrapped() + 2, basis_vec_number)};
		if (error == std::errc{}) 
		{
			permutation.push_back(basis_vec_number);
			if (basis_vec_number > 9) 
			{ 
				position++; 
			}
		}
		position++;
	}
	return permutation;
}

//**********************************STATIC_SUPPORT_OPERATORS**********************************\\

//inner product of two basis blades operator. Usefull in inner product of two general multivectors
QCGA QCGA::operator||(const QCGA& other) const
{
	return QCGA((*this * other)[abs(this->grade(this->m_mapLabelToCoefficient.begin()->first) - other.grade(other.m_mapLabelToCoefficient.begin()->first))]);
}

//outer product of two basis blades operator
QCGA QCGA::operator&&(const QCGA& other) const
{
	return QCGA((*this * other)[abs(this->grade(this->m_mapLabelToCoefficient.begin()->first) + other.grade(other.m_mapLabelToCoefficient.begin()->first))]);
}

//grade projection of basis blade
QCGA QCGA::operator()(int grade) const
{
	return (this->grade(this->m_mapLabelToCoefficient.begin()->first) == grade) ? *this : QCGA("0");
}

void QCGA::deleteZeroFromVector()
{
	if (this->m_mapLabelToCoefficient.find("1") != this->m_mapLabelToCoefficient.end())
	{
		if (this->m_mapLabelToCoefficient.at("1") == 0 && this->m_mapLabelToCoefficient.size() > 1)
		{
			m_mapLabelToCoefficient.erase("1");
		}
	}
}


//**********************************NON-MEMBER_OPERATORS**********************************\\

//multiplying by scalar from the left
QCGA operator*(const long double& scalar, const QCGA& other)
{
	return other * scalar;
}

// Operator for printing
std::ostream& operator<<(std::ostream& stream, const QCGA& multivector)
{
	stream << multivector.log();
	return stream;
}

//returns vector of basis blades in linear combination of general multivector
std::vector<QCGA> makeQCGAFromBasisBlades(const QCGA& multivector)
{
	std::vector<QCGA> basisBlades;
	basisBlades.reserve(20);
	for (const auto& [basisBlade, coef] : multivector.getSTDmapLabelToCoefficient())
	{
		basisBlades.emplace_back(std::make_pair(basisBlade, coef));
	}
	return basisBlades;
}