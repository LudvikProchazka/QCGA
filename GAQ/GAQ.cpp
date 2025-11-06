#include "GAQ.h"
#include <cmath>
#include <numeric>

using namespace std::string_literals;

GAQ GAQ::generatingBlades[GENERATING_BASIS_SIZE + 1];

template <typename T>
constexpr inline T abs_diff(T a, T b) {	return (a > b) ? (a - b) : (b - a); }
template <typename T>
constexpr inline T abs_sum(T a, T b) 
{
	if constexpr (std::is_signed_v<T>) 
	{
		return (a > 0 ? a : -a) + (b > 0 ? b : -b);
	}
	return a + b;
}

void GAQ::GenerateGeneratingBlades()
{
	GAQ::generatingBlades[0] = GAQ("1");
	for (char i = 1; i < GENERATING_BASIS_SIZE + 1; i++)
	{
		GAQ::generatingBlades[i] = GAQ("e" + std::to_string(i));
	}
}

GAQ::GAQ()
{
	m_mapLabelToCoefficient["1"] = 0;
}

GAQ::GAQ(const std::string& input)
{
	input == "0" ?
		  m_mapLabelToCoefficient["1"] = 0
		: m_mapLabelToCoefficient[input] = 1;
}

GAQ::GAQ(std::string&& input) noexcept
{
	input == "0" ?
		  m_mapLabelToCoefficient["1"] = 0
		: m_mapLabelToCoefficient[input] = 1;
}

GAQ::GAQ(const std::map<std::string, long double>& map)
{
	m_mapLabelToCoefficient = map;
	for (const auto& [basisBlade, coef] : map)
	{
		if (abs(coef) <= long double(1) / PRECISION)
		{
			m_mapLabelToCoefficient.erase(basisBlade);
		}
	}
	if (m_mapLabelToCoefficient.empty())
	{
		m_mapLabelToCoefficient["1"] = 0;
	}
}

GAQ::GAQ(std::map<std::string, long double>&& map)
{
	m_mapLabelToCoefficient = std::move(map);
}

GAQ::GAQ(const std::pair<std::string, long double>& basis_blade)
{
	m_mapLabelToCoefficient.emplace(basis_blade);
}

GAQ::GAQ(std::pair<std::string, long double>&& basis_blade)
{
	m_mapLabelToCoefficient.emplace(std::move(basis_blade));
}

GAQ::GAQ(const GAQ& instance)
{
	m_mapLabelToCoefficient = instance.m_mapLabelToCoefficient;
}

GAQ::GAQ(GAQ&& instance) noexcept
{
	m_mapLabelToCoefficient = std::move(instance.m_mapLabelToCoefficient);
}

const std::map<std::string, long double>& GAQ::GetSTDmapLabelToCoefficient() const
{
	return m_mapLabelToCoefficient;
}

GAQ gaq::Up(long double _x, long double _y, long double _z)
{
	GAQ x{(_x * e1) + (_y * e2) + (_z * e3)};
	x = eo1 + x + 0.5 * (_x * _x + _y * _y + _z * _z) * ei1 + 0.5 * (_x * _x - _y * _y) * ei2 + 0.5 * (_x * _x - _z * _z) * ei3 + _x * _y * ei4 + _x * _z * ei5 + _y * _z * ei6;
	return x;
}

GAQ gaq::MakeQuadric(long double vo6, long double vo5, long double vo4, long double vo3, long double vo2, long double vo1, long double ve1, long double ve2, long double ve3, long double vi1)
{
	return vo6 * eo6 + vo5 * eo5 + vo4 * eo4 + vo3 * eo3 + vo2 * eo2 + vo1 * eo1 + ve1 * e1 + ve2 * e2 + ve3 * e3 + vi1 * ei1;
}

long double GAQ::ToNumeric() // returns coefficient at basis vector 1
{
	if (m_mapLabelToCoefficient.find("1") != m_mapLabelToCoefficient.end())
	{
		return m_mapLabelToCoefficient.at("1");
	}
	std::cout << "WARNING! ToNumeric found no zero-degree basis blade, returned number 1";
	return 1.0;
}

bool GAQ::IsEqual(const GAQ& other, double precision) const
{
	bool equal{true};
	for (const auto& [basisBlade, value] : m_mapLabelToCoefficient)
	{
		if (other.m_mapLabelToCoefficient.find(basisBlade) == other.m_mapLabelToCoefficient.end()) //if there is on the right not the same basis blade as on the left
		{
			equal = false;
			break;
		}
		//if there is, check for coefs if they are the same
		if (abs(other.m_mapLabelToCoefficient.at(basisBlade) - value) > 1.0 / precision)
		{
			equal = false;
			break;
		}
	}
	for (const auto& [basisBlade, value] : other.m_mapLabelToCoefficient) //now vice versa
	{
		if (m_mapLabelToCoefficient.find(basisBlade) == m_mapLabelToCoefficient.end())
		{
			equal = false;
			break;
		}
		//if there is, check for coefs if they are the same
		if (abs(other.m_mapLabelToCoefficient.at(basisBlade) - value) > 1.0 / precision)
		{
			equal = false;
			break;
		}
	}
	return equal;
}

//************************************MEMBER_OPERATOR************************************\\

GAQ GAQ::RotorExponential(unsigned int degree, long double phi) const
{
	GAQ res{one};
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

GAQ GAQ::TranslatorExponential(unsigned int degree, long double distance) const
{
	GAQ res{one};
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

GAQ& GAQ::operator=(const GAQ& other)
{
	m_mapLabelToCoefficient = other.m_mapLabelToCoefficient;
	return *this;
}

GAQ& GAQ::operator=(GAQ&& other) noexcept
{
	m_mapLabelToCoefficient = std::move(other.m_mapLabelToCoefficient);
	return *this;
}

bool GAQ::operator==(const GAQ& other) const
{
	return this->IsEqual(other, 1e9);
}

bool GAQ::operator!=(const GAQ& other) const
{
	return !(*this == other);
}

//Grade projection operator
GAQ GAQ::operator[](size_t _grade) const
{
	const std::vector<GAQ> left{MakeQCGAFromBasisBlades(*this)};

	GAQ res;
	for (size_t i = 0; i < m_mapLabelToCoefficient.size(); i++)
	{
		res = res + (left[i])(_grade); //equivalent to standard formula 
	}
	res.DeleteZeroFromVector();
	return res;
}

GAQ GAQ::operator[](const GAQ& other) const
{
	const auto it = std::ranges::find_if(m_mapLabelToCoefficient,
		[&](const auto& pair)
		{
			const auto& [key, coeff] = pair;
			return other.m_mapLabelToCoefficient.contains(key);
		});
	if (it != m_mapLabelToCoefficient.end())
	{
		return it->second * GAQ(it->first);
	}
	return {};
}

//geometric product operator
GAQ GAQ::operator*(const GAQ& other) const
{
	std::map<std::string, long double> map;

	for (const auto& [thisBasisBlade, thisCoef] : m_mapLabelToCoefficient) //e1+e1*e2*e3
	{
		for (const auto& [rightBasisBlade, rightCoef] : other.m_mapLabelToCoefficient) //2e2-e3*e4
		{
				map[thisBasisBlade + "*" + rightBasisBlade] = thisCoef * rightCoef; // e1*e2 | e1*e3*e4 | e1*e2*e3*e2 | e1*e2*e3*e3*e4
		}																			//     2 |  	 -1 |			2 |				-1
	}
																				
	GAQ res;
	// try to simplify individual label
	for (const auto& [basisBlade, coef] : map)
	{ 
		std::string_view labelView{basisBlade};

		const char* data = basisBlade.data();
		size_t size = basisBlade.size();

		if (basisBlade.find('e') == std::string::npos)
		{
			labelView = "1";
		}
		else
		{
			const char* start = data;
			const char* end = data + size;
			if (size >= 2 && data[0] == '1' && data[1] == '*')
			{
				start += 2;
			}
			const char* star = static_cast<const char*>(std::memchr(data, '*', size));
			if (!star) star = end;
			if (star + 1 < end && *(star + 1) == '1')
			{
				end = star;
			}
			labelView = std::string_view(start, end - start);
		}

		int sign = 1;
		std::string simplified;
		if (labelView != "1")
			simplified = SimplifyBasisBlade(labelView, sign);
		else
			simplified = "1";

		res = res + GAQ(std::make_pair(std::move(simplified), coef * sign));
	}
	res.DeleteZeroFromVector();
	return res;
}

GAQ GAQ::operator*(GAQ&& other) const
{
	std::map<std::string, long double> map;
	for (const auto& [thisBasisBlade, thisCoef] : m_mapLabelToCoefficient)
	{
		for (auto&& [rightBasisBlade, rightCoef] : other.m_mapLabelToCoefficient)
		{
			map[std::move(thisBasisBlade + "*" + rightBasisBlade)] = std::move(thisCoef * rightCoef);
		}
	}

	GAQ res;
	// try to simplify individual label
	for (const auto& [basisBlade, coef] : map) //each label in newLabel will be modified
	{
		std::string_view labelView{basisBlade};

		const char* data = basisBlade.data();
		size_t size = basisBlade.size();

		if (basisBlade.find('e') == std::string::npos)
		{
			labelView = "1";
		}
		else
		{
			const char* start = data;
			const char* end = data + size;
			if (size >= 2 && data[0] == '1' && data[1] == '*')
			{
				start += 2;
			}
			const char* star = static_cast<const char*>(std::memchr(data, '*', size));
			if (!star) star = end;
			if (star + 1 < end && *(star + 1) == '1')
			{
				end = star;
			}
			labelView = std::string_view(start, end - start);
		}

		int sign = 1;
		std::string simplified;
		if (labelView != "1")
			simplified = SimplifyBasisBlade(labelView, sign);
		else
			simplified = "1";

		res = res + GAQ(std::make_pair(std::move(simplified), coef * sign));
	}
	res.DeleteZeroFromVector();
	return res;
}

//multiplying by scalar from the right
GAQ GAQ::operator*(long double scalar) const
{
	std::map<std::string, long double> map{m_mapLabelToCoefficient};
	for (auto& [_, coef] : map)
	{
		coef *= scalar;
	}
	return GAQ(std::move(map));
}

//reverse operator
GAQ GAQ::operator~() const
{
	std::map<std::string, long double> map;
	for (const auto& [basisBlade, coef] : m_mapLabelToCoefficient) //it is sufficient to check sign of reversed permutation of basis blades
	{
		if (basisBlade == "1")
		{
			map["1"] = coef;
		}
		else
		{
			int count{0};
			int permutation[15]{};
			ExtractIntegersFromBasisBlades(basisBlade, permutation, count);
			std::reverse(permutation, permutation + count);
			map[basisBlade] = coef * CalculateSign(permutation, count);
		}
	}
	return GAQ(std::move(map));
}

//addition operator
GAQ GAQ::operator+(const GAQ& other) const
{
	std::map<std::string, long double> map{m_mapLabelToCoefficient};
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
	return GAQ(map);
}

//substraction operator
GAQ GAQ::operator-(const GAQ& other) const
{
	return *this + ((-1) * other);
}

//inner product operator
GAQ GAQ::operator|(const GAQ& other) const
{
	//Make vectors from left and right operadns, they will store basis blades as CGA object and they will participate in product
	const std::vector<GAQ> left{MakeQCGAFromBasisBlades(*this)};
	const std::vector<GAQ> right{MakeQCGAFromBasisBlades(other)};
	GAQ res;
	for (size_t i = 0; i < m_mapLabelToCoefficient.size(); i++)
	{
		for (size_t j = 0; j < other.m_mapLabelToCoefficient.size(); j++)
		{
			res = res + (left[i] || right[j]); //equivalent to standard formula 
		}
	}
	res.DeleteZeroFromVector();
	return res;
}

//outer product operator
GAQ GAQ::operator^(const GAQ& other) const
{
	//Make vectors from left and right operadns, they will store basis blades as CGA object and they will participate in product
	const std::vector<GAQ> left{MakeQCGAFromBasisBlades(*this)};
	const std::vector<GAQ> right{MakeQCGAFromBasisBlades(other)};

	GAQ res;
	for (size_t i = 0; i < m_mapLabelToCoefficient.size(); i++)
	{
		for (size_t j = 0; j < other.m_mapLabelToCoefficient.size(); j++)
		{
			res = res + (left[i] && right[j]); //equivalent to standard formula 
		}
	}
	res.DeleteZeroFromVector();
	return res;
}

GAQ GAQ::operator^(int exponent) const
{
	if (exponent < 0)
	{
		std::cout << "Warning, calling an inverse of GAQ, not of a blade, this might fail!\n";
		const long double denominator{(((*this) * (~(*this))).ToNumeric())};
		return denominator * (~(*this));;
	}
	GAQ res{*this};
	for (int i = 0; i < exponent - 1; i++)
	{
		res = res * *this;
	}
	res.DeleteZeroFromVector();
	return res;
}

GAQ GAQ::operator/(long double divider) const
{
	if (divider == 0)
	{
		std::cout << "Warning, division by zero detected!" << std::endl;
		return *this;
	}
	return *this * (1.0 / divider);
}

GAQ GAQ::ScalarProduct(const GAQ& b) const
{
	return (*this * b)[0];
}


GAQ GAQ::Rotate(const GAQ& point, rotation_planes plane, long double angle)
{
	GAQ rotor;
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

GAQ GAQ::Translate(const GAQ& point, translation_directions direction, long double distance)
{
	switch (direction)
	{
	case x:
		return Tx(distance) * point * ~Tx(distance);
	case y:
		return Ty(distance) * point * ~Ty(distance);
	case z:
		return Tz(distance) * point * ~Tz(distance);
	default:
		std::cout << "Wrong direction for translating" << std::endl;
		return point;
	}
}

//returns Grade of basis blade (if we give it appropriate label...)
size_t GAQ::Grade(std::string_view label) const
{
	return std::ranges::count(label, 'e');
}

//returns multivector in e_inf, e_o basis, used in << operator
std::string GAQ::Log() const
{
	if (m_mapLabelToCoefficient.find("1") != m_mapLabelToCoefficient.end() && m_mapLabelToCoefficient.at("1") == 0)
	{
		return "0";
	}
	std::string coef;
	std::string s;
	for (auto& a : m_mapLabelToCoefficient)
	{
		if (abs(a.second) < (1000000.0 / PRECISION))
		{
			continue;
		}
		coef = std::to_string(a.second);
		size_t j{coef.length() - 1};
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
	return s.substr(0, s.size() - 3);
}

//************************************PRROTECTED************************************\\

//calculates sign of permutation
int GAQ::CalculateSign(int* permutation, int count) 
{
	bool neg{false};
	for (int i = 0; i < count; ++i)
	{
		for (int j = i + 1; j < count; ++j)
		{
			if (permutation[i] > permutation[j])
			{
				neg = !neg;
			}
		}
	}
	return neg ? -1 : 1;
}

//simplifies label in a form of for example  e1e2e3e2e3 into e1
std::string GAQ::SimplifyBasisBlade(std::string_view label, int& sign)
{
	std::vector<int> permutations; //keeps numbers next to individual e's
	permutations.reserve(30);
	

	size_t pos = 0;
	while (pos < label.size()) 
	{
		size_t epos = label.find('e', pos);
		size_t start = epos + 1;
		size_t end = start;
		while (end < label.size() && std::isdigit(static_cast<unsigned char>(label[end])))
			++end;

		int val = 0;
		auto [ptr, ec] = std::from_chars(label.data() + start, label.data() + end, val);
		if (ec != std::errc()) {
			break;
		}
		permutations.push_back(val);

		pos = end;
		if (pos < label.size() && label[pos] == '*')
			++pos;
	}


	ProcessVector(permutations, sign); //e1e2e5e2e3e4e5 -> e1e5e3e4e5 -> e1e3e4 represented by numbers (1252345 -> 15345 ...)

	std::string res;
	for (size_t i = 0; i < permutations.size(); i++) //creates new proper label
	{
		res += "e" + std::to_string(permutations[i]) + "*";
	}
	if (res.size() == 0)
	{
		res += "1";
	}
	else
	{
		res = res.substr(0, res.length() - 1);
	}
	return res;
}

//used when simplifying results of geometric product: e1e2e5e2e3e4e5 -> e1e5e3e4e5 -> e1e3e4 represented byjust numbers (1252345 -> 15345 ...)
void GAQ::ProcessVector(std::vector<int>& vec, int& sign)
{
	for (size_t i = 0; i < vec.size() - 1; ++i)
	{
		for (size_t j = i + 1; j < vec.size(); ++j)
		{
			if (vec[i] != vec[j])
			{
				continue;
			}
			if (vec[i] > ALGEBRA_P)
			{
				sign *= -1;
			}
			for (size_t k = 0; k < j - i - 1; ++k)
			{
				std::swap(vec[j - k - 1], vec[j - k]);
				sign *= -1;
			}
			vec.erase(vec.begin() + i, vec.begin() + i + 2);
			if (vec.empty())
			{
				return;
			}
			ProcessVector(vec, sign);
		}
		if (vec.empty())
		{
			return;
		}
	}
	if (vec.empty())
	{
		return;
	}
	//bubble sort, easy to track swaps... e3e2e1 -> e1e2e3
	size_t i, j;
	for (i = 0; i < vec.size() - 1; ++i)
	{
		for (j = 0; j < vec.size() - i - 1; ++j)
		{
			if (vec[j] > vec[j + 1])
			{
				std::swap(vec[j], vec[j + 1]);
				sign *= -1;
			}
		}
	}
}

//from a given label, for example e1*e2*e3, returns {1,2,3}
void GAQ::ExtractIntegersFromBasisBlades(std::string_view label, int out_buffer[15], int& out_count)
{
	int pos{0};
	out_count = 0;
	const size_t labelSize{label.size()};
	while (pos < labelSize) // we are guearanteed to have less than 15 elements. A bit dangerous ngl
	{
		++pos;
		int value{0};
		while (pos < labelSize && std::isdigit(label[pos]))
		{
			value = value * 10 + (label[pos++] - '0');
		}
		out_buffer[out_count++] = value;

		if (pos < labelSize && label[pos] == '*')
		{
			++pos;
		}
	}
}

//**********************************STATIC_SUPPORT_OPERATORS**********************************\\

//inner product of two basis blades operator. Usefull in inner product of two general multivectors
GAQ GAQ::operator||(const GAQ& other) const
{
	return GAQ((*this * other)[abs_diff(this->Grade(m_mapLabelToCoefficient.begin()->first), other.Grade(other.m_mapLabelToCoefficient.begin()->first))]);
}

//outer product of two basis blades operator
GAQ GAQ::operator&&(const GAQ& other) const
{
	return GAQ((*this * other)[abs_sum(this->Grade(m_mapLabelToCoefficient.begin()->first), other.Grade(other.m_mapLabelToCoefficient.begin()->first))]);
}

//Grade projection of basis blade
GAQ GAQ::operator()(size_t grade) const
{
	return (this->Grade(m_mapLabelToCoefficient.begin()->first) == grade) ? *this : GAQ("0");
}

void GAQ::DeleteZeroFromVector()
{
	if (auto it = m_mapLabelToCoefficient.find("1");
		it != m_mapLabelToCoefficient.end() &&
		it->second == 0 &&
		m_mapLabelToCoefficient.size() > 1)
	{
		m_mapLabelToCoefficient.erase(it);
	}
}

//**********************************NON-MEMBER_OPERATORS**********************************\\

//multiplying by scalar from the left
GAQ operator*(long double scalar, const GAQ& other)
{
	return other * scalar;
}

// Operator for printing
std::ostream& operator<<(std::ostream& stream, const GAQ& multivector)
{
	stream << multivector.Log();
	return stream;
}

//returns vector of basis blades in linear combination of general multivector
std::vector<GAQ> MakeQCGAFromBasisBlades(const GAQ& multivector)
{
	std::vector<GAQ> basisBlades;
	basisBlades.reserve(20);
	for (const auto& [basisBlade, coef] : multivector.GetSTDmapLabelToCoefficient())
	{
		basisBlades.emplace_back(std::make_pair(basisBlade, coef));
	}
	return basisBlades;
}