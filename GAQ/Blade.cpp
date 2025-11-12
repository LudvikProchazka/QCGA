#include "Blade.h"
#include <unordered_set>

bool Blade::IsBlade() const
{
	if ((*this)[0] == *this)
	{
		return true;
	}
	return ((*this * ~(*this)) == (*this * ~(*this))[0]);
}

//creates blade from given Multivector
Blade::Blade(const GAQ& Multivector) : GAQ(Multivector)
{
	if (IsBlade())
	{
		m_grade = GAQ::Grade(m_mapLabelToCoefficient.begin()->first);
	}
	else
	{
		m_grade = -1;
		std::cout << "WARNING: Multivector: " << Multivector << " is not a blade!, program will likely crash if it is being used as blade" << std::endl;
	}
	m_isNullBlade = (Multivector | Multivector) == zero_vector;
}

size_t Blade::GetGrade() const
{
	return m_grade;
}

bool Blade::IsNullBlade() const
{
	return m_isNullBlade;
}

//returns inverse and exponent
Blade Blade::operator^(const int exponent) const
{
	if (exponent < 0)
	{
		if (m_isNullBlade)
		{
			std::cout << "WARNING, Blade:" << *this << " is a null-blade, cant make inversion! returned with positive exponent" << std::endl;
			return static_cast<GAQ>(*this) ^ (-1 * exponent);
		}
		// inputting -exp results in an inversion and that is exponentiated to the exp power
		Blade res{(~*this) / ((*this * ~(*this)).ToNumeric())};
		res = static_cast<GAQ>(res) ^ (-1 * exponent);
		return res;
	}
	return static_cast<GAQ>(*this) ^ exponent;
}

Blade Blade::Dual() const
{
	return Blade(*this * (I ^ (-1)));
}

Blade Blade::Normalize() const
{
	const long double multiplicator{((*this | (ei1 * ei2 * ei3 * ei4 * ei5 * ei6)) | (eo2 * eo3 * eo4 * eo5 * eo6)).ToNumeric()};
	return (1.0 / multiplicator) * *this;
}

Blade Blade::Down() const
{
	const Blade& res{(*this).Normalize()};
	return res[e1] + res[e2] + res[e3];
}
