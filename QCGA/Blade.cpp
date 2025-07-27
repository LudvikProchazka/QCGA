#include "Blade.h"
#include <unordered_set>

bool Blade::IsBlade() 
{
	if ((*this)[0] == *this)
	{
		return true;
	}
	QCGA l(*this * ~(*this));
	QCGA p = l[0];
	return ((*this * ~(*this)) == (*this * ~(*this))[0]);
}

//creates blade from given Multivector
Blade::Blade(const QCGA& Multivector) : QCGA(Multivector)
{

	if (IsBlade())
	{
		m_grade = QCGA::grade(m_mapLabelToCoefficient.begin()->first);
	}
	else
	{
		m_grade = -1;
		std::cout << "WARNING: Multivector: " << Multivector << " is not a blade!, program will likely crash if it is being used as blade" << std::endl;
	}
	m_isNullBlade = (Multivector | Multivector) == zero_vector;
}

int Blade::getGrade() const
{
	return m_grade;
}

bool Blade::isNullBlade() const
{
	return m_isNullBlade;
}

Blade Blade::operator^(const Blade& other) const
{
	return Blade(static_cast<QCGA>(*this) ^ other);
}

//returns inverse and exponent
Blade Blade::operator^(const int exponent) const
{
	if (exponent < 0)
	{
		if (this->m_isNullBlade)
		{
			std::cout << "WARNING, Blade:" << *this << " is a null-blade, cant make inversion! returned with positive exponent" << std::endl;
			return static_cast<QCGA>(*this) ^ (-1 * exponent);
		}
		else
		{
			Blade res = (~*this) / ((*this * ~(*this)).ToNumeric());
			res = static_cast<QCGA>(res) ^ (-1 * exponent);
			return res;
		}
	}
	else
	{
		Blade _res = static_cast<QCGA>(*this) ^ exponent;
		return _res;
	}
}

Blade Blade::dual() const
{
	return Blade(*this * (I ^ (-1)));
}

Blade Blade::normalize() const
{
	double multiplicator = ((*this | (ei1 * ei2 * ei3 * ei4 * ei5 * ei6)) | (eo2 * eo3 * eo4 * eo5 * eo6)).ToNumeric();
	Blade res = (1.0 / multiplicator) * *this;
	return res;
}

Blade Blade::down() const
{
	Blade res = (*this).normalize();
	return res[e1] + res[e2] + res[e3];
}

//creates QCGA object as embedded 3D point
Blade up(long double _x, long double _y, long double _z)
{
	Blade x = (_x * e1) + (_y * e2) + (_z * e3); //eucledian point
	x = eo1 + x + 0.5 * (_x * _x + _y * _y + _z * _z) * ei1 + 0.5 * (_x * _x - _y * _y) * ei2 + 0.5 * (_x * _x - _z * _z) * ei3 + _x * _y * ei4 + _x * _z * ei5 + _y * _z * ei6;
	return x;
}

Blade makeQuadric(long double vo6, long double vo5, long double vo4, long double vo3, long double vo2, long double vo1, long double ve1, long double ve2, long double ve3, long double vi1)
{
	return Blade(vo6 * eo6 + vo5 * eo5 + vo4 * eo4 + vo3 * eo3 + vo2 * eo2 + vo1 * eo1 + ve1 * e1 + ve2 * e2 + ve3 * e3 + vi1 * ei1);
}

