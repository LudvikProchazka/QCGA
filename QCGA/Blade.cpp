#include "Blade.h"
#include <unordered_set>

// Used in constructor, A is blade <=> A*~A is scalar (i think)
bool Blade::isBlade(QCGA& Multivector) {
	if (Multivector[0] == Multivector)
	{
		return true;
	}
	else {
		QCGA l(Multivector * ~Multivector);
		QCGA p = l[0];
		return ((Multivector * ~Multivector) == (Multivector * ~Multivector)[0]);
	}
}

//creates blade from given Multivector
Blade::Blade(const QCGA& Multivector) : QCGA(Multivector)
{

	if (isBlade(*this))
	{
		this->grade = QCGA::grade(this->STDmapLabelToCoefficient.begin()->first);
	}
	else
	{
		this->grade = -1;
		std::cout << "WARNING: Multivector: " << Multivector << " is not a blade!, program will likely crash if it is being used as blade" << std::endl;
	}
	this->nullBlade = (Multivector | Multivector) == zero;
}

Blade Blade::operator^(const Blade& other) const
{
	return Blade((QCGA)*this ^ other);
}

//returns inverse and exponent
Blade Blade::operator^(const int exponent) const
{
	if (exponent < 0)
	{
		if (this->nullBlade)
		{
			std::cout << "WARNING, Blade:" << *this << " is a null-blade, cant make inversion! returned with positive exponent" << std::endl;
			return (QCGA)*this ^ (-1 * exponent);
		}
		else
		{
			Blade res = (~*this) / ((*this * ~(*this)).toNumeric());
			res = (QCGA)res ^ (-1 * exponent);
			return res;
		}
	}
	else
	{
		Blade _res = (QCGA)*this ^ exponent;
		return _res;
	}
}

// returns dual blade
Blade Blade::dual()
{
	return Blade(*this * (I ^ (-1)));
}

//creates CGA object as embedded 3D point
Blade up(float _x, float _y, float _z)
{
	Blade x = (_x * e1) + (_y * e2) + (_z * e3); //eucledian point
	x = x + 0.5 * ((_x * _x) * ei1 + (_y * _y) * ei2 + (_z * _z) * ei3) + (_x * _y) * ei4 + (_x * _z) * ei5 + (_y * _z) * ei6 + eo1 + eo2 + eo3;
	return x;
}