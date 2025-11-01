#pragma once
#include "GAQ.h"

class Blade : public GAQ
{
public:
	Blade(const GAQ& Multivector); //creates blade from given Multivector, if possible. If given multivector is not blade, warning occurs and program might crash
	virtual ~Blade() = default;

	size_t GetGrade() const;
	bool IsNullBlade() const;

	//**********************************OPERATORS**********************************\\

	Blade operator^(int exponent) const; //exponent operator, mainly for inverse: A^(-1)
	Blade Dual() const; //Dual: A.Dual() = A * I^(-1)
	Blade Normalize() const;	
	Blade Down() const;	
	bool IsBlade() const; 

private:
	size_t m_grade;
	bool m_isNullBlade;
};

GAQ Up(long double x, long double y, long double z); //embedding of a 3D point.
GAQ MakeQuadric(long double vo6, long double vo5, long double vo4, long double vo3, long double vo2, long double vo1, long double ve1, long double ve2, long double ve3, long double vi1); //IPNS representation of a quadric