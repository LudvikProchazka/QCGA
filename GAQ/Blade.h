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