#ifndef _BLADE_H_
#define _BLADE_H_

#include "QCGA.h"

class Blade : public QCGA
{
public:
	Blade(const QCGA& Multivector); //creates blade from given Multivector, if possible. If given multivector is not blade, warning occurs and program might crash

	int getGrade() { return this->grade; };
	bool isNullBlade() { return this->nullBlade; };

	//**********************************OPERATORS**********************************\\

	Blade operator^(const Blade& other) const;
	Blade operator ^(const int exponent) const; //exponent operator, mainly for inverse: A^(-1)
	Blade dual(); //dual: A.dual() = A * I^(-1)
	Blade down() const;

private:
	int grade;
	bool nullBlade;
	static bool isBlade(QCGA& Multivector);  // Used in constructor, A is blade <=> A*~A is scalar. Well, i hope so
};
Blade up(long double x, long double y, long double z); //Embedding of a 3D point. This has to be modified when changing algebra
Blade MujUp(long double x, long double y, long double z); //Embedding of a 3D point. This has to be modified when changing algebra
#endif