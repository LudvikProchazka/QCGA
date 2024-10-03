#ifndef _BLADE_H_
#define _BLADE_H_

#include "QCGA.h"

class Blade : public QCGA
{
public:
	Blade(const QCGA& Multivector); //creates blade from given Multivector, if possible. If given multivector is not blade, warning occurs and program might crash
	virtual ~Blade() = default;
	int getGrade() const { return this->grade; };
	bool isNullBlade() const { return this->nullBlade; };

	//**********************************OPERATORS**********************************\\

	Blade operator^(const Blade& other) const;
	Blade operator^(const int exponent) const; //exponent operator, mainly for inverse: A^(-1)
	Blade dual() const; //dual: A.dual() = A * I^(-1)
	Blade down() const;

private:
	int grade;
	bool nullBlade;
	static bool isBlade(QCGA& Multivector);  // Used in constructor, A is blade <=> A*~A is scalar. Well, i hope so
};

Blade up(long double x, long double y, long double z); //embedding of a 3D point.
Blade makeQuadric(long double vo6, long double vo5,long double vo4,long double vo3,long double vo2,long double vo1,long double ve1,long double ve2,long double ve3,long double vi1); //IPNS representation of a quadric

#endif