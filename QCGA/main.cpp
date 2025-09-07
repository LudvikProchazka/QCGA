#include "Examples.h"

#include <chrono>

int main()
{
	QCGA::GenerateGeneratingBlades(); //create an array of basis vectors in R^{9,6}... one, e1,e2,...,e15 - has to be called first!
	
	/* --------------------------USAGE--------------------------
	QCGA multivecorExample1 = 2 * one + 3 * e1 * e2 - e3 + e11 - 3 * ei1 + 2 * eo6;
	QCGA multivecorExample2 = -1 * one + 2 * e1 * e5 - e9 + e15 - 3 * ei3 + 2 * eo2;
	QCGA multivecorExample3 = multivecorExample1 * multivecorExample2;
	QCGA multivecorExample4 = 3 * multivecorExample2 + eo1;
	*/
	
	RotorXY();
	RotorXZ();
	RotorYZ();
	TranslatorX();
	TranslatorY();
	TranslatorZ();
}
