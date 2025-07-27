#include "Examples.h"

int main()
{
	QCGA::generateGeneratingBlades(); //create an array of basis vectors in R^{9,6}... one, e1,e2,...,e15
	
	RotorXY();
	RotorXZ();
	RotorYZ();
	TranslatorX();
	TranslatorY();
	TranslatorZ();

}
