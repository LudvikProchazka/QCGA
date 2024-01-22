#include "QCGA.h"
#include "Blade.h"

#include <iostream>




int main()
{
	QCGA::generateGeneratingBlades();
	Blade a = up(0, 0, 0);
	Blade b = up(1, 1, 1);
	std::cout << (a | b) << std::endl;
	std::cout << (b | a) << std::endl;

	return 0;
}
