#include "QCGA.h"
#include "Blade.h"

#include <iostream>




int main()
{
	QCGA::generateGeneratingBlades();
	Blade a = e1 + e2;
	std::cout << a << std::endl;
	std::cout << (a^(-1)) << std::endl;
	std::cout << (a* (a ^ (-1))) << std::endl;
	return 0;
}
