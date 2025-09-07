#include "Examples.h"

#include <chrono>

int main()
{
	auto start = std::chrono::high_resolution_clock::now();
	QCGA::GenerateGeneratingBlades(); //create an array of basis vectors in R^{9,6}... one, e1,e2,...,e15

	RotorXY();
	RotorXZ();
	RotorYZ();
	TranslatorX();
	TranslatorY();
	TranslatorZ();

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = duration_cast<std::chrono::milliseconds>(stop - start);

	std::cout << "Execution duration: " << duration.count() << " ms" << std::endl;
}
