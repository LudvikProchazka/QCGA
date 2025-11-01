#include "Examples.h"

#include <chrono>

void SpeedTest()
{
	const GAQ r1 = e1 ^ e2;
	const GAQ r2 = eo6 ^ ei5;
	const GAQ r3 = ei6 ^ eo5;
	const GAQ r4 = 2 * (eo4 ^ ei2);
	const GAQ r5 = 2 * (ei4 ^ eo2);
	const GAQ r6 = eo4 ^ ei3;

	const GAQ r = r1 + r2 + r3 + r4 + r5 + r6;

	const double phi = std::numbers::pi / 4.0;
	const GAQ R1 = cos(phi / 2) * one + sin(phi / 2) * r1;
	const GAQ R2 = cos(phi / 2) * one + sin(phi / 2) * r2;
	const GAQ R3 = cos(phi / 2) * one + sin(phi / 2) * r3;
	const GAQ R4 = cos(phi) * one + sin(phi) * (0.5 * r4);
	const GAQ R5 = cos(phi) * one + sin(phi) * (0.5 * r5);


	const GAQ C = Up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1-1.5*ei2-4*ei3+2*ei4+3*ei5+6*ei6
	const GAQ c_cga = eo1 + e1 + 2 * e2 + 3 * e3 + 7 * ei1;
	const GAQ c_24 = -1.5 * ei2 + 2 * ei4;
	const GAQ c_56 = 3 * ei5 + 6 * ei6;
	const GAQ c_3 = -4 * ei3;

	const double x = 1;
	const double y = 2;
	const double z = 3;
	const double theta = atan(2) - phi;
	const GAQ target = Up(sqrt(5) * cos(theta), sqrt(5) * sin(theta), 3);
	const GAQ _rotated = (R1 * c_cga * ~R1) + ((R2 ^ R3) * c_56 * (~R3 ^ ~R2)) + c_3 + ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + (-0.5 * sin(phi) * sin(phi) * (x * x - y * y) + sin(phi) * cos(phi) * x * y) * ei3;

	const GAQ rotor = r.RotorExponential(20, phi);
	const GAQ rotated = (rotor * C * ~rotor)[1];

	std::cout << (rotated == target) << std::endl;
}

static size_t allocations{0};

void* operator new(size_t bytes)
{
	++allocations;
	return malloc(bytes);
}

int main()
{
	GAQ::GenerateGeneratingBlades(); //create an array of basis vectors in R^{9,6}... one, e1,e2,...,e15 - has to be called first!
	
	/* --------------------------USAGE--------------------------
	GAQ multivecorExample1 = 2 * one + 3 * e1 * e2 - e3 + e11 - 3 * ei1 + 2 * eo6;
	GAQ multivecorExample2 = -1 * one + 2 * e1 * e5 - e9 + e15 - 3 * ei3 + 2 * eo2;
	GAQ multivecorExample3 = multivecorExample1 * multivecorExample2;
	GAQ multivecorExample4 = 3 * multivecorExample2 + eo1;
	*/
	using namespace std::chrono;

	auto start = high_resolution_clock::now();


	const double phi = std::numbers::pi / 4.0;
	const double theta = atan(2) - phi;
	const GAQ XY = Up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1-1.5*ei2-4*ei3+2*ei4+3*ei5+6*ei6
	const GAQ XZ = Up(1, 3, 2); //eo1+e1+3*e2+2*e3+7*ei1-2*ei2-3*ei3+3*ei4+2*ei5+6*ei6
	const GAQ YZ = Up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1-1.5*ei2-4*ei3+2*ei4+3*ei5+6*ei6
	const GAQ rotatedXY = GAQ::Rotate(XY, xy, phi);
	const GAQ rotatedXZ = GAQ::Rotate(XZ, xz, phi);
	const GAQ rotatedYZ = GAQ::Rotate(YZ, yz, phi);
	const GAQ targetXY = Up(sqrt(5) * cos(theta), sqrt(5) * sin(theta), 3);
	const GAQ targetXZ = Up(sqrt(5) * cos(theta), 3, sqrt(5) * sin(theta));
	const GAQ targetYZ = Up(1, 0.5 * 5 * sqrt(2), 0.5 * sqrt(2));

	auto end = high_resolution_clock::now();

	auto duration = duration_cast<microseconds>(end - start);

	std::cout << duration.count() << "ms = " << duration.count()/1000000.0 << " sec, Allocations: " << allocations << std::endl;
	//RotorXZ();
	//RotorYZ();
	//TranslatorX();
	//TranslatorY();
	//TranslatorZ();
}
