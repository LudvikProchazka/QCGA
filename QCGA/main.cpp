#include "QCGA.h"
#include "Blade.h"

#include <iostream>
#include <numbers>


void com(const QCGA& first, const QCGA& second)
{
	std::cout << 0.5 * ((first * second) - (second * first)) << std::endl;
}

int main()
{
	//TEST FRO REPO
	QCGA::generateGeneratingBlades();
	QCGA r1 = e1 ^ e2;
	QCGA r2 = eo6 ^ ei5;
	QCGA r3 = ei6 ^ eo5;
	QCGA r4 = 2 * (eo4 ^ ei2);
	QCGA r5 = 2 * (ei4 ^ eo2);
	QCGA r6 = ei4 ^ eo3;
	QCGA r7 = eo1 ^ ei2;

	double phi = std::numbers::pi / 4;

	QCGA R1 = cos(phi / 2) * one + sin(phi / 2) * r1;
	QCGA R2 = cos(phi / 2) * one + sin(phi / 2) * r2;
	QCGA R3 = cos(phi / 2) * one + sin(phi / 2) * r3;
	QCGA R4 = cos(phi) * one + sin(phi) * r4;
	QCGA R5 = cos(phi) * one + sin(phi) * r5;
	QCGA R6 = cos(phi / 2) * one + sin(phi / 2) * r6;
	QCGA R7 = cos(phi / 2) * one + sin(phi / 2) * r7;

	QCGA R = R1 * R2 * R3 * R4 * R5 * R6 * R7;

	QCGA A = up(1, 1, 0);
	QCGA B = up(sqrt(2), 0, 0);

	QCGA Aeuc = eo1 + e1 + e2 + ei1;
	QCGA Aneeuc = ei3+ei4;

	std::cout << Blade(R1 * (R2 * (R3 * (R4 * (R5 * (R6 * (R7 * A * (~R7)) * (~R6)) * (~R5)) * (~R4)) * (~R3)) * (~R2)) * (~R1)).down() << std::endl;
	std::cout << (R1*Aeuc* (~R1)) + (R2 * (R3 * (R4 * (R5 * (R6 * (R7 * Aneeuc * (~R7)) * (~R6)) * (~R5)) * (~R4)) * (~R3)) * (~R2)) << std::endl;
	std::cout << B << std::endl;
	return 0;
}
