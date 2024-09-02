#include "QCGA.h"
#include "Blade.h"

#include <iostream>
#include <numbers>


void com(const QCGA& a, const QCGA& b)
{
	std::cout << 0.5 * ((a * b) - (b * a)) << std::endl;
}

void case1()
{
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
	QCGA R4 = cos(phi) * one + sin(phi) * (0.5 * r4);
	QCGA R5 = cos(phi) * one + sin(phi) * (0.5 * r5);
	QCGA R6 = cos(phi / 2) * one + sin(phi / 2) * r6;
	QCGA R7 = cos(phi / 2) * one + sin(phi / 2) * r7;
	
	
	QCGA A = up(1, 1, 0);
	QCGA C = up(1, 1, 1);
	QCGA D = up(1, -1, 1);
	QCGA B = up(sqrt(2), 0, 0);
	QCGA CC = up(sqrt(2), 0, 1);
	
	std::cout << "r2r3: " << r2 * r3 - r3 * r2 << std::endl;
	std::cout << "r4r5: " << r4 * r5 - r5 * r4 << std::endl;
	std::cout << "r4r6: " << r4 * r6 - r6 * r4 << std::endl;
	std::cout << "r5r7: " << r5 * r7 - r7 * r5 << std::endl;
	
	std::cout << "************************" << std::endl;
	std::cout << (2 * e5 ^ e7) * (r4 + r5) << std::endl;
	std::cout << (2 * e13 ^ e11) * (r4 + r5) << std::endl;
	std::cout << -1 * one + 0.5 * ((r4 + r5) ^ (r4 + r5)) << std::endl;
	
	QCGA r1_ = e1 ^ e2;
	QCGA r2_ = e8 ^ e9;
	QCGA r3_ = e15 ^ e14;
	QCGA r4_ = 2 * (e5 ^ e7);
	QCGA r5_ = 2 * (e13 ^ e11);
	QCGA r6_ = ei4 ^ eo3;
	QCGA r7_ = eo1 ^ ei2;
	
	QCGA R1_ = cos(phi / 2) * one + sin(phi / 2) * r1;
	QCGA R2_ = cos(phi / 2) * one + sin(phi / 2) * r2_;
	QCGA R3_ = cos(phi / 2) * one + sin(phi / 2) * r3_;
	QCGA R4_ = cos(phi) * one + sin(phi) * (0.5 * r4_);
	QCGA R5_ = cos(phi) * one + sin(phi) * (0.5 * r5_);
	QCGA R6_ = cos(phi / 2) * one + sin(phi / 2) * r6_;
	QCGA R7_ = cos(phi / 2) * one + sin(phi / 2) * r7_;
	
	QCGA bod = ei1 + e1 + e2 + e3 + 1.5 * ei1 + 0.5 * ei2 + 0.5 * ei3 + ei4 + ei5 + ei6;
	
	std::cout << "************************" << std::endl;
	//std::cout << "rotA: " << R1*(R2*(R3*(R4*(R5*(R6*(R7*(eo1+1.5*ei1)*(~R7))*(~R6))*(~R5))*(~R4))*(~R3))*(~R2))*(~R1) << std::endl;
	//std::cout << "rotA: " << R1*((R2^R3)*((R4^R5)*((R6^R7)*C*((~R7)^(~R6)))*((~R5)^(~R4)))*((~R3)^(~R2)))*(~R1) << std::endl;
	//std::cout << ": " << R1* (ei1 + e1 + e2 + e3 + 1.5 * ei1) *(~R1) << std::endl;
	std::cout << ": " << ei2 << std::endl;
	std::cout << ": " << ((R6) ^ (R7)) * (ei2) * ((~R7) ^ (~R6)) << std::endl;
	std::cout << "CC: " << CC << std::endl;
}	

void case2()
{
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
	QCGA R4 = cos(phi) * one + sin(phi) * (0.5 * r4);
	QCGA R5 = cos(phi) * one + sin(phi) * (0.5 * r5);
	QCGA R6 = cos(phi / 2) * one + sin(phi / 2) * r6;
	QCGA R7 = cos(phi / 2) * one + sin(phi / 2) * r7;
	
	
	QCGA A = up(1, 1, 0);
	QCGA a_euc = eo1 + e1 + e2 + ei1;
	QCGA a_24 = ei4;
	QCGA a_3 = ei3;
	QCGA a_56 = zero;

	QCGA C = up(1, 1, 1); //eo1+e1+e2+e3+1.5*ei1+0.5*ei2+0.5*ei3+ei4+ei5+ei6
	QCGA c_euc = eo1 + e1 + e2 + e3 + 1.5 * ei1;
	QCGA c_24 =0.5 * ei2 + 1*ei4;
	QCGA c_3 = 0.5 * ei3;
	QCGA c_56 = ei5 + ei6;

	QCGA B = up(sqrt(2), 0, 0);
	QCGA CC = up(sqrt(2), 0, 1);
	std::cout << "      C: " << C << std::endl;
	//std::cout << "Rotated: " << (R1 * a_euc * ~R1) + ((R4 ^ R5) * a_24 * (~R5 ^ ~R4)) + a_3 << std::endl;
	//std::cout << " Target: " << B << std::endl;
	std::cout << "Rotated: " << (R1 * c_euc * ~R1) + ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + c_3 + ((R2^R3)*c_56*(~R3 ^ ~R2))<< std::endl;
	std::cout << "Rotated: " << (R1 * c_euc * ~R1) + ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + c_3 + ((R2 ^ R3) * c_56 * (~R3 ^ ~R2)) + (0.5*(one -cos(2*phi)*one)*1*ei2) + ((0.5*sin(2*phi)*one)*ei4)<< std::endl;
	
	std::cout << " Target: " << CC << std::endl;

}

void case3()
{
	QCGA r1 = e1 ^ e2;
	QCGA r2 = eo6 ^ ei5;
	QCGA r3 = ei6 ^ eo5;
	QCGA r4 = 2 * (eo4 ^ ei2);
	QCGA r5 = 2 * (ei4 ^ eo2);
	QCGA r6 = ei4 ^ eo3;
	QCGA r7 = eo1 ^ ei2;

	double phi = std::numbers::pi / 2;

	QCGA R1 = cos(phi / 2) * one + sin(phi / 2) * r1;
	QCGA R2 = cos(phi / 2) * one + sin(phi / 2) * r2;
	QCGA R3 = cos(phi / 2) * one + sin(phi / 2) * r3;
	QCGA R4 = cos(phi) * one + sin(phi) * (0.5 * r4);
	QCGA R5 = cos(phi) * one + sin(phi) * (0.5 * r5);
	QCGA R6 = cos(phi / 2) * one + sin(phi / 2) * r6;
	QCGA R7 = cos(phi / 2) * one + sin(phi / 2) * r7;


	QCGA C = up(1, 2, 3); //eo1+e1+e2+e3+1.5*ei1+0.5*ei2+0.5*ei3+ei4+ei5+ei6
	QCGA c_euc = eo1 + e1 + 2*e2 + 3*e3 + 7 * ei1;
	QCGA c_24 = 3 * ei2 + 2 * ei4;
	QCGA c_3 = -2 * ei3;
	QCGA c_56 = 3*ei5 + 6*ei6;
	
	double z = 3;
	QCGA B = up(sqrt(2), 0, 0);
	QCGA CC = up(2, -1, 3);
	std::cout << "      C: " << C << std::endl;
	QCGA rot1 = (R1 * c_euc * ~R1) + ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + c_3 + ((R2 ^ R3) * c_56 * (~R3 ^ ~R2));
	QCGA rot2 = (R1 * c_euc * ~R1) + ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + c_3 + ((R2 ^ R3) * c_56 * (~R3 ^ ~R2)) + (0.5 * (one - cos(2 * phi) * one) * z * z * ei2) + ((0.5 * sin(2 * phi) * one) * z * z * ei4);
	std::cout << "Rotated: " << rot1 << std::endl;
	std::cout << "Rotated: " << rot2 << std::endl;
	
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rot2 == CC) << std::endl;

}	

int main()
{	
	//TEST FRO REPO
	QCGA::generateGeneratingBlades();
	case3();

	return 0;
}
