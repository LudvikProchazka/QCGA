#include "QCGA.h"
#include "Blade.h"

#include <iostream>
#include <numbers>
#include <iomanip>      // std::setprecision


QCGA com(const QCGA& a, const QCGA& b)
{
	return 0.5 * ((a * b) - (b * a));
}

void Rotor1XY()
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

void Rotor2XY()
{
	QCGA r1 = e1 ^ e2;
	QCGA r2 = eo6 ^ ei5;
	QCGA r3 = ei6 ^ eo5;
	QCGA r4 = 2 * (eo4 ^ ei2);
	QCGA r5 = 2 * (ei4 ^ eo2);
	QCGA r6 = ei4 ^ eo3;
	QCGA r7 = eo1 ^ ei2;

	QCGA r = r2 + r3;

	double phi = std::numbers::pi / 4;

	QCGA R1 = cos(phi / 2) * one + sin(phi / 2) * r1;
	QCGA R2 = cos(phi / 2) * one + sin(phi / 2) * r2;
	QCGA R3 = cos(phi / 2) * one + sin(phi / 2) * r3;
	QCGA R4 = cos(phi) * one + sin(phi) * (0.5 * r4);
	QCGA R5 = cos(phi) * one + sin(phi) * (0.5 * r5);
	QCGA R6 = cos(phi / 2) * one + sin(phi / 2) * r6;
	QCGA R7 = cos(phi / 2) * one + sin(phi / 2) * r7;


	QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA c_euc = eo1 + e1 + 2 * e2 + 3 * e3 + 7 * ei1;
	QCGA c_24 = 3 * ei2 + 2 * ei4;
	QCGA c_3 = -2 * ei3;
	QCGA c_56 = 3 * ei5 + 6 * ei6;

	double z = 3;
	double theta = atan(2) - phi;
	//QCGA CC = up(2, -1, 3);
	QCGA CC = up(sqrt(5) * cos(theta), sqrt(5) * sin(theta), 3);
	std::cout << "      C: " << C << std::endl;
	QCGA rot1 = (R1 * c_euc * ~R1) + ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + c_3 + ((R2 ^ R3) * c_56 * (~R3 ^ ~R2));
	QCGA rot2 = (R1 * c_euc * ~R1) + ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + c_3 + ((R2 ^ R3) * c_56 * (~R3 ^ ~R2)) + (0.5 * (one - cos(2 * phi) * one) * z * z * ei2) + ((0.5 * sin(2 * phi) * one) * z * z * ei4);
	std::cout << "Rotated: " << rot1 << std::endl;
	std::cout << "Rotated: " << rot2 << std::endl;

	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rot2 == CC) << std::endl;

}

void Rotor3XY()
{
	QCGA r1 = e1 ^ e2;
	QCGA r2 = eo6 ^ ei5;
	QCGA r3 = ei6 ^ eo5;
	QCGA r4 = 2 * (eo4 ^ ei2);
	QCGA r5 = 2 * (ei4 ^ eo2);
	QCGA r6 = ei4 ^ eo3;
	QCGA r7 = eo1 ^ ei2;

	QCGA r = r2 + r3;

	double phi = std::numbers::pi / 4;

	QCGA R1 = cos(phi / 2) * one + sin(phi / 2) * r1;
	QCGA R2 = cos(phi / 2) * one + sin(phi / 2) * r2;
	QCGA R3 = cos(phi / 2) * one + sin(phi / 2) * r3;
	QCGA R4 = cos(phi) * one + sin(phi) * (0.5 * r4);
	QCGA R5 = cos(phi) * one + sin(phi) * (0.5 * r5);
	QCGA R6 = cos(phi / 2) * one + sin(phi / 2) * r6;
	QCGA R7 = cos(phi / 2) * one + sin(phi / 2) * r7;


	QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA c_euc = eo1 + e1 + 2*e2 + 3*e3 + 7 * ei1;
	QCGA c_24 = 3 * ei2 + 2 * ei4;
	QCGA c_3 = -2 * ei3;
	QCGA c_56 = 3*ei5 + 6*ei6;
	
	double z = 3;
	double theta = atan(2) - phi;
	//QCGA CC = up(2, -1, 3);
	QCGA CC = up(sqrt(5)*cos(theta), sqrt(5) * sin(theta), 3);
	std::cout << "      C: " << C << std::endl;
	QCGA rot1 = (R1 * c_euc * ~R1) + ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + c_3 + ((R2 ^ R3) * c_56 * (~R3 ^ ~R2));
	QCGA rot2 = (R1 * c_euc * ~R1) + ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + c_3 + ((R2 ^ R3) * c_56 * (~R3 ^ ~R2)) + (0.5 * (one - cos(2 * phi) * one) * z * z * ei2) + ((0.5 * sin(2 * phi) * one) * z * z * ei4);
	std::cout << "Rotated: " << rot1 << std::endl;
	std::cout << "Rotated: " << rot2 << std::endl;
	
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rot2 == CC) << std::endl;
#pragma region MyRegion
	//com(r1, r2);
	//com(r1, r3);
	//com(r1, r4);
	//com(r1, r5);
	//com(r1, r6);
	//com(r1, r7);
	//std::cout << std::endl;
	//com(r2, r3);
	//com(r2, r4);
	//com(r2, r5);
	//com(r2, r6);
	//com(r2, r7);
	//std::cout << std::endl;
	//com(r3, r4);
	//com(r3, r5);
	//com(r3, r6);
	//com(r3, r7);
	//std::cout << std::endl;
	//com(r4, r5);
	//com(r4, r6);
	//com(r4, r7);
	//std::cout << std::endl;
	//com(r5, r6);
	//com(r5, r7);
	//std::cout << std::endl;
	//com(r6, r7);
	//std::cout << std::endl;
	//std::cout << r1 * r1 << std::endl;
	//std::cout << r2 * r2 << std::endl;
	//std::cout << r3 * r3 << std::endl;
	//std::cout << r4 * r4 << std::endl;
	//std::cout << r5 * r5 << std::endl;
	//std::cout << r6 * r6 << std::endl;
	//std::cout << r7 * r7 << std::endl;
#pragma endregion

	

}	

void Rotor3XTaylor()
{
	QCGA r1 = 0.5*e1 ^ e2;
	QCGA r2 = 0.5*eo6 ^ ei5;
	QCGA r3 = 0.5*ei6 ^ eo5;
	QCGA r4 = 0.5*2 * (eo4 ^ ei2);
	QCGA r5 = 0.5*2 * (ei4 ^ eo2);
	QCGA r6 = 0.5*ei4 ^ eo3;
	QCGA r7 = 0.5*eo1 ^ ei2;

	QCGA r234567 = r2 + r3 + r4 + r5 + r6 + r7;

	double phi = std::numbers::pi / 4;

	QCGA R1 = cos(phi / 2) * one + sin(phi / 2) *2* r1;
	QCGA R2 = cos(phi / 2) * one + sin(phi / 2) *2* r2;
	QCGA R3 = cos(phi / 2) * one + sin(phi / 2) *2* r3;
	QCGA R4 = cos(phi) * one + sin(phi) * (0.5 *2* r4);
	QCGA R5 = cos(phi) * one + sin(phi) * (0.5 *2* r5);
	QCGA R6 = cos(phi / 2) * one + sin(phi / 2) *2* r6;
	QCGA R7 = cos(phi / 2) * one + sin(phi / 2) *2* r7;


	QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA c_euc = eo1 + e1 + 2 * e2 + 3 * e3 + 7 * ei1;
	QCGA c_24 = 3 * ei2 + 2 * ei4;
	QCGA c_3 = -2 * ei3;
	QCGA c_56 = 3 * ei5 + 6 * ei6;


	double z = 3;
	double theta = atan(2) - phi;
	//QCGA CC = up(2, -1, 3);
	QCGA CC = up(sqrt(5) * cos(theta), sqrt(5) * sin(theta), 3);
	std::cout << "      C: " << C << std::endl;
	QCGA rot1 = (R1 * c_euc * ~R1) + ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + c_3 + ((R2 ^ R3) * c_56 * (~R3 ^ ~R2)) + (0.5 * (one - cos(2 * phi) * one) * z * z * ei2) + ((0.5 * sin(2 * phi) * one) * z * z * ei4);
	std::cout << "Rotated: " << rot1 << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rot1 == CC) << std::endl << std::endl;
	QCGA taylor_r1 = r1.rotorExponential(20, phi);
	QCGA taylor_r45 = (r4+r5).rotorExponential(20, phi);
	QCGA taylor_23 = (r2+r3).rotorExponential(20, phi);
	QCGA rot2 = (taylor_r1 * c_euc * ~taylor_r1) + (taylor_23 * c_56 * ~taylor_23)+(taylor_r45 * (c_24) * ~taylor_r45) + c_3 + (0.5 * (one - cos(2 * phi) * one) * z * z * ei2) + ((0.5 * sin(2 * phi) * one) * z * z * ei4);
	std::cout << "Rotated: " << rot2 << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rot2 == CC) << std::endl << std::endl;
	std::cout << taylor_r1 << std::endl << std::endl;
	std::cout << R1 << std::endl << std::endl;
}

void Rotor3XTaylor2()
{
	QCGA r1 = 0.5 * e1 ^ e2;
	QCGA r2 = 0.5 * eo6 ^ ei5;
	QCGA r3 = 0.5 * ei6 ^ eo5;
	QCGA r4 = 0.5 * 2 * (eo4 ^ ei2);
	QCGA r5 = 0.5 * 2 * (ei4 ^ eo2);
	QCGA r6 = 0.5 * ei4 ^ eo3;
	QCGA r7 = 0.5 * eo1 ^ ei2;

	double phi = std::numbers::pi / 4;

	QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA c_euc = eo1 + e1 + 2 * e2 + 3 * e3;
	QCGA c_24 = 3 * ei2 + 2 * ei4;
	QCGA c_3 = -2 * ei3;
	QCGA c_56 = 3 * ei5 + 6 * ei6;


	double z = 3;
	double theta = atan(2) - phi;
	//QCGA CC = up(2, -1, 3);
	QCGA CC = up(sqrt(5) * cos(theta), sqrt(5) * sin(theta), 3);
	std::cout << "      C: " << C << std::endl;
	QCGA taylor_r1 = r1.rotorExponential(20, phi);
	QCGA taylor_r1234 = (r4 + r5 + r6+r7).rotorExponential(20, phi);
	QCGA taylor_56 = (r2 + r3).rotorExponential(20, phi);
	QCGA rot2 = (taylor_r1 * c_euc * ~taylor_r1) + (taylor_56 * c_56 * ~taylor_56) + (taylor_r1234 * (c_24 + (7 * ei1) + c_3) * ~taylor_r1234)[1];
	std::cout << "Rotated: " << rot2 << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rot2 == CC) << std::endl << std::endl;
}

void Rotor1XZ()
{
	QCGA r1 = -0.5*1*e3 ^ e1; //rotating in opposite direction, thus -1*, viz corresponding matrices
	QCGA r2 = -0.5*1*ei4 ^ eo6;
	QCGA r3 = -0.5*2 * (ei3 ^ eo5);
	QCGA r4 = -0.5*1*eo4 ^ ei6;
	QCGA r5 = -0.5*2 * (eo3 ^ ei5);
	QCGA r6 = -0.5*1*eo2 ^ ei5;
	QCGA r7 = -0.5*1*ei5 ^ eo1;

	//double phi = std::numbers::pi / 2;
	double phi = std::numbers::pi / 4;

	QCGA R1 = (cos(phi / 2) * one + sin(phi / 2) *2* r1);
	QCGA R2 = (cos(phi / 2) * one + sin(phi / 2) *2* r2);
	QCGA R3 = (cos(phi) * one + sin(phi) * (0.5 *2* r3));
	QCGA R4 = (cos(phi/2) * one + sin(phi/2) *2* r4);
	QCGA R5 = (cos(phi) * one + sin(phi) * (0.5 *2* r5));
	QCGA R6 = (cos(phi / 2) * one + sin(phi / 2) *2* r6);
	QCGA R7 = (cos(phi / 2) * one + sin(phi / 2) *2* r7);


	QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA c_euc = eo1 + e1 + 2 * e2 + 3 * e3 + 7 * ei1;
	QCGA c_35 = -2*ei3 + 3*ei5;
	QCGA c_2 = 3 * ei2;
	QCGA c_46 = 2 * ei4 + 6 * ei6;

	double y = 2;
	//QCGA CC = up(3, 2, -1);
	QCGA CC = up(2*sqrt(2), 2, sqrt(2));
	std::cout << "      C: " << C << std::endl;
	QCGA rot2 = (R1 * c_euc * ~R1) + ((R3 ^ R5) * c_35 * (~R5 ^ ~R3)) + c_2 + ((R2 ^ R4) * c_46 * (~R4 ^ ~R2)) + (0.5 * (one - cos(2 * phi) * one) * y * y * ei3) + ((0.5 * sin(2 * phi) * one) * y * y * ei5);
	std::cout << "Rotated: " << rot2 << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rot2 == CC) << std::endl;

	QCGA taylor_r = r1.rotorExponential(20, phi);
	QCGA taylor_r35 = (r3+r5).rotorExponential(20, phi);
	QCGA taylor_r46 = (r2+r4).rotorExponential(20, phi);
	QCGA rot3 = (taylor_r * c_euc * ~taylor_r) + (taylor_r35 * c_35 * ~taylor_r35) + c_2 + (taylor_r46 * c_46 * ~taylor_r46) + (0.5 * (one - cos(2 * phi) * one) * y * y * ei3) + ((0.5 * sin(2 * phi) * one) * y * y * ei5);
	std::cout << "Rotated: " << rot3 << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rot3 == CC) << std::endl;

}

void Rotor1YZ()
{
	QCGA r1 = 0.5*e2 ^ e3;
	QCGA r2 = 0.5*2*(eo6 ^ ei3);
	QCGA r6 = 0.5 * (ei6 ^ eo3);
	QCGA r7 = 0.5 * (eo2 ^ ei6);
	QCGA r3 = 0.5*2*(ei2 ^ eo6);
	QCGA r4 = 0.5*eo5 ^ ei4;
	QCGA r5 = 0.5*ei5 ^ eo4;

	//long double phi = std::numbers::pi / 4;
	long double phi = std::numbers::pi / 2;

	//QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA C = up(2, 3, -1); //eo1+2*e1+3*e2-1*e3+7*ei1-2*ei2+6*ei3+6*ei4-2*ei5-3*ei6
	//QCGA c_euc = eo1 + e1 + 2 * e2 + 3 * e3 + 7 * ei1;
	QCGA c_euc = eo1 + 2*e1 + 3 * e2 -1 * e3 + 7 * ei1;
	//QCGA c_45 = 2 * ei4 + 3 * ei5;
	QCGA c_45 = 6 * ei4 + -2 * ei5;
	//QCGA c_236 = 3*ei2-2*ei3 + 6 * ei6;
	QCGA c_236 = -2*ei2+6*ei3 -3* ei6;


	//QCGA target = up(1, 0.5*5*sqrt(2), 0.5*sqrt(2));
	QCGA target = up(2, -1, -3);
	QCGA taylor_r = r1.rotorExponential(20, phi);
	QCGA taylor_r45 = (r4+r5).rotorExponential(20, phi);
	QCGA taylor_r2367 = (r2+r3+r6+r7).rotorExponential(20, phi);
	std::cout << "      C: " << C << std::endl;
	QCGA rot1 = (taylor_r * c_euc * ~taylor_r);
	QCGA rot2 = (taylor_r45 * c_45 * ~taylor_r45);
	QCGA rot3 = (taylor_r2367 * c_236 * ~taylor_r2367);
	QCGA rot = (rot1 + rot2 + rot3)[1];
	//QCGA rot2 = (R1 * c_euc * ~R1) + ((R3 ^ R5) * c_45 * (~R5 ^ ~R3)) + c_2 + ((R2 ^ R4) * c_46 * (~R4 ^ ~R2)) + (0.5 * (one - cos(2 * phi) * one) * y * y * ei3) + ((0.5 * sin(2 * phi) * one) * y * y * ei5);
	std::cout << "rot1: " << rot1 << std::endl << std::endl;
	std::cout << "rot2: " << rot2 << std::endl << std::endl;
	std::cout << "rot3_1: " << rot3 << std::endl << std::endl;
	//std::cout << "Rotated: " << rot1+rot2+rot3 << std::endl << std::endl;
	std::cout << "Rotated: " << rot << std::endl;
	std::cout << " Target: " << target << std::endl;
	std::cout << "   Good: " << (rot == target) << std::endl << std::endl;
	//std::cout << "    c_45: " << c_45 << std::endl;
	//std::cout << "taylor45: " << taylor_r45 << std::endl;
	//std::cout << "~taylr45: " << ~taylor_r45 << std::endl;
	std::cout << "    3*ei2: " << 3*ei2 << std::endl;
	std::cout << "taylor2367: " << taylor_r2367 << std::endl;
}

void TranslatorX()
{
	QCGA t1 = -0.5 * e1 ^ ei1;
	QCGA t2 = -0.5 * e1 ^ ei2;
	QCGA t3 = -0.5 * e1 ^ ei3;
	QCGA t4 = -0.5 * e2 ^ ei4;
	QCGA t5 = -0.5 * e3 ^ ei5;

	QCGA t = t1 + t2 + t3 + t4 + t5;

	int distance = 3;
	QCGA translatorX = t.translatorExponential(20,distance);

	QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA target = up(1 + distance, 2, 3);

	QCGA translated = t * C * ~t;

	std::cout << "Translated: " << translated <<std::endl;
	std::cout << "    Target: " << target <<std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
}

void TranslatorY()
{
	QCGA t1 = -0.5 * e2 ^ ei1;
	QCGA t2 = 0.5 * e2 ^ ei2;
	QCGA t3 = -0.5 * e2 ^ ei3;
	QCGA t4 = -0.5 * e1 ^ ei4;
	QCGA t5 = -0.5 * e3 ^ ei6;

	QCGA t = t1 + t2 + t3 + t4 + t5;

	int distance = 3;
	QCGA translatorY = t.translatorExponential(20, distance);

	QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA target = up(1, 2 + distance, 3);

	QCGA translated = t * C * ~t;

	std::cout << "Translated: " << translated << std::endl;
	std::cout << "    Target: " << target << std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
}

void TranslatorZ()
{
	QCGA t1 = -0.5 * e3 ^ ei1;
	QCGA t2 = -0.5 * e3 ^ ei2;
	QCGA t3 =  0.5 * e3 ^ ei3;
	QCGA t4 = -0.5 * e1 ^ ei5;
	QCGA t5 = -0.5 * e2 ^ ei6;

	QCGA t = t1 + t2 + t3 + t4 + t5;

	int distance = 3;
	QCGA translatorZ = t.translatorExponential(20, distance);

	QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA target = up(1, 2, 3 + distance);

	QCGA translated = t * C * ~t;

	std::cout << "Translated: " << translated << std::endl;
	std::cout << "    Target: " << target << std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
}

void Scalor()
{
}

int main()
{	
	//TEST FRO REPO
	QCGA::generateGeneratingBlades();
	//Rotor3XTaylor();
	Rotor1XZ();
	//QCGA r4 = 2 * (eo4 ^ ei2);
	//QCGA r5 = 2 * (ei4 ^ eo2);
	//QCGA R4 = cos(1) * one + sin(1) * (0.5 * r4);
	//QCGA R5 = cos(1) * one + sin(1) * (0.5 * r5);
	//std::cout << (R4 ^ R5) << std::endl;
	//std::cout << (r4+r5).rotorExponential(20,1) << std::endl;
	//QCGA r2 = eo6 ^ ei5;
	//QCGA r3 = ei6 ^ eo5;
	//std::cout << ((ei5 ^ eo5) + (eo6 ^ ei6)) << std::endl;
	//std::cout << ((ei5 ^ eo5) ^ (eo6 ^ ei6)) << std::endl;
	//std::cout << com(r2, r3) << std::endl;
	//std::cout << (com(r2, r3) ^ 2) << std::endl;
	//std::cout << (com(r2, r3)^3) << std::endl;
	//std::cout << (com(r2, r3)^4) << std::endl;
	//std::cout << (com(r2, r3)^5) << std::endl;
	//std::cout << (com(r2, r3)^5) << std::endl;
	//std::cout << (com(r2, r3)^6) << std::endl;
	//std::cout << (com(r2, r3)^7) << std::endl;
	//std::cout << (com(r2, r3)^8) << std::endl;
	//std::cout << (com(r2, r3)^9) << std::endl;
	//QCGA R2 = cos(1) * one + sin(1) * r2;
	//QCGA R3 = cos(1) * one + sin(1) * r3;
	//std::cout << (R2 ^ R3) << std::endl;
	//std::cout << (r2+r3).rotorExponential(20,1) << std::endl;
	//Rotor3XY();
	return 0;
}
