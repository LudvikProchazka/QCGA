#include "QCGA.h"
#include "Blade.h"

#include <iostream>
#include <numbers>
#include <iomanip>      // std::setprecision
#include <chrono>

enum rotation_planes
{
	xy = 1,
	xz,
	yz
};
enum translation_directions
{
	x = 1,
	y,
	z
};

QCGA com(const QCGA& a, const QCGA& b)
{
	return 0.5*((a * b) - (b * a));
}

void RotorXY()
{
	QCGA r1 = e1 ^ e2;
	QCGA r2 = eo6 ^ ei5;
	QCGA r3 = ei6 ^ eo5;
	QCGA r4 = 2 * (eo4 ^ ei2);
	QCGA r5 = 2 * (ei4 ^ eo2);
	QCGA r6 = ei4 ^ eo3;
	QCGA r7 = eo1 ^ ei2;

	QCGA r = r1 + r2 + r3 + r4 + r5 + r6 + r7;

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
	QCGA rot1 = (R1 * c_euc * ~R1) + ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + c_3 + ((R2 ^ R3) * c_56 * (~R3 ^ ~R2)) + (0.5 * (one - cos(2 * phi) * one) * z * z * ei2) + ((0.5 * sin(2 * phi) * one) * z * z * ei4);
	std::cout << "Rotated: " << rot1 << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rot1 == CC) << std::endl << std::endl;
	QCGA taylor_r1 = r1.rotorExponential(20, phi);
	QCGA taylor_r45 = (r4+r5).rotorExponential(20, phi);
	QCGA taylor_23 = (r2+r3).rotorExponential(20, phi);
	QCGA rot2 = (taylor_r1 * c_euc * ~taylor_r1) + (taylor_23 * c_56 * ~taylor_23) + (taylor_r45 * (c_24) * ~taylor_r45) + c_3 + (0.5 * (one - cos(2 * phi) * one) * z * z * ei2) + ((0.5 * sin(2 * phi) * one) * z * z * ei4);
	std::cout << "Rotated: " << rot2 << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rot2 == CC) << std::endl << std::endl;
	//std::cout << "       : " << ((0.5 * sin(2 * phi) * one) * z * z * ei4) << std::endl << std::endl;

	QCGA rotor = r.rotorExponential(20, phi);
	QCGA rot3 = (rotor * C * ~rotor)[1];
	std::cout << "Rotated: " << rot3 << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rot3==CC) << std::endl;
}
void RotorXYMuj()
{
	QCGA r1 = e1 ^ e2;
	QCGA r2 = eo6 ^ ei5;
	QCGA r3 = ei6 ^ eo5;
	QCGA r4 = 2 * (eo4 ^ ei2);
	QCGA r5 = 2 * (ei4 ^ eo2);
	QCGA r6 = eo4 ^ ei3;

	double phi = std::numbers::pi / 4;

	QCGA R1 = cos(phi / 2) * one + sin(phi / 2) * r1;
	QCGA R2 = cos(phi / 2) * one + sin(phi / 2) * r2;
	QCGA R3 = cos(phi / 2) * one + sin(phi / 2) * r3;
	QCGA R4 = cos(phi) * one + sin(phi) * (0.5 * r4);
	QCGA R5 = cos(phi) * one + sin(phi) * (0.5 * r5);

	QCGA r = r1 + r2 + r3 + r4 + r5 + r6;

	QCGA C = MujUp(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1-1.5*ei2-4*ei3+2*ei4+3*ei5+6*ei6
	QCGA c_cga = eo1 + e1 + 2 * e2 + 3 * e3 + 7 * ei1;
	QCGA c_24 = -1.5 * ei2 + 2 * ei4;
	QCGA c_56 = 3 * ei5 + 6 * ei6;
	QCGA c_3 = -4 * ei3;

	double x = 1;
	double y = 2;
	double z = 3;
	double theta = atan(2) - phi;
	QCGA CC = MujUp(sqrt(5) * cos(theta), sqrt(5) * sin(theta), 3);
	QCGA _rotated = (R1 * c_cga * ~R1) + ((R2 ^ R3) * c_56 * (~R3 ^ ~R2)) +c_3+ ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + (-0.5 * sin(phi) * sin(phi) * (x * x - y * y) + sin(phi) * cos(phi) * x * y) * ei3;
	std::cout << "Rotated: " << _rotated << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (_rotated == CC) << std::endl;

	QCGA rotor = r.rotorExponential(20, phi);
	QCGA rotated = (rotor * C * ~rotor)[1];
	std::cout << "Rotated: " << rotated << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rotated == CC) << std::endl;
}

void RotorXZ()
{
	QCGA r1 = -1*e3 ^ e1; //rotating in opposite direction, thus -1*, viz corresponding matrices
	QCGA r2 = -1*ei4 ^ eo6;
	QCGA r3 = -2 * (ei3 ^ eo5);
	QCGA r4 = -1*eo4 ^ ei6;
	QCGA r5 = -2 * (eo3 ^ ei5);
	QCGA r6 = -1*eo2 ^ ei5;
	QCGA r7 = -1*ei5 ^ eo1;

	QCGA r = r1 + r2 + r3 + r4 + r5 + r6 + r7;
	//double phi = std::numbers::pi / 2;
	double phi = std::numbers::pi / 4;

	QCGA R1 = (cos(phi / 2) * one + sin(phi / 2) * r1);
	QCGA R2 = (cos(phi / 2) * one + sin(phi / 2) * r2);
	QCGA R3 = (cos(phi) * one + sin(phi) * (0.5 * r3));
	QCGA R4 = (cos(phi/2) * one + sin(phi/2) * r4);
	QCGA R5 = (cos(phi) * one + sin(phi) * (0.5 * r5));
	QCGA R6 = (cos(phi / 2) * one + sin(phi / 2) * r6);
	QCGA R7 = (cos(phi / 2) * one + sin(phi / 2) * r7);


	//QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA C = up(1, 3, 2); //eo1+e1+3*e2+2*e3+7*ei1-2*ei2-3*ei3+3*ei4+2*ei5+6*ei6
	QCGA c_euc = eo1 + e1 + 3 * e2 + 2 * e3 + 7 * ei1;
	QCGA c_35 = -3*ei3 + 2*ei5;
	QCGA c_2 = -2 * ei2;
	QCGA c_46 = 3 * ei4 + 6 * ei6;

	//double y = 2;
	double y = 3;
	//QCGA CC = up(2*sqrt(2), 2, sqrt(2));
	double theta = atan(2) - phi;
	QCGA CC = up(sqrt(5)*cos(theta), 3, sqrt(5) * sin(theta));
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


	QCGA rotor = r.rotorExponential(20, phi);
	QCGA rot4 = (rotor * C * ~rotor)[1];
	std::cout << "Rotated: " << rot4 << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rot4==CC) << std::endl;
}
void RotorXZMuj()
{
	QCGA r1 = -1 * e3 ^ e1; //rotating in opposite direction, thus -1*, viz corresponding matrices
	QCGA r2 = -1 * ei4 ^ eo6;
	QCGA r3 = -2 * (ei3 ^ eo5);
	QCGA r4 = -1 * (ei2 ^ eo5);
	QCGA r5 = -1 * eo4 ^ ei6;
	QCGA r6 = -2 * (eo3 ^ ei5);

	QCGA r = r1 + r2 + r3 + r4 + r5 + r6;
	//double phi = std::numbers::pi / 2;
	double phi = std::numbers::pi / 4;


	//QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA C = MujUp(1, 3, 2); //eo1+e1+3*e2+2*e3+7*ei1-2*ei2-3*ei3+3*ei4+2*ei5+6*ei6
	QCGA c_euc = eo1 + e1 + 3 * e2 + 2 * e3 + 7 * ei1;
	QCGA c_35 = -3 * ei3 + 2 * ei5;
	QCGA c_2 = -2 * ei2;
	QCGA c_46 = 3 * ei4 + 6 * ei6;

	//double y = 2;
	double y = 3;
	//QCGA CC = up(2*sqrt(2), 2, sqrt(2));
	double theta = atan(2) - phi;
	QCGA CC = MujUp(sqrt(5) * cos(theta), 3, sqrt(5) * sin(theta));
	QCGA rotor = r.rotorExponential(20, phi);
	QCGA rotated = (rotor * C * ~rotor)[1];
	std::cout << "Rotated: " << rotated << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rotated == CC) << std::endl;
}

void RotorYZ()
{
	QCGA r1 = e2 ^ e3;
	QCGA r2 = 2*(eo6 ^ ei3);
	QCGA r6 = (ei6 ^ eo3);
	QCGA r7 = (eo2 ^ ei6);
	QCGA r3 = 2*(ei2 ^ eo6);
	QCGA r4 = eo5 ^ ei4;
	QCGA r5 = ei5 ^ eo4;

	QCGA r = r1 + r2 + r3 + r4 + r5 + r6 + r7;
	long double phi = std::numbers::pi / 4;
	//long double phi = std::numbers::pi / 2;

	QCGA R4 = cos(phi / 2)* one + sin(phi / 2) * r4;
	QCGA R5 = cos(phi / 2)* one + sin(phi / 2) * r5;

	QCGA R2 = cos(phi / 2)* one + sin(phi / 2) * r2;
	QCGA R3 = cos(phi / 2)* one + sin(phi / 2) * r3;
	QCGA R6 = cos(phi / 2)* one + sin(phi / 2) * r6;
	QCGA R7 = cos(phi / 2)* one + sin(phi / 2) * r7;

	QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	//QCGA C = up(2, 3, -1); //eo1+2*e1+3*e2-1*e3+7*ei1-2*ei2+6*ei3+6*ei4-2*ei5-3*ei6
	QCGA c_euc = eo1 + e1 + 2 * e2 + 3 * e3 + 7 * ei1;
	//QCGA c_euc = eo1 + 2*e1 + 3 * e2 -1 * e3 + 7 * ei1;
	QCGA c_45 = 2 * ei4 + 3 * ei5;
	//QCGA c_45 = 6 * ei4 + -2 * ei5;
	QCGA c_236 = 3*ei2-2*ei3 + 6 * ei6;
	//QCGA c_236 = -2*ei2+6*ei3 -3* ei6;


	QCGA target = up(1, 0.5*5*sqrt(2), 0.5*sqrt(2));
	//QCGA target = up(2, -1, -3);
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
	std::cout << "rot3: " << rot3 << std::endl << std::endl;
	//std::cout << "Rotated: " << rot1+rot2+rot3 << std::endl << std::endl;
	std::cout << "Rotated: " << rot << std::endl;
	std::cout << " Target: " << target << std::endl;
	std::cout << "   Good: " << (rot == target) << std::endl << std::endl;

	QCGA rotor = r.rotorExponential(20, phi);
	QCGA rot4 = (rotor * C * ~rotor)[1];
	std::cout << "Rotated: " << rot4 << std::endl;
	std::cout << " Target: " << target << std::endl;
	std::cout << "   Good: " << (rot4 == target) << std::endl;
}
void RotorYZMuj()
{
	QCGA r1 = e2 ^ e3;
	QCGA r2 = eo6 ^ ei3;
	QCGA r3 = ei5 ^ eo4;
	QCGA r4 = 2 * (ei6 ^ eo3);
	QCGA r5 = 2 * (eo2 ^ ei6);
	QCGA r6 = ei2 ^ eo6;
	QCGA r7 = eo5 ^ ei4;

	QCGA r = r1 + r2 + r3 + r4 + r5 + r6 + r7;
	long double phi = std::numbers::pi / 4;
	//long double phi = std::numbers::pi / 2;

	QCGA R1 = cos(phi / 2) * one + sin(phi / 2) * r1;
	QCGA R2 = cos(phi / 2) * one + sin(phi / 2) * r2;
	QCGA R3 = cos(phi / 2) * one + sin(phi / 2) * r3;
	QCGA R4 = cos(phi / 2) * one + sin(phi / 2) * r4;
	QCGA R5 = cos(phi / 2) * one + sin(phi / 2) * r5;
	QCGA R6 = cos(phi / 2) * one + sin(phi / 2) * r6;
	QCGA R7 = cos(phi / 2) * one + sin(phi / 2) * r7;

	QCGA C = MujUp(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1-1.5*ei2-4*ei3+2*ei4+3*ei5+6*ei6
	QCGA c_cga = eo1 + e1 + 2 * e2 + 3 * e3 + 7 * ei1;
	QCGA c_45 = 2 * ei4 + 3 * ei5;
	QCGA c_236 = -1.5 * ei2 - 4 * ei3 + 6 * ei6;
	//QCGA C = MujUp(2, 3, -1); 


	QCGA target = MujUp(1, 0.5 * 5 * sqrt(2), 0.5 * sqrt(2));
	//QCGA target = MujUp(2, -1, -3);
	QCGA rot_236 = (r2 + r4 + r5 + r6).rotorExponential(20, phi);

	QCGA _rotated = ((R1 * c_cga * ~R1) + ((R3 ^ R7) * c_45 * (~R7 ^ ~R3)) + (rot_236 * c_236 * ~rot_236))[1];
	std::cout << "Rotated: " << _rotated << std::endl;
	std::cout << " Target: " << target << std::endl;
	std::cout << "   Good: " << (_rotated == target) << std::endl;


	QCGA rotor = r.rotorExponential(20, phi);
	QCGA rotated = (rotor * C * ~rotor)[1];
	std::cout << "Rotated: " << rotated << std::endl;
	std::cout << " Target: " << target << std::endl;
	std::cout << "   Good: " << (rotated == target) << std::endl;
}

void TranslatorX()
{
	QCGA t1 = -1*e1 ^ ei1;
	QCGA t2 = -1*e1 ^ ei2;
	QCGA t3 = -1*e1 ^ ei3;
	QCGA t4 = -1*e2 ^ ei4;
	QCGA t5 = -1*e3 ^ ei5;

	QCGA t = t1 + t2 + t3 + t4 + t5;

	int distance = 3;
	QCGA translatorX = t.translatorExponential(20,distance);

	QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA target = up(1 + distance, 2, 3);

	QCGA translated = translatorX * C * ~translatorX;

	std::cout << "Translated: " << translated <<std::endl;
	std::cout << "    Target: " << target <<std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
}
void TranslatorXMuj()
{
	QCGA t1 = -1 * e1 ^ ei1;
	QCGA t2 = -1 * e1 ^ ei2;
	QCGA t3 = -1 * e1 ^ ei3;
	QCGA t4 = -1 * e2 ^ ei4;
	QCGA t5 = -1 * e3 ^ ei5;

	QCGA t = t1 + t2 + t3 + t4 + t5;
	int distance = -1/double(3);

	QCGA T1 = one - 0.5 * distance * (e1 ^ ei1);
	QCGA T2 = one - 0.5 * distance * (e1 ^ ei2) + 0.25 * pow(distance, 2) * (ei1 ^ ei2);
	QCGA T3 = one - 0.5 * distance * (e1 ^ ei3) + 0.25 * pow(distance, 2) * (ei1 ^ ei3) + 0.25 * pow(distance, 2) * (ei2 ^ ei3);
	QCGA T4 = one - 0.5 * distance * (e2 ^ ei4);
	QCGA T5 = one - 0.5 * distance * (e3 ^ ei5);

	QCGA translatorX = t.translatorExponential(20, distance);
	QCGA translatorX2 = T1 * T2	* T3 * T4 * T5;

	QCGA C = MujUp(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA target = MujUp(1 + distance, 2, 3);

	QCGA translated = translatorX * C * ~translatorX;
	QCGA translated2 = translatorX2 * C * ~translatorX2;

	std::cout << "Translated: " << translated << std::endl;
	std::cout << "Translated: " << translated2 << std::endl;
	std::cout << "    Target: " << target << std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
	std::cout << "      Good: " << (translated2 == target) << std::endl;
}

void TranslatorY()
{
	QCGA t1 = -1* e2 ^ ei1;
	QCGA t2 =  1* e2 ^ ei2;
	QCGA t3 = -1* e2 ^ ei3;
	QCGA t4 = -1* e1 ^ ei4;
	QCGA t5 = -1* e3 ^ ei6;

	QCGA t = t1 + t2 + t3 + t4 + t5;

	int distance = 3;
	QCGA translatorY = t.translatorExponential(20, distance);

	QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA target = up(1, 2 + distance, 3);

	QCGA translated = translatorY * C * ~translatorY;

	std::cout << "Translated: " << translated << std::endl;
	std::cout << "    Target: " << target << std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
}
void TranslatorYMuj()
{
	QCGA t1 = -1 * e2 ^ ei1;
	QCGA t2 = 1 * e2 ^ ei2;
	QCGA t3 = -1 * e1 ^ ei4;
	QCGA t4 = -1 * e3 ^ ei6;

	QCGA t = t1 + t2 + t3 + t4;
	int distance = 7;

	QCGA T1 = one - 0.5 * distance * (e2 ^ ei1);
	QCGA T2 = one + 0.5 * distance * (e2 ^ ei2) - 0.25 * pow(distance, 2) * (ei1 ^ ei2);
	QCGA T3 = one - 0.5 * distance * (e1 ^ ei4);
	QCGA T4 = one - 0.5 * distance * (e3 ^ ei6);

	QCGA translatorY = t.translatorExponential(20, distance);
	QCGA translatorY2 = T1*T2*T3*T4;

	QCGA C = MujUp(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA target = MujUp(1, 2 + distance, 3);

	QCGA translated = translatorY * C * ~translatorY;
	QCGA translated2 = translatorY2 * C * ~translatorY2;

	std::cout << "Translated: " << translated << std::endl;
	std::cout << "Translated: " << translated2 << std::endl;
	std::cout << "    Target: " << target << std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
	std::cout << "      Good: " << (translated2 == target) << std::endl;
}

void TranslatorZ()
{
	QCGA t1 = -1 * e3 ^ ei1;
	QCGA t2 = -1 * e3 ^ ei2;
	QCGA t3 =  1 * e3 ^ ei3;
	QCGA t4 = -1 * e1 ^ ei5;
	QCGA t5 = -1 * e2 ^ ei6;

	QCGA t = t1 + t2 + t3 + t4 + t5;

	int distance = 3;
	QCGA translatorZ = t.translatorExponential(20, distance);

	QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA target = up(1, 2, 3 + distance);

	QCGA translated = translatorZ * C * ~translatorZ;

	std::cout << "Translated: " << translated << std::endl;
	std::cout << "    Target: " << target << std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
}
void TranslatorZMuj()
{
	QCGA t1 = -1 * e3 ^ ei1;
	QCGA t2 = 1 * e3 ^ ei3;
	QCGA t3 = -1 * e1 ^ ei5;
	QCGA t4 = -1 * e2 ^ ei6;

	QCGA t = t1 + t2 + t3 + t4;
	int distance = -4;

	QCGA T1 = one - 0.5 * distance * (e3 ^ ei1);
	QCGA T2 = one + 0.5 * distance * (e3 ^ ei3) - 0.25 * pow(distance, 2) * (ei1 ^ ei3);
	QCGA T3 = one - 0.5 * distance * (e1 ^ ei5);
	QCGA T4 = one - 0.5 * distance * (e2 ^ ei6);

	QCGA translatorZ = t.translatorExponential(20, distance);
	QCGA translatorZ2 = T1*T2*T3*T4;

	QCGA C = MujUp(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	QCGA target = MujUp(1, 2, 3 + distance);

	QCGA translated = translatorZ * C * ~translatorZ;
	QCGA translated2 = translatorZ2 * C * ~translatorZ2;

	std::cout << "Translated: " << translated << std::endl;
	std::cout << "Translated: " << translated2 << std::endl;
	std::cout << "    Target: " << target << std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
	std::cout << "      Good: " << (translated2 == target) << std::endl;
}

static int s_allocation_count = 0; //for debugging purposes

//void* operator new(size_t size) //for debugging purposes
//{
//	s_allocation_count++;
//	return malloc(size);
//}

int main()
{
	auto start = std::chrono::high_resolution_clock::now();
	QCGA::generateGeneratingBlades(); //Generates generating basis, 1,e1,e2,e3,e4,...,e15

	RotorXYMuj();

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - start;
	std::cout << "Execution time: " << elapsed.count() << " ms" << std::endl;
	return 0;
}
