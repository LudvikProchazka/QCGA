#pragma once
#include "QCGA.h"
#include "Blade.h"

#include <numbers>

QCGA com(const QCGA& a, const QCGA& b)
{
	return 0.5 * ((a * b) - (b * a));
}

void RotationExample() //Parabolas
{
	QCGA r1 = e2 ^ e3;
	QCGA r2 = eo6 ^ ei3;
	QCGA r3 = ei5 ^ eo4;
	QCGA r4 = 2 * (ei6 ^ eo3);
	QCGA r5 = 2 * (eo2 ^ ei6);
	QCGA r6 = ei2 ^ eo6;
	QCGA r7 = eo5 ^ ei4;

	double phi = std::numbers::pi / 2.0;

	QCGA Q = makeQuadric(0, 0, 0, 2.0 / 3.0, -4.0 / 3.0, -4.0 / 3.0, 0, 1, 0, 0);

	QCGA r = r1 + r2 + r3 + r4 + r5 + r6 + r7; //create generator

	QCGA rotor = r.rotorExponential(30, phi);
	Blade rotated = (rotor * Q * ~rotor)[1];
	std::cout << "Quadric: " << Q << std::endl;
	std::cout << "Rotated: " << rotated << std::endl;
}

void OPNS_IPNS_Duality()
{
	QCGA p1 = up(0.0, 0.0, 0.0);
	QCGA p2 = up(1.0, -2.0, 1.0);
	QCGA p3 = up(-1.0, -2.0, -1.0);
	QCGA p4 = up(1.0, -2.0, -1.0);
	QCGA p5 = up(-1.0, -2.0, 1.0);
	QCGA p6 = up(1.0, -1.0, 0.0);
	QCGA p7 = up(-1.0, -1.0, 0.0);
	QCGA p8 = up(0.0, -1.0, 1.0);
	QCGA p9 = up(0.0, -1.0, -1.0);

	Blade OPNS = p1 ^ p2 ^ p3 ^ p4 ^ p5 ^ p6 ^ p7 ^ p8 ^ p9 ^ eo2 ^ eo3 ^ eo4 ^ eo5 ^ eo6;
	Blade OPNS_dual = OPNS.dual();

	Blade IPNS = makeQuadric(0, 0, 0, 2.0 / 3.0, -4.0 / 3.0, -4.0 / 3.0, 0, 1, 0, 0);

	std::cout << "-1/48*OPNS: " << ((-1.0 / 48.0) * OPNS_dual) << std::endl;
	std::cout << "      IPNS: " << IPNS << std::endl;
}

void Elipsoid()
{
	Blade Q = makeQuadric(0.0, 0.0, 0.0, 5.0 / 3.0, -1.0 / 3.0, -7.0 / 3.0, 5.0, 3.0, -5.0, 5.0);

	double distance = -10.0;
	QCGA T1 = one - 0.5 * distance * (e1 ^ ei1);
	QCGA T2 = one - 0.5 * distance * (e1 ^ ei2) + 0.25 * pow(distance, 2) * (ei1 ^ ei2);
	QCGA T3 = one - 0.5 * distance * (e1 ^ ei3) + 0.25 * pow(distance, 2) * (ei1 ^ ei3) + 0.25 * pow(distance, 2) * (ei2 ^ ei3);
	QCGA T4 = one - 0.5 * distance * (e2 ^ ei4);
	QCGA T5 = one - 0.5 * distance * (e3 ^ ei5);
	QCGA T = T1 * T2 * T3 * T4 * T5; //Translator in x direction

	double phi = std::numbers::pi / 4.0;
	QCGA r1 = e1 ^ e2;
	QCGA r2 = eo6 ^ ei5;
	QCGA r3 = ei6 ^ eo5;
	QCGA r4 = 2 * (eo4 ^ ei2);
	QCGA r5 = 2 * (ei4 ^ eo2);
	QCGA r6 = eo4 ^ ei3;
	QCGA r = r1 + r2 + r3 + r4 + r5 + r6;
	QCGA R = r.rotorExponential(30, phi); //Rotor in the xy-plane

	Blade transformed = (R * (T * Q * ~T)[1] * ~R)[1];

	std::cout << "Quadric: " << Q << std::endl;
	std::cout << "transformed: " << transformed << std::endl;
}

void Elipsoid2()
{
	Blade Q = makeQuadric(0.0, 0.0, 0.0, 5.0 / 3.0, -1.0 / 3.0, -7.0 / 3.0, 5.0, 3.0, -5.0, 5.0);

	double distance1x = 5.0;
	double distance1y = 3.0 / 2;
	double distance2x = 10.0;
	double distance2y = 6.0;
	double phi = 3.0 * std::numbers::pi / 4.0;

	QCGA T1_x = one - 0.5 * distance1x * (e1 ^ ei1);
	QCGA T2_x = one - 0.5 * distance1x * (e1 ^ ei2) + 0.25 * pow(distance1x, 2) * (ei1 ^ ei2);
	QCGA T3_x = one - 0.5 * distance1x * (e1 ^ ei3) + 0.25 * pow(distance1x, 2) * (ei1 ^ ei3) + 0.25 * pow(distance1x, 2) * (ei2 ^ ei3);
	QCGA T4_x = one - 0.5 * distance1x * (e2 ^ ei4);
	QCGA T5_x = one - 0.5 * distance1x * (e3 ^ ei5);
	QCGA T_x = T1_x * T2_x * T3_x * T4_x * T5_x; //Translator in x direction into origin

	QCGA T1_y = one - 0.5 * distance1y * (e2 ^ ei1);
	QCGA T2_y = one + 0.5 * distance1y * (e2 ^ ei2) - 0.25 * pow(distance1y, 2) * (ei1 ^ ei2);
	QCGA T3_y = one - 0.5 * distance1y * (e1 ^ ei4);
	QCGA T4_y = one - 0.5 * distance1y * (e3 ^ ei6);
	QCGA T_y = T1_y * T2_y * T3_y * T4_y; //Translator in y direction into origin

	QCGA r1 = e1 ^ e2;
	QCGA r2 = eo6 ^ ei5;
	QCGA r3 = ei6 ^ eo5;
	QCGA r4 = 2 * (eo4 ^ ei2);
	QCGA r5 = 2 * (ei4 ^ eo2);
	QCGA r6 = eo4 ^ ei3;
	QCGA r = r1 + r2 + r3 + r4 + r5 + r6;
	QCGA R = r.rotorExponential(40, phi); //Rotor in the xy-plane

	Blade rotatedInOrigin = (R * (T_y * (T_x * Q * ~T_x)[1] * ~T_y)[1] * ~R)[1];

	T1_x = one - 0.5 * distance2x * (e1 ^ ei1);
	T2_x = one - 0.5 * distance2x * (e1 ^ ei2) + 0.25 * pow(distance2x, 2) * (ei1 ^ ei2);
	T3_x = one - 0.5 * distance2x * (e1 ^ ei3) + 0.25 * pow(distance2x, 2) * (ei1 ^ ei3) + 0.25 * pow(distance2x, 2) * (ei2 ^ ei3);
	T4_x = one - 0.5 * distance2x * (e2 ^ ei4);
	T5_x = one - 0.5 * distance2x * (e3 ^ ei5);
	T_x = T1_x * T2_x * T3_x * T4_x * T5_x;

	T1_y = one - 0.5 * distance2y * (e2 ^ ei1);
	T2_y = one + 0.5 * distance2y * (e2 ^ ei2) - 0.25 * pow(distance2y, 2) * (ei1 ^ ei2);
	T3_y = one - 0.5 * distance2y * (e1 ^ ei4);
	T4_y = one - 0.5 * distance2y * (e3 ^ ei6);
	T_y = T1_y * T2_y * T3_y * T4_y;

	Blade transformed = (T_y * (T_x * rotatedInOrigin * ~T_x)[1] * ~T_y)[1];

	std::cout << "Quadric: " << Q << std::endl;
	std::cout << "transformed: " << transformed << std::endl;
}

void RotorXY()
{
	QCGA r1 = e1 ^ e2;
	QCGA r2 = eo6 ^ ei5;
	QCGA r3 = ei6 ^ eo5;
	QCGA r4 = 2 * (eo4 ^ ei2);
	QCGA r5 = 2 * (ei4 ^ eo2);
	QCGA r6 = eo4 ^ ei3;

	QCGA r = r1 + r2 + r3 + r4 + r5 + r6;

	double phi = std::numbers::pi / 4.0;
	QCGA R1 = cos(phi / 2) * one + sin(phi / 2) * r1;
	QCGA R2 = cos(phi / 2) * one + sin(phi / 2) * r2;
	QCGA R3 = cos(phi / 2) * one + sin(phi / 2) * r3;
	QCGA R4 = cos(phi) * one + sin(phi) * (0.5 * r4);
	QCGA R5 = cos(phi) * one + sin(phi) * (0.5 * r5);


	QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1-1.5*ei2-4*ei3+2*ei4+3*ei5+6*ei6
	QCGA c_cga = eo1 + e1 + 2 * e2 + 3 * e3 + 7 * ei1;
	QCGA c_24 = -1.5 * ei2 + 2 * ei4;
	QCGA c_56 = 3 * ei5 + 6 * ei6;
	QCGA c_3 = -4 * ei3;

	double x = 1;
	double y = 2;
	double z = 3;
	double theta = atan(2) - phi;
	QCGA CC = up(sqrt(5) * cos(theta), sqrt(5) * sin(theta), 3);
	QCGA _rotated = (R1 * c_cga * ~R1) + ((R2 ^ R3) * c_56 * (~R3 ^ ~R2)) + c_3 + ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + (-0.5 * sin(phi) * sin(phi) * (x * x - y * y) + sin(phi) * cos(phi) * x * y) * ei3;
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
	const QCGA r1 = -1 * e3 ^ e1; //rotating in opposite direction, thus -1*, viz corresponding matrices
	const QCGA r2 = -1 * ei4 ^ eo6;
	const QCGA r3 = -2 * (ei3 ^ eo5);
	const QCGA r4 = -1 * (ei2 ^ eo5);
	const QCGA r5 = -1 * eo4 ^ ei6;
	const QCGA r6 = -2 * (eo3 ^ ei5);

	const QCGA r = r1 + r2 + r3 + r4 + r5 + r6;
	//double phi = std::numbers::pi / 2;
	const double phi = std::numbers::pi / 4.0;


	//QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	const QCGA C = up(1, 3, 2); //eo1+e1+3*e2+2*e3+7*ei1-2*ei2-3*ei3+3*ei4+2*ei5+6*ei6
	const QCGA c_euc = eo1 + e1 + 3 * e2 + 2 * e3 + 7 * ei1;
	const QCGA c_35 = -3 * ei3 + 2 * ei5;
	const QCGA c_2 = -2 * ei2;
	const QCGA c_46 = 3 * ei4 + 6 * ei6;

	//double y = 2;
	const double y = 3;
	//QCGA CC = up(2*sqrt(2), 2, sqrt(2));
	const double theta = atan(2) - phi;
	const QCGA CC = up(sqrt(5) * cos(theta), 3, sqrt(5) * sin(theta));
	const QCGA rotor = r.rotorExponential(20, phi);
	const QCGA rotated = (rotor * C * ~rotor)[1];
	std::cout << "Rotated: " << rotated << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rotated == CC) << std::endl;
}
void RotorYZ()
{
	const QCGA r1 = e2 ^ e3;
	const QCGA r2 = eo6 ^ ei3;
	const QCGA r3 = ei5 ^ eo4;
	const QCGA r4 = 2 * (ei6 ^ eo3);
	const QCGA r5 = 2 * (eo2 ^ ei6);
	const QCGA r6 = ei2 ^ eo6;
	const QCGA r7 = eo5 ^ ei4;

	const QCGA r = r1 + r2 + r3 + r4 + r5 + r6 + r7;
	const long double phi = std::numbers::pi / 4;
	//long double phi = std::numbers::pi / 2;

	const QCGA R1 = cos(phi / 2) * one + sin(phi / 2) * r1;
	const QCGA R2 = cos(phi / 2) * one + sin(phi / 2) * r2;
	const QCGA R3 = cos(phi / 2) * one + sin(phi / 2) * r3;
	const QCGA R4 = cos(phi / 2) * one + sin(phi / 2) * r4;
	const QCGA R5 = cos(phi / 2) * one + sin(phi / 2) * r5;
	const QCGA R6 = cos(phi / 2) * one + sin(phi / 2) * r6;
	const QCGA R7 = cos(phi / 2) * one + sin(phi / 2) * r7;

	const QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1-1.5*ei2-4*ei3+2*ei4+3*ei5+6*ei6
	const QCGA c_cga = eo1 + e1 + 2 * e2 + 3 * e3 + 7 * ei1;
	const QCGA c_45 = 2 * ei4 + 3 * ei5;
	const QCGA c_236 = -1.5 * ei2 - 4 * ei3 + 6 * ei6;
	//QCGA C = MujUp(2, 3, -1); 


	const QCGA target = up(1, 0.5 * 5 * sqrt(2), 0.5 * sqrt(2));
	//QCGA target = MujUp(2, -1, -3);
	const QCGA rot_236 = (r2 + r4 + r5 + r6).rotorExponential(20, phi);

	const QCGA _rotated = ((R1 * c_cga * ~R1) + ((R3 ^ R7) * c_45 * (~R7 ^ ~R3)) + (rot_236 * c_236 * ~rot_236))[1];
	std::cout << "Rotated: " << _rotated << std::endl;
	std::cout << " Target: " << target << std::endl;
	std::cout << "   Good: " << (_rotated == target) << std::endl;


	const QCGA rotor = r.rotorExponential(20, phi);
	const QCGA rotated = (rotor * C * ~rotor)[1];
	std::cout << "Rotated: " << rotated << std::endl;
	std::cout << " Target: " << target << std::endl;
	std::cout << "   Good: " << (rotated == target) << std::endl;
}
void TranslatorX()
{
	const QCGA t1 = -1 * e1 ^ ei1;
	const QCGA t2 = -1 * e1 ^ ei2;
	const QCGA t3 = -1 * e1 ^ ei3;
	const QCGA t4 = -1 * e2 ^ ei4;
	const QCGA t5 = -1 * e3 ^ ei5;

	const QCGA t = t1 + t2 + t3 + t4 + t5;
	const double distance = -1 / double(3);

	const QCGA T1 = one - 0.5 * distance * (e1 ^ ei1);
	const QCGA T2 = one - 0.5 * distance * (e1 ^ ei2) + 0.25 * pow(distance, 2) * (ei1 ^ ei2);
	const QCGA T3 = one - 0.5 * distance * (e1 ^ ei3) + 0.25 * pow(distance, 2) * (ei1 ^ ei3) + 0.25 * pow(distance, 2) * (ei2 ^ ei3);
	const QCGA T4 = one - 0.5 * distance * (e2 ^ ei4);
	const QCGA T5 = one - 0.5 * distance * (e3 ^ ei5);

	const QCGA translatorX = t.translatorExponential(20, distance);
	const QCGA translatorX2 = T1 * T2 * T3 * T4 * T5;

	const QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	const QCGA target = up(1 + distance, 2, 3);

	const QCGA translated = translatorX * C * ~translatorX;
	const QCGA translated2 = translatorX2 * C * ~translatorX2;

	std::cout << "Translated: " << translated << std::endl;
	std::cout << "Translated: " << translated2 << std::endl;
	std::cout << "    Target: " << target << std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
	std::cout << "      Good: " << (translated2 == target) << std::endl;
}
void TranslatorY()
{
	const QCGA t1 = -1 * e2 ^ ei1;
	const QCGA t2 = 1 * e2 ^ ei2;
	const QCGA t3 = -1 * e1 ^ ei4;
	const QCGA t4 = -1 * e3 ^ ei6;

	const QCGA t = t1 + t2 + t3 + t4;
	const int distance = 7;

	const QCGA T1 = one - 0.5 * distance * (e2 ^ ei1);
	const QCGA T2 = one + 0.5 * distance * (e2 ^ ei2) - 0.25 * pow(distance, 2) * (ei1 ^ ei2);
	const QCGA T3 = one - 0.5 * distance * (e1 ^ ei4);
	const QCGA T4 = one - 0.5 * distance * (e3 ^ ei6);

	const QCGA translatorY = t.translatorExponential(20, distance);
	const QCGA translatorY2 = T1 * T2 * T3 * T4;

	const QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	const QCGA target = up(1, 2 + distance, 3);

	const QCGA translated = translatorY * C * ~translatorY;
	const QCGA translated2 = translatorY2 * C * ~translatorY2;

	std::cout << "Translated: " << translated << std::endl;
	std::cout << "Translated: " << translated2 << std::endl;
	std::cout << "    Target: " << target << std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
	std::cout << "      Good: " << (translated2 == target) << std::endl;
}
void TranslatorZ()
{
	const QCGA t1 = -1 * e3 ^ ei1;
	const QCGA t2 = 1 * e3 ^ ei3;
	const QCGA t3 = -1 * e1 ^ ei5;
	const QCGA t4 = -1 * e2 ^ ei6;

	const QCGA t = t1 + t2 + t3 + t4;
	const int distance = -4;

	const QCGA T1 = one - 0.5 * distance * (e3 ^ ei1);
	const QCGA T2 = one + 0.5 * distance * (e3 ^ ei3) - 0.25 * pow(distance, 2) * (ei1 ^ ei3);
	const QCGA T3 = one - 0.5 * distance * (e1 ^ ei5);
	const QCGA T4 = one - 0.5 * distance * (e2 ^ ei6);

	const QCGA translatorZ = t.translatorExponential(20, distance);
	const QCGA translatorZ2 = T1 * T2 * T3 * T4;

	const QCGA C = up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	const QCGA target = up(1, 2, 3 + distance);

	const QCGA translated = translatorZ * C * ~translatorZ;
	const QCGA translated2 = translatorZ2 * C * ~translatorZ2;

	std::cout << "Translated: " << translated << std::endl;
	std::cout << "Translated: " << translated2 << std::endl;
	std::cout << "    Target: " << target << std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
	std::cout << "      Good: " << (translated2 == target) << std::endl;
}