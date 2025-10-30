#pragma once
#include "GAQ.h"
#include "Blade.h"
#include <numbers>

void RotationExample() //Parabolas
{
	GAQ r1 = e2 ^ e3;
	GAQ r2 = eo6 ^ ei3;
	GAQ r3 = ei5 ^ eo4;
	GAQ r4 = 2 * (ei6 ^ eo3);
	GAQ r5 = 2 * (eo2 ^ ei6);
	GAQ r6 = ei2 ^ eo6;
	GAQ r7 = eo5 ^ ei4;

	double phi = std::numbers::pi / 2.0;

	GAQ Q = MakeQuadric(0, 0, 0, 2.0 / 3.0, -4.0 / 3.0, -4.0 / 3.0, 0, 1, 0, 0);

	GAQ r = r1 + r2 + r3 + r4 + r5 + r6 + r7; //create generator

	GAQ rotor = r.RotorExponential(30, phi);
	Blade rotated = (rotor * Q * ~rotor)[1];
	std::cout << "Quadric: " << Q << std::endl;
	std::cout << "Rotated: " << rotated << std::endl;
}

void OPNS_IPNS_Duality()
{
	GAQ p1 = Up(0.0, 0.0, 0.0);
	GAQ p2 = Up(1.0, -2.0, 1.0);
	GAQ p3 = Up(-1.0, -2.0, -1.0);
	GAQ p4 = Up(1.0, -2.0, -1.0);
	GAQ p5 = Up(-1.0, -2.0, 1.0);
	GAQ p6 = Up(1.0, -1.0, 0.0);
	GAQ p7 = Up(-1.0, -1.0, 0.0);
	GAQ p8 = Up(0.0, -1.0, 1.0);
	GAQ p9 = Up(0.0, -1.0, -1.0);

	Blade OPNS = p1 ^ p2 ^ p3 ^ p4 ^ p5 ^ p6 ^ p7 ^ p8 ^ p9 ^ eo2 ^ eo3 ^ eo4 ^ eo5 ^ eo6;
	Blade OPNS_dual = OPNS.Dual();

	Blade IPNS = MakeQuadric(0, 0, 0, 2.0 / 3.0, -4.0 / 3.0, -4.0 / 3.0, 0, 1, 0, 0);

	std::cout << "-1/48*OPNS: " << ((-1.0 / 48.0) * OPNS_dual) << std::endl;
	std::cout << "      IPNS: " << IPNS << std::endl;
}

void Elipsoid()
{
	Blade Q = MakeQuadric(0.0, 0.0, 0.0, 5.0 / 3.0, -1.0 / 3.0, -7.0 / 3.0, 5.0, 3.0, -5.0, 5.0);

	double distance = -10.0;
	GAQ T1 = one - 0.5 * distance * (e1 ^ ei1);
	GAQ T2 = one - 0.5 * distance * (e1 ^ ei2) + 0.25 * pow(distance, 2) * (ei1 ^ ei2);
	GAQ T3 = one - 0.5 * distance * (e1 ^ ei3) + 0.25 * pow(distance, 2) * (ei1 ^ ei3) + 0.25 * pow(distance, 2) * (ei2 ^ ei3);
	GAQ T4 = one - 0.5 * distance * (e2 ^ ei4);
	GAQ T5 = one - 0.5 * distance * (e3 ^ ei5);
	GAQ T = T1 * T2 * T3 * T4 * T5; //Translator in x direction

	double phi = std::numbers::pi / 4.0;
	GAQ r1 = e1 ^ e2;
	GAQ r2 = eo6 ^ ei5;
	GAQ r3 = ei6 ^ eo5;
	GAQ r4 = 2 * (eo4 ^ ei2);
	GAQ r5 = 2 * (ei4 ^ eo2);
	GAQ r6 = eo4 ^ ei3;
	GAQ r = r1 + r2 + r3 + r4 + r5 + r6;
	GAQ R = r.RotorExponential(30, phi); //Rotor in the xy-plane

	Blade transformed = (R * (T * Q * ~T)[1] * ~R)[1];

	std::cout << "Quadric: " << Q << std::endl;
	std::cout << "transformed: " << transformed << std::endl;
}

void Elipsoid2()
{
	Blade Q = MakeQuadric(0.0, 0.0, 0.0, 5.0 / 3.0, -1.0 / 3.0, -7.0 / 3.0, 5.0, 3.0, -5.0, 5.0);

	double distance1x = 5.0;
	double distance1y = 3.0 / 2;
	double distance2x = 10.0;
	double distance2y = 6.0;
	double phi = 3.0 * std::numbers::pi / 4.0;

	GAQ T1_x = one - 0.5 * distance1x * (e1 ^ ei1);
	GAQ T2_x = one - 0.5 * distance1x * (e1 ^ ei2) + 0.25 * pow(distance1x, 2) * (ei1 ^ ei2);
	GAQ T3_x = one - 0.5 * distance1x * (e1 ^ ei3) + 0.25 * pow(distance1x, 2) * (ei1 ^ ei3) + 0.25 * pow(distance1x, 2) * (ei2 ^ ei3);
	GAQ T4_x = one - 0.5 * distance1x * (e2 ^ ei4);
	GAQ T5_x = one - 0.5 * distance1x * (e3 ^ ei5);
	GAQ T_x = T1_x * T2_x * T3_x * T4_x * T5_x; //Translator in x direction into origin

	GAQ T1_y = one - 0.5 * distance1y * (e2 ^ ei1);
	GAQ T2_y = one + 0.5 * distance1y * (e2 ^ ei2) - 0.25 * pow(distance1y, 2) * (ei1 ^ ei2);
	GAQ T3_y = one - 0.5 * distance1y * (e1 ^ ei4);
	GAQ T4_y = one - 0.5 * distance1y * (e3 ^ ei6);
	GAQ T_y = T1_y * T2_y * T3_y * T4_y; //Translator in y direction into origin

	GAQ r1 = e1 ^ e2;
	GAQ r2 = eo6 ^ ei5;
	GAQ r3 = ei6 ^ eo5;
	GAQ r4 = 2 * (eo4 ^ ei2);
	GAQ r5 = 2 * (ei4 ^ eo2);
	GAQ r6 = eo4 ^ ei3;
	GAQ r = r1 + r2 + r3 + r4 + r5 + r6;
	GAQ R = r.RotorExponential(40, phi); //Rotor in the xy-plane

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
	GAQ r1 = e1 ^ e2;
	GAQ r2 = eo6 ^ ei5;
	GAQ r3 = ei6 ^ eo5;
	GAQ r4 = 2 * (eo4 ^ ei2);
	GAQ r5 = 2 * (ei4 ^ eo2);
	GAQ r6 = eo4 ^ ei3;

	GAQ r = r1 + r2 + r3 + r4 + r5 + r6;

	double phi = std::numbers::pi / 4.0;
	GAQ R1 = cos(phi / 2) * one + sin(phi / 2) * r1;
	GAQ R2 = cos(phi / 2) * one + sin(phi / 2) * r2;
	GAQ R3 = cos(phi / 2) * one + sin(phi / 2) * r3;
	GAQ R4 = cos(phi) * one + sin(phi) * (0.5 * r4);
	GAQ R5 = cos(phi) * one + sin(phi) * (0.5 * r5);


	GAQ C = Up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1-1.5*ei2-4*ei3+2*ei4+3*ei5+6*ei6
	GAQ c_cga = eo1 + e1 + 2 * e2 + 3 * e3 + 7 * ei1;
	GAQ c_24 = -1.5 * ei2 + 2 * ei4;
	GAQ c_56 = 3 * ei5 + 6 * ei6;
	GAQ c_3 = -4 * ei3;

	double x = 1;
	double y = 2;
	double z = 3;
	double theta = atan(2) - phi;
	GAQ CC = Up(sqrt(5) * cos(theta), sqrt(5) * sin(theta), 3);
	GAQ _rotated = (R1 * c_cga * ~R1) + ((R2 ^ R3) * c_56 * (~R3 ^ ~R2)) + c_3 + ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + (-0.5 * sin(phi) * sin(phi) * (x * x - y * y) + sin(phi) * cos(phi) * x * y) * ei3;
	std::cout << "Rotated: " << _rotated << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (_rotated == CC) << std::endl;

	GAQ rotor = r.RotorExponential(20, phi);
	GAQ rotated = (rotor * C * ~rotor)[1];
	std::cout << "Rotated: " << rotated << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rotated == CC) << std::endl;
}

void RotorXZ()
{
	const GAQ r1 = -1 * e3 ^ e1; //rotating in opposite direction, thus -1*, viz corresponding matrices
	const GAQ r2 = -1 * ei4 ^ eo6;
	const GAQ r3 = -2 * (ei3 ^ eo5);
	const GAQ r4 = -1 * (ei2 ^ eo5);
	const GAQ r5 = -1 * eo4 ^ ei6;
	const GAQ r6 = -2 * (eo3 ^ ei5);

	const GAQ r = r1 + r2 + r3 + r4 + r5 + r6;
	const double phi = std::numbers::pi / 4.0;


	const GAQ C = Up(1, 3, 2); //eo1+e1+3*e2+2*e3+7*ei1-2*ei2-3*ei3+3*ei4+2*ei5+6*ei6
	const GAQ c_euc = eo1 + e1 + 3 * e2 + 2 * e3 + 7 * ei1;
	const GAQ c_35 = -3 * ei3 + 2 * ei5;
	const GAQ c_2 = -2 * ei2;
	const GAQ c_46 = 3 * ei4 + 6 * ei6;

	const double y = 3;
	const double theta = atan(2) - phi;
	const GAQ CC = Up(sqrt(5) * cos(theta), 3, sqrt(5) * sin(theta));
	const GAQ rotor = r.RotorExponential(20, phi);
	const GAQ rotated = (rotor * C * ~rotor)[1];
	std::cout << "Rotated: " << rotated << std::endl;
	std::cout << " Target: " << CC << std::endl;
	std::cout << "   Good: " << (rotated == CC) << std::endl;
}

void RotorYZ()
{
	const GAQ r1 = e2 ^ e3;
	const GAQ r2 = eo6 ^ ei3;
	const GAQ r3 = ei5 ^ eo4;
	const GAQ r4 = 2 * (ei6 ^ eo3);
	const GAQ r5 = 2 * (eo2 ^ ei6);
	const GAQ r6 = ei2 ^ eo6;
	const GAQ r7 = eo5 ^ ei4;

	const GAQ r = r1 + r2 + r3 + r4 + r5 + r6 + r7;
	const long double phi = std::numbers::pi / 4;

	const GAQ R1 = cos(phi / 2) * one + sin(phi / 2) * r1;
	const GAQ R2 = cos(phi / 2) * one + sin(phi / 2) * r2;
	const GAQ R3 = cos(phi / 2) * one + sin(phi / 2) * r3;
	const GAQ R4 = cos(phi / 2) * one + sin(phi / 2) * r4;
	const GAQ R5 = cos(phi / 2) * one + sin(phi / 2) * r5;
	const GAQ R6 = cos(phi / 2) * one + sin(phi / 2) * r6;
	const GAQ R7 = cos(phi / 2) * one + sin(phi / 2) * r7;

	const GAQ C = Up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1-1.5*ei2-4*ei3+2*ei4+3*ei5+6*ei6
	const GAQ c_cga = eo1 + e1 + 2 * e2 + 3 * e3 + 7 * ei1;
	const GAQ c_45 = 2 * ei4 + 3 * ei5;
	const GAQ c_236 = -1.5 * ei2 - 4 * ei3 + 6 * ei6;


	const GAQ target = Up(1, 0.5 * 5 * sqrt(2), 0.5 * sqrt(2));
	const GAQ rot_236 = (r2 + r4 + r5 + r6).RotorExponential(20, phi);

	const GAQ _rotated = ((R1 * c_cga * ~R1) + ((R3 ^ R7) * c_45 * (~R7 ^ ~R3)) + (rot_236 * c_236 * ~rot_236))[1];
	std::cout << "Rotated: " << _rotated << std::endl;
	std::cout << " Target: " << target << std::endl;
	std::cout << "   Good: " << (_rotated == target) << std::endl;


	const GAQ rotor = r.RotorExponential(20, phi);
	const GAQ rotated = (rotor * C * ~rotor)[1];
	std::cout << "Rotated: " << rotated << std::endl;
	std::cout << " Target: " << target << std::endl;
	std::cout << "   Good: " << (rotated == target) << std::endl;
}

void TranslatorX()
{
	const GAQ t1 = -1 * e1 ^ ei1;
	const GAQ t2 = -1 * e1 ^ ei2;
	const GAQ t3 = -1 * e1 ^ ei3;
	const GAQ t4 = -1 * e2 ^ ei4;
	const GAQ t5 = -1 * e3 ^ ei5;

	const GAQ t = t1 + t2 + t3 + t4 + t5;
	const double distance = -1 / double(3);

	const GAQ T1 = one - 0.5 * distance * (e1 ^ ei1);
	const GAQ T2 = one - 0.5 * distance * (e1 ^ ei2) + 0.25 * pow(distance, 2) * (ei1 ^ ei2);
	const GAQ T3 = one - 0.5 * distance * (e1 ^ ei3) + 0.25 * pow(distance, 2) * (ei1 ^ ei3) + 0.25 * pow(distance, 2) * (ei2 ^ ei3);
	const GAQ T4 = one - 0.5 * distance * (e2 ^ ei4);
	const GAQ T5 = one - 0.5 * distance * (e3 ^ ei5);

	const GAQ translatorX = t.TranslatorExponential(3, distance);
	const GAQ translatorX2 = T1 * T2 * T3 * T4 * T5;

	const GAQ C = Up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	const GAQ target = Up(1 + distance, 2, 3);

	const GAQ translated = translatorX * C * ~translatorX;
	const GAQ translated2 = translatorX2 * C * ~translatorX2;

	std::cout << "Translated: " << translated << std::endl;
	std::cout << "Translated: " << translated2 << std::endl;
	std::cout << "    Target: " << target << std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
	std::cout << "      Good: " << (translated2 == target) << std::endl;
}

void TranslatorY()
{
	const GAQ t1 = -1 * e2 ^ ei1;
	const GAQ t2 = 1 * e2 ^ ei2;
	const GAQ t3 = -1 * e1 ^ ei4;
	const GAQ t4 = -1 * e3 ^ ei6;

	const GAQ t = t1 + t2 + t3 + t4;
	const int distance = 7;

	const GAQ T1 = one - 0.5 * distance * (e2 ^ ei1);
	const GAQ T2 = one + 0.5 * distance * (e2 ^ ei2) - 0.25 * pow(distance, 2) * (ei1 ^ ei2);
	const GAQ T3 = one - 0.5 * distance * (e1 ^ ei4);
	const GAQ T4 = one - 0.5 * distance * (e3 ^ ei6);

	const GAQ translatorY = t.TranslatorExponential(2, distance);
	const GAQ translatorY2 = T1 * T2 * T3 * T4;

	const GAQ C = Up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	const GAQ target = Up(1, 2 + distance, 3);

	const GAQ translated = translatorY * C * ~translatorY;
	const GAQ translated2 = translatorY2 * C * ~translatorY2;

	std::cout << "Translated: " << translated << std::endl;
	std::cout << "Translated: " << translated2 << std::endl;
	std::cout << "    Target: " << target << std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
	std::cout << "      Good: " << (translated2 == target) << std::endl;
}

void TranslatorZ()
{
	const GAQ t1 = -1 * e3 ^ ei1;
	const GAQ t2 = 1 * e3 ^ ei3;
	const GAQ t3 = -1 * e1 ^ ei5;
	const GAQ t4 = -1 * e2 ^ ei6;

	const GAQ t = t1 + t2 + t3 + t4;
	const int distance = -4;

	const GAQ T1 = one - 0.5 * distance * (e3 ^ ei1);
	const GAQ T2 = one + 0.5 * distance * (e3 ^ ei3) - 0.25 * pow(distance, 2) * (ei1 ^ ei3);
	const GAQ T3 = one - 0.5 * distance * (e1 ^ ei5);
	const GAQ T4 = one - 0.5 * distance * (e2 ^ ei6);

	const GAQ translatorZ = t.TranslatorExponential(2, distance);
	const GAQ translatorZ2 = T1 * T2 * T3 * T4;

	const GAQ C = Up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
	const GAQ target = Up(1, 2, 3 + distance);

	const GAQ translated = translatorZ * C * ~translatorZ;
	const GAQ translated2 = translatorZ2 * C * ~translatorZ2;

	std::cout << "Translated: " << translated << std::endl;
	std::cout << "Translated: " << translated2 << std::endl;
	std::cout << "    Target: " << target << std::endl;
	std::cout << "      Good: " << (translated == target) << std::endl;
	std::cout << "      Good: " << (translated2 == target) << std::endl;
}

GAQ Com(const GAQ& a, const GAQ& b)
{
	return 0.5 * ((a * b) - (b * a));
}