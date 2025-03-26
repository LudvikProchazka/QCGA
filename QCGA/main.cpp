#include "QCGA.h"
#include "Blade.h"
#include <iostream>
#include <numbers>
#include <iomanip>      
#include <chrono>



QCGA com(const QCGA& a, const QCGA& b)
{
	return 0.5*((a * b) - (b * a));
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
	QCGA T = T1*T2*T3*T4*T5; //Translator in x direction

	double phi = std::numbers::pi / 4.0;
	QCGA r1 = e1 ^ e2;
	QCGA r2 = eo6 ^ ei5;
	QCGA r3 = ei6 ^ eo5;
	QCGA r4 = 2 * (eo4 ^ ei2);
	QCGA r5 = 2 * (ei4 ^ eo2);
	QCGA r6 = eo4 ^ ei3;
	QCGA r = r1 + r2 + r3 + r4 + r5 + r6;
	QCGA R = r.rotorExponential(30, phi);

	Blade transformed = (R * (T * Q * ~T) * ~R)[1];

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
	const double distance = -1/double(3);

	const QCGA T1 = one - 0.5 * distance * (e1 ^ ei1);
	const QCGA T2 = one - 0.5 * distance * (e1 ^ ei2) + 0.25 * pow(distance, 2) * (ei1 ^ ei2);
	const QCGA T3 = one - 0.5 * distance * (e1 ^ ei3) + 0.25 * pow(distance, 2) * (ei1 ^ ei3) + 0.25 * pow(distance, 2) * (ei2 ^ ei3);
	const QCGA T4 = one - 0.5 * distance * (e2 ^ ei4);
	const QCGA T5 = one - 0.5 * distance * (e3 ^ ei5);

	const QCGA translatorX = t.translatorExponential(20, distance);
	const QCGA translatorX2 = T1 * T2	* T3 * T4 * T5;

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
	const QCGA translatorY2 = T1*T2*T3*T4;

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
	const QCGA translatorZ2 = T1*T2*T3*T4;

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

int main()
{
	auto start = std::chrono::high_resolution_clock::now();
	QCGA::generateGeneratingBlades(); //creater an array of basis vectors in R^{9,6}... one, e1,e2,...,e15

	//RotationExample();
	//OPNS_IPNS_Duality();
	Elipsoid();

	//RotorXY();
	//RotorXZ();
	//RotorYZ();
	//TranslatorX();
	//TranslatorY();
	//TranslatorZ();

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
}
