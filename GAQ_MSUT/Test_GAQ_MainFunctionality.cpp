#include "CppUnitTest.h"
#include "Examples.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

// #define SKIP_MAIN_TESTS

namespace GAQ_MSUT
{
	TEST_CLASS(GAQ_MSUT)
	{
	public:
#ifndef SKIP_MAIN_TESTS
		
		TEST_METHOD(Test_RotationXY)
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

			Assert::IsTrue(_rotated == target);
			Assert::IsTrue(rotated == target);
		}

		TEST_METHOD(Test_RotationXZ)
		{
			const GAQ r1 = -1 * e3 ^ e1;
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
			const GAQ target = Up(sqrt(5) * cos(theta), 3, sqrt(5) * sin(theta));
			const GAQ rotor = r.RotorExponential(20, phi);
			const GAQ rotated = (rotor * C * ~rotor)[1];

			Assert::IsTrue(rotated == target);
		}

		TEST_METHOD(Test_RotationYZ)
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

			const GAQ rotor = r.RotorExponential(20, phi);
			const GAQ rotated = (rotor * C * ~rotor)[1];

			Assert::IsTrue(_rotated == target);
			Assert::IsTrue(rotated == target);
		}

		TEST_METHOD(Test_TranslationX)
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

			Assert::IsTrue(translated == target);
			Assert::IsTrue(translated2 == target);
		}

		TEST_METHOD(Test_TranslationY)
		{
			const GAQ t1 = -1 * e2 ^ ei1;
			const GAQ t2 = 1 * e2 ^ ei2;
			const GAQ t3 = -1 * e1 ^ ei4;
			const GAQ t4 = -1 * e3 ^ ei6;

			const GAQ t = t1 + t2 + t3 + t4;
			const double distance = 7;

			const GAQ T1 = one - 0.5 * distance * (e2 ^ ei1);
			const GAQ T2 = one + 0.5 * distance * (e2 ^ ei2) - 0.25 * pow(distance, 2) * (ei1 ^ ei2);
			const GAQ T3 = one - 0.5 * distance * (e1 ^ ei4);
			const GAQ T4 = one - 0.5 * distance * (e3 ^ ei6);

			const GAQ translatorY = t.TranslatorExponential(3, distance);
			const GAQ translatorY2 = T1 * T2 * T3 * T4;

			const GAQ C = Up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
			const GAQ target = Up(1, 2 + distance, 3);

			const GAQ translated = translatorY * C * ~translatorY;
			const GAQ translated2 = translatorY2 * C * ~translatorY2;

			Assert::IsTrue(translated == target);
			Assert::IsTrue(translated2 == target);
		}

		TEST_METHOD(Test_TranslationZ)
		{
			const GAQ t1 = -1 * e3 ^ ei1;
			const GAQ t2 = 1 * e3 ^ ei3;
			const GAQ t3 = -1 * e1 ^ ei5;
			const GAQ t4 = -1 * e2 ^ ei6;

			const GAQ t = t1 + t2 + t3 + t4;
			const double distance = -4;

			const GAQ T1 = one - 0.5 * distance * (e3 ^ ei1);
			const GAQ T2 = one + 0.5 * distance * (e3 ^ ei3) - 0.25 * pow(distance, 2) * (ei1 ^ ei3);
			const GAQ T3 = one - 0.5 * distance * (e1 ^ ei5);
			const GAQ T4 = one - 0.5 * distance * (e2 ^ ei6);

			const GAQ translatorZ = t.TranslatorExponential(3, distance);
			const GAQ translatorZ2 = T1 * T2 * T3 * T4;

			const GAQ C = Up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
			const GAQ target = Up(1, 2, 3 + distance);

			const GAQ translated = translatorZ * C * ~translatorZ;
			const GAQ translated2 = translatorZ2 * C * ~translatorZ2;

			Assert::IsTrue(translated == target);
			Assert::IsTrue(translated2 == target);
		}

		TEST_METHOD(Test_Rotate)
		{
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

			Assert::IsTrue(rotatedXY == targetXY);
			Assert::IsTrue(rotatedXZ == targetXZ);
			Assert::IsTrue(rotatedYZ == targetYZ);
		}

		TEST_METHOD(Test_Translate)
		{
			const double distanceX = -1.0 / 3.0;
			const double distanceY = 7.0;
			const double distanceZ = -4.0;

			const GAQ point = Up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1-1.5*ei2-4*ei3+2*ei4+3*ei5+6*ei6
			const GAQ translatedXY = GAQ::Translate(point, x, distanceX);
			const GAQ translatedXZ = GAQ::Translate(point, y, distanceY);
			const GAQ translatedYZ = GAQ::Translate(point, z, distanceZ);
			const GAQ targetX = Up(1 + distanceX, 2, 3);
			const GAQ targetY = Up(1, 2 + distanceY, 3);
			const GAQ targetZ = Up(1, 2, 3 + distanceZ);

			Assert::IsTrue(translatedXY == targetX);
			Assert::IsTrue(translatedXZ == targetY);
			Assert::IsTrue(translatedYZ == targetZ);
		}
#endif
	};
}
