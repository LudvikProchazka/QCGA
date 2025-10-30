#include "CppUnitTest.h"
#include "Examples.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace GAQ_MSUT
{
	TEST_CLASS(GAQ_MSUT)
	{
	public:
		
		TEST_METHOD(Test_RotationXY)
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
			GAQ target = Up(sqrt(5) * cos(theta), sqrt(5) * sin(theta), 3);
			GAQ _rotated = (R1 * c_cga * ~R1) + ((R2 ^ R3) * c_56 * (~R3 ^ ~R2)) + c_3 + ((R4 ^ R5) * c_24 * (~R5 ^ ~R4)) + (-0.5 * sin(phi) * sin(phi) * (x * x - y * y) + sin(phi) * cos(phi) * x * y) * ei3;

			GAQ rotor = r.RotorExponential(20, phi);
			GAQ rotated = (rotor * C * ~rotor)[1];

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

			const GAQ translatorX = t.TranslatorExponential(20, distance);
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
			const int distance = 7;

			const GAQ T1 = one - 0.5 * distance * (e2 ^ ei1);
			const GAQ T2 = one + 0.5 * distance * (e2 ^ ei2) - 0.25 * pow(distance, 2) * (ei1 ^ ei2);
			const GAQ T3 = one - 0.5 * distance * (e1 ^ ei4);
			const GAQ T4 = one - 0.5 * distance * (e3 ^ ei6);

			const GAQ translatorY = t.TranslatorExponential(20, distance);
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
			const int distance = -4;

			const GAQ T1 = one - 0.5 * distance * (e3 ^ ei1);
			const GAQ T2 = one + 0.5 * distance * (e3 ^ ei3) - 0.25 * pow(distance, 2) * (ei1 ^ ei3);
			const GAQ T3 = one - 0.5 * distance * (e1 ^ ei5);
			const GAQ T4 = one - 0.5 * distance * (e2 ^ ei6);

			const GAQ translatorZ = t.TranslatorExponential(20, distance);
			const GAQ translatorZ2 = T1 * T2 * T3 * T4;

			const GAQ C = Up(1, 2, 3); //eo1+e1+2*e2+3*e3+7*ei1+3*ei2-2*ei3+2*ei4+3*ei5+6*ei6
			const GAQ target = Up(1, 2, 3 + distance);

			const GAQ translated = translatorZ * C * ~translatorZ;
			const GAQ translated2 = translatorZ2 * C * ~translatorZ2;

			Assert::IsTrue(translated == target);
			Assert::IsTrue(translated2 == target);
		}
	};
}
