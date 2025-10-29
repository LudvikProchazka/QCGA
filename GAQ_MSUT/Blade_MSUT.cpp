#include "CppUnitTest.h"
#include "Blade.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace GAQ_MSUT
{
	TEST_CLASS(Blade_MSUT)
	{
	public:

		TEST_METHOD(Test_Up)
		{
			GAQ::GenerateGeneratingBlades();

			double x{1.0};
			double y{2.0};
			double z{-3.0};
			Blade b = Up(x, y, z);
			Blade expected = eo1 + x * e1 + y * e2 + z * e3
				+ 1.0 / 2.0 * (x * x + y * y + z * z) * ei1
				+ 1.0 / 2.0 * (x * x - y * y) * ei2
				+ 1.0 / 2.0 * (x * x - z * z) * ei3
				+ x * y * ei4
				+ x * z * ei5
				+ y * z * ei6;

			Assert::IsTrue(b == expected);
		}

		TEST_METHOD(Test_MakeQuadric)
		{
			GAQ::GenerateGeneratingBlades();

			double v06{3.0};
			double v05{3.14};
			double v04{-11.9};
			double v03{7.43};
			double v02{0.0128};
			double v01{1.035};
			double ve1{1.0};
			double ve2{2.0};
			double ve3{-3.0};
			double vi1{-4.69};

			Blade Q = MakeQuadric(v06, v05, v04, v03, v02, v01, ve1, ve2, ve3, vi1);
			Blade expected = v06 * eo6 + v05 * eo5 + v04 * eo4 + v03 * eo3 + v02 * eo2 + v01 * eo1 + ve1 * e1 + ve2 * e2 + ve3 * e3 + vi1 * ei1;

			Assert::IsTrue(Q == expected);
		}

		TEST_METHOD(Test_IsBlade)
		{
			GAQ::GenerateGeneratingBlades();

			Blade a = e1*e2 + e1*e3;
			Blade b = e1*e2 + e3*e4;

			Assert::IsTrue(a.IsBlade());
			Assert::IsFalse(b.IsBlade());
		}

		TEST_METHOD(Test_Down)
		{
			GAQ::GenerateGeneratingBlades();

			Blade euclidean = 1.41421356237 * e1 + 2.71828 * e2 -3.14 * e3;
			Blade b = Up(1.41421356237, 2.71828, -3.14);
			Blade res = b.Down();

			Assert::IsTrue(res == euclidean);
			Assert::IsTrue(b.IsBlade());
		}

		TEST_METHOD(Test_Normalize)
		{
			GAQ::GenerateGeneratingBlades();

			Blade euclidean = 1.41421356237 * e1 + 2.71828 * e2 - 3.14 * e3;
			Blade b = log(13) * Up(1.41421356237, 2.71828, -3.14);
			Blade res = b.Normalize();

			Assert::IsTrue(res == Up(1.41421356237, 2.71828, -3.14));
			Assert::IsTrue(b.IsBlade());
		}

		TEST_METHOD(Test_Dual)
		{
			GAQ::GenerateGeneratingBlades();

			Blade a = e1 * e2 + e1 * e3 + e1 * e8;
			Blade a_dual = a * (I ^ (-1));
			Blade res = a.Dual();

			Assert::IsTrue(res == a_dual);
			Assert::IsTrue(a_dual.IsBlade());
		}

		TEST_METHOD(Test_ExponentPositive)
		{
			GAQ::GenerateGeneratingBlades();

			Blade a = e1 * e2 + e1 * e3 + e1 * e8;
			Blade res1 = a * a * a;
			Blade res2 = a^3;

			Assert::IsTrue(res1 == res2);
		}

		TEST_METHOD(Test_ExponentNegative)
		{
			GAQ::GenerateGeneratingBlades();

			Blade a = e1 * e2 + e1 * e3 + e1 * e8;
			Blade inversion = a ^ -1;
			Blade res = a ^ -3;

			Assert::IsTrue(inversion * a == one);
			Assert::IsTrue((res * inversion ^ 3) == inversion ^ 3);
		}

		TEST_METHOD(Test_ExponentInversion) // Propper usage
		{
			GAQ::GenerateGeneratingBlades();

			Blade a = e1 * e2 + e1 * e3 + e1 * e8;
			Blade res = a ^ -1;

			Assert::IsTrue(a * res == one);
		}

		TEST_METHOD(Test_NullBlade)
		{
			GAQ::GenerateGeneratingBlades();

			Blade euclidean = 1.41421356237 * e1 + 2.71828 * e2 - 3.14 * e3;
			Blade b = Up(1.41421356237, 2.71828, -3.14);

			Assert::IsTrue(b.IsBlade());
			Assert::IsTrue(b.IsNullBlade());
			Assert::IsFalse(euclidean.IsNullBlade());
		}

		TEST_METHOD(Test_Grade)
		{
			GAQ::GenerateGeneratingBlades();

			Blade zero1 = 1.41421356237 * zero_vector;
			Blade zero2 = 1.41421356237 * one;
			Blade _1 = e1;
			Blade _2 = e1 ^ e2;
			Blade _3 = e1 ^ e2 ^ e3;
			Blade _4 = e1 ^ e2 ^ e3 ^ e4;
			Blade _5 = e1 ^ e2 ^ e3 ^ e4 ^ e5;
			Blade _6 = e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6;
			Blade _7 = e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7;
			Blade _8 = e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8;
			Blade _9 = e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8 ^ e9;
			Blade _10 = e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8 ^ e9 ^ e10;
			Blade _11 = e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8 ^ e9 ^ e10 ^ e11;
			Blade _12 = e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8 ^ e9 ^ e10 ^ e11 ^ e12;
			Blade _13 = e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8 ^ e9 ^ e10 ^ e11 ^ e12 ^ e13;
			Blade _14 = e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8 ^ e9 ^ e10 ^ e11 ^ e12 ^ e13 ^ e14;
			Blade _15 = e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8 ^ e9 ^ e10 ^ e11 ^ e12 ^ e13 ^ e14 ^ e15;

			Assert::IsTrue(zero1.GetGrade() == 0);
			Assert::IsTrue(zero2.GetGrade() == 0);
			Assert::IsTrue(_1.GetGrade() == 1);
			Assert::IsTrue(_2.GetGrade() == 2);
			Assert::IsTrue(_3.GetGrade() == 3);
			Assert::IsTrue(_4.GetGrade() == 4);
			Assert::IsTrue(_5.GetGrade() == 5);
			Assert::IsTrue(_6.GetGrade() == 6);
			Assert::IsTrue(_7.GetGrade() == 7);
			Assert::IsTrue(_8.GetGrade() == 8);
			Assert::IsTrue(_9.GetGrade() == 9);
			Assert::IsTrue(_10.GetGrade() == 10);
			Assert::IsTrue(_11.GetGrade() == 11);
			Assert::IsTrue(_12.GetGrade() == 12);
			Assert::IsTrue(_13.GetGrade() == 13);
			Assert::IsTrue(_14.GetGrade() == 14);
			Assert::IsTrue(_15.GetGrade() == 15);
		}

		TEST_METHOD(Test_GradeCalculated)
		{
			GAQ::GenerateGeneratingBlades();

			Blade _11 = e1 * e2 * e3 * e4 * e5 * e6 * e7 * e8 * e9 * e10 * e11 * e12 * e13 * e14 * e15 * e1 * e3 * e4 * e5 * e11;

			Assert::IsTrue(_11.GetGrade() == 10);
		}
	};
}
