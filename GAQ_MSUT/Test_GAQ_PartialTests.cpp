#include "CppUnitTest.h"
#include "GAQ.h"
#include <numbers>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace GAQ_MSUT
{
	TEST_CLASS(GAQ_Partial)
	{
	public:
		TEST_METHOD(Test_GradeProjection)
		{
			GAQ::GenerateGeneratingBlades();

			GAQ zero1 = 1.41421356237 * zero_vector;
			GAQ zero2 = 1.41421356237 * one;
			GAQ _1 = e1;
			GAQ _2 = 2 * e1 ^ e2;
			GAQ _3 = 3 * e1 ^ e2 ^ e3;
			GAQ _4 = 4 * e1 ^ e2 ^ e3 ^ e4;
			GAQ _5 = 5 * e1 ^ e2 ^ e3 ^ e4 ^ e5;
			GAQ _6 = 6 * e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6;
			GAQ _7 = 7 * e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7;
			GAQ _8 = 8 * e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8;
			GAQ _9 = 9 * e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8 ^ e9;
			GAQ _10 = 10 * e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8 ^ e9 ^ e10;
			GAQ _11 = 11 * e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8 ^ e9 ^ e10 ^ e11;
			GAQ _12 = 12 * e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8 ^ e9 ^ e10 ^ e11 ^ e12;
			GAQ _13 = 13 * e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8 ^ e9 ^ e10 ^ e11 ^ e12 ^ e13;
			GAQ _14 = 14 * e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8 ^ e9 ^ e10 ^ e11 ^ e12 ^ e13 ^ e14;
			GAQ _15 = 15 * e1 ^ e2 ^ e3 ^ e4 ^ e5 ^ e6 ^ e7 ^ e8 ^ e9 ^ e10 ^ e11 ^ e12 ^ e13 ^ e14 ^ e15;
			GAQ allGrades = zero1 + zero2 + _1 + _2 + _3 + _4 + _5 + _6 + _7 + _8 + _9 + _10 + _11 + _12 + _13 + _14 + _15;

			Assert::IsTrue(allGrades[0] == 1.41421356237 * one);
			Assert::IsTrue(allGrades[1] == _1);
			Assert::IsTrue(allGrades[2] == _2);
			Assert::IsTrue(allGrades[3] == _3);
			Assert::IsTrue(allGrades[4] == _4);
			Assert::IsTrue(allGrades[5] == _5);
			Assert::IsTrue(allGrades[6] == _6);
			Assert::IsTrue(allGrades[7] == _7);
			Assert::IsTrue(allGrades[8] == _8);
			Assert::IsTrue(allGrades[9] == _9);
			Assert::IsTrue(allGrades[10] == _10);
			Assert::IsTrue(allGrades[11] == _11);
			Assert::IsTrue(allGrades[12] == _12);
			Assert::IsTrue(allGrades[13] == _13);
			Assert::IsTrue(allGrades[14] == _14);
			Assert::IsTrue(allGrades[15] == _15);
		}

		TEST_METHOD(Test_NotEqualOperator)
		{
			GAQ::GenerateGeneratingBlades();

			double distance = std::numbers::pi;
			GAQ precise = Tx(distance);
			GAQ closeEnough = precise * (1 + 1e-10);
			GAQ outOfPrecision = precise * (1 + 1e-8);

			Assert::IsFalse(precise != closeEnough);
			Assert::IsTrue(precise != outOfPrecision);
		}

		TEST_METHOD(Test_EqualEqualOperator)
		{
			GAQ::GenerateGeneratingBlades();

			double distance = std::numbers::pi;
			GAQ precise = Tx(distance);
			GAQ closeEnough = precise * (1 + 1e-10);
			GAQ outOfPrecision = precise * (1 + 1e-8);

			Assert::IsTrue(precise == closeEnough);
			Assert::IsFalse(precise == outOfPrecision);
		}

		TEST_METHOD(Test_ToNumeric)
		{
			GAQ::GenerateGeneratingBlades();

			GAQ a = 2 * e1 * e2 + 3 * e1 * e2 * e3;
			GAQ b = -7 * one + 2 * e1 * e2 + 3 * e1 * e2 * e3;
			GAQ c;

			Assert::IsTrue(a.ToNumeric() == 1.0);
			Assert::IsTrue(b.ToNumeric() == - 7.0);
			Assert::IsTrue(c.ToNumeric() == 0.0);
		}

		TEST_METHOD(Test_RotorExponential)
		{
			GAQ::GenerateGeneratingBlades();

			GAQ rotorXY = rxy.RotorExponential(20, std::numbers::pi / 4);

			GAQ expected =
				0.39429025373762072 * one
				+ 0.16332037060753438 * (e1 * e2)
				- 0.16332037060929533 * (e1 * e2 * e11 * e13)
				+ 0.067649512520283811 * (e1 * e2 * e11 * e13 * e14 * e15)
				- 0.040830092652323832 * (e1 * e2 * e12 * e13)
				+ 0.016912378130070953 * (e1 * e2 * e12 * e13 * e14 * e15)
				- 0.067649512518023203 * (e1 * e2 * e14 * e15)
				+ 0.040830092652889921 * (e1 * e2 * e5 * e6 * e7 * e13)
				- 0.016912378129631506 * (e1 * e2 * e5 * e6 * e7 * e13 * e14 * e15)
				+ 0.016912378129631506 * (e1 * e2 * e5 * e6 * e7 * e8 * e9 * e13)
				- 0.0070053363927480138 * (e1 * e2 * e5 * e6 * e7 * e8 * e9 * e13 * e14 * e15)
				+ 0.16332037060929533 * (e1 * e2 * e5 * e7)
				- 0.16332037061155971 * (e1 * e2 * e5 * e7 * e11 * e13)
				+ 0.067649512518526023 * (e1 * e2 * e5 * e7 * e11 * e13 * e14 * e15)
				- 0.040830092652889928 * (e1 * e2 * e5 * e7 * e12 * e13)
				+ 0.016912378129631506 * (e1 * e2 * e5 * e7 * e12 * e13 * e14 * e15)
				- 0.067649512520283811 * (e1 * e2 * e5 * e7 * e14 * e15)
				+ 0.067649512520283825 * (e1 * e2 * e5 * e7 * e8 * e9)
				- 0.067649512518526023 * (e1 * e2 * e5 * e7 * e8 * e9 * e11 * e13)
				+ 0.028021345570992055 * (e1 * e2 * e5 * e7 * e8 * e9 * e11 * e13 * e14 * e15)
				- 0.016912378129631506 * (e1 * e2 * e5 * e7 * e8 * e9 * e12 * e13)
				+ 0.0070053363927480138 * (e1 * e2 * e5 * e7 * e8 * e9 * e12 * e13 * e14 * e15)
				- 0.028021345573248927 * (e1 * e2 * e5 * e7 * e8 * e9 * e14 * e15)
				- 0.040830092652323832 * (e1 * e2 * e6 * e13)
				+ 0.016912378130070953 * (e1 * e2 * e6 * e13 * e14 * e15)
				+ 0.040830092652323832 * (e1 * e2 * e6 * e7)
				- 0.040830092652889928 * (e1 * e2 * e6 * e7 * e11 * e13)
				+ 0.016912378129631506 * (e1 * e2 * e6 * e7 * e11 * e13 * e14 * e15)
				- 0.016912378130070953 * (e1 * e2 * e6 * e7 * e14 * e15)
				+ 0.016912378130070956 * (e1 * e2 * e6 * e7 * e8 * e9)
				- 0.016912378129631506 * (e1 * e2 * e6 * e7 * e8 * e9 * e11 * e13)
				+ 0.0070053363927480138 * (e1 * e2 * e6 * e7 * e8 * e9 * e11 * e13 * e14 * e15)
				- 0.0070053363933122318 * (e1 * e2 * e6 * e7 * e8 * e9 * e14 * e15)
				- 0.016912378130070956 * (e1 * e2 * e6 * e8 * e9 * e13)
				+ 0.0070053363933122318 * (e1 * e2 * e6 * e8 * e9 * e13 * e14 * e15)
				- 0.040830092652889928 * (e1 * e2 * e7 * e11 * e12 * e13)
				+ 0.016912378129631506 * (e1 * e2 * e7 * e11 * e12 * e13 * e14 * e15)
				- 0.040830092652323832 * (e1 * e2 * e7 * e12)
				+ 0.016912378130070953 * (e1 * e2 * e7 * e12 * e14 * e15)
				- 0.016912378129631506 * (e1 * e2 * e7 * e8 * e9 * e11 * e12 * e13)
				+ 0.0070053363927480138 * (e1 * e2 * e7 * e8 * e9 * e11 * e12 * e13 * e14 * e15)
				- 0.016912378130070956 * (e1 * e2 * e7 * e8 * e9 * e12)
				+ 0.0070053363933122318 * (e1 * e2 * e7 * e8 * e9 * e12 * e14 * e15)
				+ 0.067649512518023203 * (e1 * e2 * e8 * e9)
				- 0.067649512520283825 * (e1 * e2 * e8 * e9 * e11 * e13)
				+ 0.028021345573248927 * (e1 * e2 * e8 * e9 * e11 * e13 * e14 * e15)
				- 0.016912378130070956 * (e1 * e2 * e8 * e9 * e12 * e13)
				+ 0.0070053363933122318 * (e1 * e2 * e8 * e9 * e12 * e13 * e14 * e15)
				- 0.028021345575003562 * (e1 * e2 * e8 * e9 * e14 * e15)
				- 0.39429025373535265 * (e11 * e13)
				+ 0.16332037060929530 * (e11 * e13 * e14 * e15)
				- 0.098572563433838162 * (e12 * e13)
				+ 0.040830092652323825 * (e12 * e13 * e14 * e15)
				- 0.16332037060753440 * (e14 * e15)
				+ 0.098572563434279198 * (e5 * e6 * e7 * e13)
				- 0.040830092652889921 * (e5 * e6 * e7 * e13 * e14 * e15)
				+ 0.040830092652889928 * (e5 * e6 * e7 * e8 * e9 * e13)
				- 0.016912378129631506 * (e5 * e6 * e7 * e8 * e9 * e13 * e14 * e15)
				+ 0.39429025373535265 * (e5 * e7)
				- 0.39429025373711679 * (e5 * e7 * e11 * e13)
				+ 0.16332037061155968 * (e5 * e7 * e11 * e13 * e14 * e15)
				- 0.098572563434279198 * (e5 * e7 * e12 * e13)
				+ 0.040830092652889921 * (e5 * e7 * e12 * e13 * e14 * e15)
				- 0.16332037060929530 * (e5 * e7 * e14 * e15)
				+ 0.16332037060929533 * (e5 * e7 * e8 * e9)
				- 0.16332037061155971 * (e5 * e7 * e8 * e9 * e11 * e13)
				+ 0.067649512518526023 * (e5 * e7 * e8 * e9 * e11 * e13 * e14 * e15)
				- 0.040830092652889928 * (e5 * e7 * e8 * e9 * e12 * e13)
				+ 0.016912378129631506 * (e5 * e7 * e8 * e9 * e12 * e13 * e14 * e15)
				- 0.067649512520283825 * (e5 * e7 * e8 * e9 * e14 * e15)
				- 0.098572563433838134 * (e6 * e13)
				+ 0.040830092652323825 * (e6 * e13 * e14 * e15)
				+ 0.098572563433838134 * (e6 * e7)
				- 0.098572563434279198 * (e6 * e7 * e11 * e13)
				+ 0.040830092652889921 * (e6 * e7 * e11 * e13 * e14 * e15)
				- 0.040830092652323825 * (e6 * e7 * e14 * e15)
				+ 0.040830092652323832 * (e6 * e7 * e8 * e9)
				- 0.040830092652889928 * (e6 * e7 * e8 * e9 * e11 * e13)
				+ 0.016912378129631506 * (e6 * e7 * e8 * e9 * e11 * e13 * e14 * e15)
				- 0.016912378130070956 * (e6 * e7 * e8 * e9 * e14 * e15)
				- 0.040830092652323832 * (e6 * e8 * e9 * e13)
				+ 0.016912378130070956 * (e6 * e8 * e9 * e13 * e14 * e15)
				- 0.098572563434279184 * (e7 * e11 * e12 * e13)
				+ 0.040830092652889921 * (e7 * e11 * e12 * e13 * e14 * e15)
				- 0.098572563433838162 * (e7 * e12)
				+ 0.040830092652323825 * (e7 * e12 * e14 * e15)
				- 0.040830092652889921 * (e7 * e8 * e9 * e11 * e12 * e13)
				+ 0.016912378129631506 * (e7 * e8 * e9 * e11 * e12 * e13 * e14 * e15)
				- 0.040830092652323832 * (e7 * e8 * e9 * e12)
				+ 0.016912378130070956 * (e7 * e8 * e9 * e12 * e14 * e15)
				+ 0.16332037060753438 * (e8 * e9)
				- 0.16332037060929533 * (e8 * e9 * e11 * e13)
				+ 0.067649512520283825 * (e8 * e9 * e11 * e13 * e14 * e15)
				- 0.040830092652323832 * (e8 * e9 * e12 * e13)
				+ 0.016912378130070956 * (e8 * e9 * e12 * e13 * e14 * e15)
				- 0.067649512518023203 * (e8 * e9 * e14 * e15);

			Assert::IsTrue(rotorXY == expected);
		}

		TEST_METHOD(Test_TranslatorExponential)
		{
			GAQ::GenerateGeneratingBlades();

			GAQ t1 = -1 * e1 ^ ei1;
			GAQ t2 = -1 * e1 ^ ei2;
			GAQ t3 = -1 * e1 ^ ei3;
			GAQ t4 = -1 * e2 ^ ei4;
			GAQ t5 = -1 * e3 ^ ei5;

			GAQ t = t1 + t2 + t3 + t4 + t5;
			double distance = std::numbers::pi;

			GAQ translatorX = t.TranslatorExponential(3, distance); // indeed only 3 because of commutation relations - zeros for higher orders
			GAQ expected = Tx(distance);

			Assert::IsTrue(translatorX == expected);
		}
	};
}
