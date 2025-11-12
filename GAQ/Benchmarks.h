#pragma once

#include "Blade.h"
#include "Instrumentor.h"
#include <numbers>

#define PROFILING 1
#if PROFILING
#define PROFILE_SCOPE(name) InstrumentationTimer timer##__LINE__(name)
#define PROFILE_FUNCTION() PROFILE_SCOPE(__FUNCSIG__)
#else
#define PROFILE_SCOPE(name)
#endif

static size_t gs_allocations{0};

void* operator new(size_t bytes)
{
	++gs_allocations;
	return malloc(bytes);
}

class Benchmarks : public GAQ
{
public:
	Benchmarks(const GAQ& other)
	{
		m_mapLabelToCoefficient = other.m_mapLabelToCoefficient;
	}

	static Benchmarks MakeQuadric(long double vo6, long double vo5, long double vo4, long double vo3, long double vo2, long double vo1, long double ve1, long double ve2, long double ve3, long double vi1)
	{
		PROFILE_FUNCTION();
		return gaq::MakeQuadric(vo6, vo5, vo4, vo3, vo2, vo1, ve1, ve2, ve3, vi1);
	}

	Benchmarks RotorExponential(unsigned int degree, long double phi) const
	{
		PROFILE_FUNCTION();
		return GAQ::RotorExponential(degree, phi);
	}

	Benchmarks operator*(const GAQ& other) const
	{
		PROFILE_SCOPE("GeometricProduct");
		return GAQ::operator*(other);
	}

	static Benchmarks TranslateRotateToOrigin(const Benchmarks& R, const Benchmarks& T_x, const Benchmarks& T_y, const Benchmarks& Q)
	{
		PROFILE_SCOPE("TranslateToOriginAndRotate");
		return (R * (T_y * (T_x * Q * ~T_x)[1] * ~T_y)[1] * ~R)[1];
	}

	static Benchmarks Transform(const Benchmarks& T_x, const Benchmarks& T_y, const Benchmarks& rotatedInOrigin)
	{
		PROFILE_SCOPE("TranslateToFinalPosition");
		return (T_y * (T_x * rotatedInOrigin * ~T_x)[1] * ~T_y)[1];
	}

	Benchmarks operator+(const GAQ& other) const
	{
		PROFILE_SCOPE("+");
		return GAQ::operator+(other);
	}

	static void Elipsoid()
	{
		PROFILE_FUNCTION();

		const Benchmarks Q = Benchmarks::MakeQuadric(0.0, 0.0, 0.0, 5.0 / 3.0, -1.0 / 3.0, -7.0 / 3.0, 5.0, 3.0, -5.0, 5.0);

		const double distance1x = 5.0;
		const double distance1y = 3.0 / 2;
		const double distance2x = 10.0;
		const double distance2y = 6.0;
		const double phi = 3.0 * std::numbers::pi / 4.0;

		Benchmarks T1_x = one - 0.5 * distance1x * (e1 ^ ei1);
		Benchmarks T2_x = one - 0.5 * distance1x * (e1 ^ ei2) + 0.25 * pow(distance1x, 2) * (ei1 ^ ei2);
		Benchmarks T3_x = one - 0.5 * distance1x * (e1 ^ ei3) + 0.25 * pow(distance1x, 2) * (ei1 ^ ei3) + 0.25 * pow(distance1x, 2) * (ei2 ^ ei3);
		Benchmarks T4_x = one - 0.5 * distance1x * (e2 ^ ei4);
		Benchmarks T5_x = one - 0.5 * distance1x * (e3 ^ ei5);
		Benchmarks T_x = T1_x * T2_x * T3_x * T4_x * T5_x;

		Benchmarks T1_y = one - 0.5 * distance1y * (e2 ^ ei1);
		Benchmarks T2_y = one + 0.5 * distance1y * (e2 ^ ei2) - 0.25 * pow(distance1y, 2) * (ei1 ^ ei2);
		Benchmarks T3_y = one - 0.5 * distance1y * (e1 ^ ei4);
		Benchmarks T4_y = one - 0.5 * distance1y * (e3 ^ ei6);
		Benchmarks T_y = T1_y * T2_y * T3_y * T4_y;

		const Benchmarks r1 = e1 ^ e2;
		const Benchmarks r2 = eo6 ^ ei5;
		const Benchmarks r3 = ei6 ^ eo5;
		const Benchmarks r4 = 2 * (eo4 ^ ei2);
		const Benchmarks r5 = 2 * (ei4 ^ eo2);
		const Benchmarks r6 = eo4 ^ ei3;
		const Benchmarks r = r1 + r2 + r3 + r4 + r5 + r6;
		const Benchmarks R = r.RotorExponential(40, phi);

		const Benchmarks rotatedInOrigin = TranslateRotateToOrigin(R, T_x, T_y, Q);

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

		Benchmarks transformed = Transform(T_y, T_y, rotatedInOrigin);
	}

	static void RunBenchmarks()
	{
		PROFILE_FUNCTION();

		Elipsoid();
	}
};