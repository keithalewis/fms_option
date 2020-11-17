// fms_option.t.cpp - test fms option
#include <cassert>
#include <functional>
#include <iostream>
#include <utility>
#include "fms_test.h"
#include "fms_variate_normal.h"
#include "fms_variate.h"

using namespace fms;
using namespace fms::variate;

int test_derivative_square()
{
	std::function<double(double)> sq = [](double x) { return x * x; };

	double dx = 1e-4;
	for (double x = -1; x < 1; x += 0.1) {
		double df = derivative<>(sq, x, dx);
		double df_ = 2 * x;
		double err = df - df_;
		assert(fabs(err) < dx * dx * dx);
	}

	return 0;
}
int test_derivative_square_ = test_derivative_square();

template<class X>
int test_variate_normal()
{
	X eps = 2*sqrt(std::numeric_limits<X>::epsilon());
	variate::normal<X, X> n;

	assert(n.cumulant(0) == 0); // true for all cumulants
	assert(n.cumulant(0, 1) == 0); // mean
	assert(n.cumulant(0, 2) == 1); // variance
	assert(n.cumulant(0, 3) == 0);

	X s = 0;
	const std::function<X(X)>& f = [s, &n](X x) { return n.cdf(x, s); };
	const std::function<X(X)>& df = [s, &n](X x) { return n.cdf(x, s, 1); };

	for (X dx : {X(0.01), X(0.001), X(0.0001)}) {
		auto [lo, hi] = test_derivative<X>(f, df, dx, X(-1), X(1), X(0.01));
		assert(fabs(lo) < std::max(eps, dx * dx));
		assert(fabs(hi) < std::max(eps, dx * dx));
	}

	return 0;
}
int test_variate_normal_f = test_variate_normal<float>();
int test_variate_normal_d = test_variate_normal<double>();
