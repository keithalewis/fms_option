// fms_option.t.cpp - test fms option
#include <cassert>
#include <functional>
#include <iostream>
#include <utility>
#include "fms_option.h"
#include "fms_variate_normal.h"
#include "fms_variate.h"

using namespace fms;
using namespace fms::variate;

// (f(x + h) - f(x + h))/2h = f'(x) + (1/6)f'''(x) h^2 + O(h^4)
template<class X = double>
inline X derivative(const std::function<X(X)>& f, X x, X dx)
{
	return (f(x + dx) - f(x - dx)) / (2 * dx);
}

int test_derivative_square()
{
	std::function<double(double)> sq = [](double x) { return x * x; };
	
	double dx = 1e-4;
	for (double x = -1; x < 1; x += 0.1) {
		double df = derivative<>(sq, x, dx);
		double df_ = 2 * x;
		double err = df - df_;
		ensure(fabs(err) < dx * dx * dx);
	}

	return 0;
}
int test_derivative_square_ = test_derivative_square();

// min and max error over [a, b) in steps of h
template<class X>
inline auto test_derivative(const std::function<X(X)>& f, const std::function<X(X)>& df, X dx, X a, X b, X h)
{
	X lo = std::numeric_limits<X>::max();
	X hi = -lo;

	for (X x = a; x < b; x += h) {
		X d = derivative(f, x, dx);
		X d_ = df(x);
		X err = d - d_;
		if (err < lo) {
			lo = err;
		}
		else if (err > hi) {
			hi = err;
		}
	}

	return std::pair(lo, hi);
}

template<class X>
int test_normal()
{
	variate::normal<X, X> n;

	ensure(n.cumulant(0) == 0); // true for all cumulants
	ensure(n.cumulant(0, 1) == 0); // mean
	ensure(n.cumulant(0, 2) == 1); // variance
	ensure(n.cumulant(0, 3) == 0);

	X dx = X(0.01);
	X s = 0;
	const std::function<X(X)>& f = [s, &n](X x) { return n.cdf(x, s); };
	const std::function<X(X)>& df = [s, &n](X x) { return n.cdf(x, s, 1); };
	X d = derivative(f, X(0), dx);
	X d_ = df(0);
	d -= d_;

	auto [lo, hi] = test_derivative<X>(f, df, dx, X(-1), X(1), X(0.01));
	ensure(fabs(lo) < dx * dx);
	ensure(fabs(hi) < dx * dx);

	return 0;
}

template<class D, class X>
int test_option_value()
{
	X eps = std::numeric_limits<X>::epsilon();
	X f = X(100);
	X s = X(0.1); // 3-month 20% volatility
	X k = X(100);
	X p = X(3.9877611676744920);

	X x;
	{
		variate::normal<X,X> n;
		option m(n);
		x = m.moneyness(f, s, k);
		x -= X(0.05);
		ensure (fabs(x) < eps);
		x = m.put_value(f, s, k);
		x -= p;
		ensure(fabs(x) <= 10 * eps);
		X cp = m.call_value(f, s, k) - m.put_value(f, s, k);
		ensure(cp == f - k);
	}
	
	{
		variate_model n(variate::normal<X, X>{});
		option m(n);
		x = m.moneyness(f, s, k);
		x -= X(0.05);
		ensure(fabs(x) < eps);
		x = m.put_value(f, s, k);
		x -= p;
		ensure(fabs(x) <= 10 * eps);
	}
	

	x = x;

	return 0;
}

int main()
{
	try {
		test_normal<float>();
		test_normal<double>();
		test_option_value<normal<float>,float>();
		test_option_value<normal<double>,double>();
		// test_option_delta
		// test_option_gamma
		// test_option_vega
		// test_option_implied

	}
	catch (const std::exception& ex) {
		std::cerr << ex.what() << std::endl;
	}

	return 0;
}