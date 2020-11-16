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

template<class X = double>
inline X derivative(const std::function<X(X)>& f, X x, X dx)
{
	return (f(x + dx) - f(x - dx)) / (2 * dx);
}

// min and max error over [a, b) in steps of h
template<class X>
inline auto test_derivative(const std::function<X(X)>& f, const std::function<X(X)>& df, X dx, X a, X b, X h)
{
	X A = std::numeric_limits<X>::max();
	X B = -A;

	for (X x = a; x < b; x += h) {
		X err = derivative(f, x, dx) - df(x);
		if (err < A) {
			A = err;
		}
		else if (err > B) {
			B = err;
		}
	}

	return std::pair(A, B);
}

template<class X>
int test_normal()
{
	variate::normal<X, X> n;
	X s = 0;
	const auto& f = [s, &n](X x) { return n.cdf(x, s); };
	const auto& df = [s, &n](X x) { return n.cdf(x, s, 1); };
	auto [lo, hi] = test_derivative<X>(f, df, X(.0001), X(-1), X(1), X(0.01));
	//!!! failing on derivative

	return 0;
}

template<class D, class X>
int test_option()
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
		test_option<normal<float>,float>();
		test_option<normal<double>,double>();
	}
	catch (const std::exception& ex) {
		std::cerr << ex.what() << std::endl;
	}

	return 0;
}