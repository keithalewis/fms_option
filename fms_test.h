// fms_test.h - test helper functions
#pragma once
#include <cassert>
#include <functional>
#include <utility>
#include "fms_variate.h"

// (f(x + h) - f(x - h))/2h = f'(x) + f'''(x) h^2/3! + O(h^4)
template<class X = double>
inline X derivative(auto& f, X x, X dx)
{
	return (f(x + dx) - f(x - dx)) / (2 * dx);
}

// min and max error over [a, b) in steps of h
template<class X>
inline auto test_derivative(auto& f, auto& df, X dx, X a, X b, X h)
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

template<class M, class X>
inline auto test_variate_derivative(const M& m, X dx, X s, X a, X b, X h, size_t n)
{
	auto f = [s, n, &m](X x) { return m.cdf(x, s, n); };
	auto df = [s, n, &m](X x) { return m.cdf(x, s, n + 1); };

	return test_derivative(f, df, dx, a, b, h);
}

template<class M, class X = M::xtype, class S = M::stype>
inline int test_variate(const M& m, X dx)
{
	X eps = 2 * sqrt(std::numeric_limits<X>::epsilon());
	S mu = m.cumulant(0, 1);
	S sigma = ::sqrt(m.cumulant(0, 2));
	for (size_t n : {0,1,2}) {
		for (S s : {S(-.1), S(0), S(1)}) {
			auto [lo, hi] = test_variate_derivative(m, dx, s, mu - 2 * sigma, mu + 2 * sigma, sigma / 10, n);
			assert(fabs(lo) < std::max(eps, dx * dx));
			assert(fabs(hi) < std::max(eps, dx * dx));
		}
	}

	return 0;
}
