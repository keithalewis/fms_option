// fms_test.h - test helper functions
#pragma once
#include <functional>
#include <utility>

namespace fms {

	// (f(x + h) - f(x + h))/2h = f'(x) + f'''(x) h^2/3! + O(h^4)
	template<class X = double>
	inline X derivative(const std::function<X(X)>& f, X x, X dx)
	{
		return (f(x + dx) - f(x - dx)) / (2 * dx);
	}

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

}