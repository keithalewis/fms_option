// fms_option.t.cpp - test fms option
#include <cassert>
#include <functional>
#include <iostream>
#include <utility>
#include "fms_test.h"
#include "fms_option.h"

namespace fms::variate {

	// Hermite polynomials H_0(x) = 1, H_1(x) = x, H_{n+1}(x) = x H_n(x) - n H_{n-1}(x)
	template<class X>
	inline constexpr X Hermite(size_t n, X x)
	{
		if (n == 0) {
			return X(1);
		}
		if (n == 1) {
			return x;
		}

		return x * Hermite(n - 1, x) - X(n - 1) * Hermite(n - 2, x);
	}


	// Normal mean 0 variance 1
	template<class X = double, class S = X>
	class normal
	{
		static constexpr X M_SQRT2 = X(1.41421356237309504880);
		static constexpr X M_SQRT2PI = X(2.50662827463100050240);
	public:
		typedef X xtype;
		typedef S stype;

		static X cdf(X x, S s = 0, size_t n = 0)
		{
			X x_ = x - s;

			if (n == 0) {
				return (1 + erf(x_ / X(M_SQRT2))) / 2;
			}

			X phi = exp(-x_ * x_ / X(2)) / X(M_SQRT2PI);

			if (n == 1) {
				return phi;
			}

			// (d/dx)^n phi(x) = (-1)^n phi(x) H_n(x)
			return phi * Hermite(n - 1, x_) * ((n & 1) ? 1 : -1);
		}

		// (d/ds) cdf(x, s, 0)
		static X edf(X x, S s)
		{
			return -cdf(x, s, 1);
		}

		static S cumulant(S s, size_t n = 0)
		{
			if (n == 0) {
				return s * s / 2;
			}
			if (n == 1) {
				return s;
			}
			if (n == 2) {
				return 1;
			}

			return S(0);
		}
	};
}

using namespace fms;
using namespace fms::variate;

template<class X>
int test_option_normal_value()
{
	X eps = std::numeric_limits<X>::epsilon();
	X f = X(100);
	X s = X(0.1); // 3-month 20% volatility
	X k = X(100);
	X p = X(3.9877611676744920);

	using payoff::put;
	using payoff::call;

	X x;
	{
		variate::normal<X> n;
		option m(n);
		x = m.moneyness(f, s, k);
		x -= X(0.05);
		ensure (fabs(x) < eps);
		x = m.value(f, s, put(k));
		x -= p;
		ensure(fabs(x) <= 10 * eps);
		X cp = m.value(f, s, call(k)) - m.value(f, s, put(k));
		ensure(cp == f - k);
	}
	{
		variate::normal<X, X> n;
		option m(n);
		auto v = [s, k, &m](X x) { return m.value(x, s, call(k)); };
		auto vf = [s, k, &m](X x) { return m.delta(x, s, call(k)); };
		X dx = X(0.01);
		auto [lo, hi] = test_derivative(v, vf, dx, X(90), X(110), X(1));
		assert(fabs(lo) < std::max(eps, 10 * dx * dx));
		assert(fabs(hi) < std::max(eps, 10 * dx * dx));
	}
	{
		variate::normal<X, X> n;
		option m(n);
		auto vf = [s, k, &m](X x) { return m.delta(x, s, call(k)); };
		auto vff = [s, k, &m](X x) { return m.gamma(x, s, k); };
		X dx = X(0.01);
		auto [lo, hi] = test_derivative(vf, vff, dx, X(90), X(110), X(1));
		assert(fabs(lo) < std::max(eps, 10 * dx * dx));
		assert(fabs(hi) < std::max(eps, 10 * dx * dx));
	}
	{
		variate::normal<X, X> n;
		option m(n);
		auto v = [f, k, &m](X x) { return m.value(f, x, call(k)); };
		auto vs = [f, k, &m](X x) { return m.vega(f, x, k); };
		X dx = X(0.01);
		auto [lo, hi] = test_derivative(v, vs, dx, X(.1), X(.2), X(.01));
		assert(fabs(lo) < std::max(eps, 10 * dx * dx));
		assert(fabs(hi) < std::max(eps, 10 * dx * dx));
	}

	return 0;
}
int test_option_normal_value_f = test_option_normal_value<float>();
int test_option_normal_value_d = test_option_normal_value<double>();

template<class X>
int test_option_payoff()
{
	X eps = std::numeric_limits<X>::epsilon();
	X f = X(100);
	X s = X(0.1); // 3-month 20% volatility
	X k = X(100);
	//X p = X(3.9877611676744920);

	{
		payoff::call c(k);
		assert(c.strike == k);
	}
	{
		variate::normal<X, X> n;
		option m(n);

		X x = m.moneyness(f, s, k);
		ensure(fabs(x - X(0.05)) < eps);

		X c = m.value(f, s, payoff::call(k));
		X p = m.value(f, s, payoff::call(k));
		assert(c - p == f - k);
		/*
		x -= p;
		ensure(fabs(x) <= 10 * eps);
		X cp = m.call_value(f, s, k) - m.put_value(f, s, k);
		ensure(cp == f - k);
		*/
	}

	return 0;
}
int test_option_payoff_d = test_option_payoff<double>();

template<class X>
int test_implied()
{
	{
		X f = 100;
		X s = 0.2;
		X k = 100;
		X s_;

		variate::normal<X> n;
		option m(n);
		X v = m.value(f, s, k);
		s_ = m.implied(f, v, k);
		s_ -= s;
	}

	return 0;
}
int test_implied_d = test_implied<double>();


int main()
{
	return 0;
}