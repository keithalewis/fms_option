// fms_option.h - option valuation and greeks
// https://keithalewis.github.io/math/op.html
/// !!! add digital options
#pragma once
#include <algorithm>
#include <cmath>
#include <concepts>
#include <functional>
#include <limits>
#include "fms_ensure.h"
#include "fms_payoff.h"

namespace fms {

	// Strike K is always scalar floating point
	template<class M,
		class F = typename M::type, class S = typename M::ctype,
		class X = std::common_type_t<F, S>>
	class option {
		const M& m;
	public:
		option(const M& m)
			: m(m)
		{ }
		option(const option&) = default;
		option& operator=(const option&) = default;
		~option()
		{ }

		template<class K>
		X moneyness(F f, S s, K k) const
		{
			ensure(f > 0);
			ensure(s > 0);
			ensure(k > 0);

			return (::log(k / f) + m.cumulant(s)) / s;
		}

		// c = p + f - k so d^2c/df^2 = d^2p/df^2
		template<class K>
		X gamma(F f, S s, K k) const
		{
			if (f == 0 or k == 0) {
				return X(0);
			}
			if (s == 0) {
				return f == k ? std::numeric_limits<X>::infinity() : X(0);
			}

			X x = moneyness(f, s, k);

			return m.cdf(x, s, 1) / (f * s);
		}

		// f - k = c - p so dc/ds = dp/ds
		template<class K>
		X vega(F f, S s, K k) const
		{
			auto x = moneyness(f, s, k);

			return f * f * m.cdf(x, s, 1) / k; // /sigma //!!! only for normal model
		}

		//
		// put
		//

		template<class K>
		X value(F f, S s, const payoff::put<K>& p) const
		{
			K k = p.strike;

			if (f == 0) {
				return X(0);
			}
			if (s == 0) { // intrinsic
				return (std::max)(k - f, X(0));
			}
			if (k == 0) {
				return X(0);
			}

			X x = moneyness(f, s, k);

			return k * m.cdf(x) - f * m.cdf(x, s);
		}
		template<class K>
		X delta(F f, S s, const payoff::put<K>& p) const
		{
			K k = p.strike;

			if (f == 0) {
				return X(0);
			}
			if (s == 0) {
				return -X(f <= k);
			}
			if (k == 0) {
				return X(0);
			}

			X x = moneyness(f, s, k);

			return -m.cdf(x, s);
		}
		template<class K>
		X gamma(F f, S s, const payoff::put<K>& p) const
		{
			return gamma(f, s, p.strike);
		}
		template<class K>
		X vega(F f, S s, const payoff::put<K>& p) const
		{
			return vega(f, s, p.strike);
		}
		// Vol matching option value using Newton-Raphson.
		template<class K>
		K implied(F f, K v, payoff::put<K> p, K s0 = 0, size_t n = 0, K eps = 0) const
		{
			static K epsilon = std::numeric_limits<K>::epsilon();

			if (s0 == 0) {
				s0 = K(0.1);
			}
			if (n == 0) {
				n = 10;
			}
			if (eps == 0) {
				eps = sqrt(epsilon);
			}
			else if (eps <= epsilon) {
				eps = 10 * epsilon;
			}

			K s_ = s0 + 2 * eps; // loop at least once
			while (fabs(s_ - s0) > eps) {
				s_ = s0 - (value(f, s0, p) - v) / vega(f, s0, p);
				std::swap(s_, s0);
				if (--n == 0) {
					break;
				}
			}
			ensure(n != 0);

			return s_;
		}

		//
		// call
		//

		template<class K>
		X value(F f, S s, const payoff::call<K>& c) const
		{
			K k = c.strike;

			if (f == 0) {
				return X(0);
			}
			if (s == 0) { // intrinsic
				return (std::max)(f - k, X(0));
			}
			if (k == 0) {
				return f;
			}

			X x = moneyness(f, s, k);

			return f * (1 - m.cdf(x, s)) - k * (1 - m.cdf(x));
		}
		template<class K>
		X delta(F f, S s, const payoff::call<K>& c) const
		{
			K k = c.strike;

			if (f == 0) {
				return X(0);
			}
			if (s == 0) {
				return X(f > k);
			}
			if (k == 0) {
				return X(1);
			}

			X x = moneyness(f, s, k);

			return X(1) - m.cdf(x, s);
		}
		template<class K>
		X gamma(F f, S s, const payoff::call<K>& c) const
		{
			return gamma(f, s, c.strike);
		}
		template<class K>
		X vega(F f, S s, const payoff::call<K>& c) const
		{
			return vega(f, s, c.strike);
		}
		template<class K>
		X implied(F f, K v, payoff::call<K> c, K s0 = 0, size_t n = 0, K eps = 0) const
		{
			K k = c.strike;

			// c - p = f - k so p = c - f + k
			return implied(f, v - f + k, payoff::put<K>(k), s0, n, eps);
		}

		//!!! digital_call, digital_put
	};

}
