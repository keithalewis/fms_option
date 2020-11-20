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

	// Strike K is always scalar
	template<class M,
		class F = typename M::type, class S = typename M::ctype,
		class X = std::common_type_t<F, S>>
	class opt {
		const M& m;
	public:
		opt(const M& m)
			: m(m)
		{ }
		opt(const opt&) = default;
		opt& operator=(const opt&) = default;
		~opt()
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
			K k = p.strike();

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
			K k = p.strike();

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
			return gamma(f, s, p.strike());
		}
		template<class K>
		X vega(F f, S s, const payoff::put<K>& p) const
		{
			return vega(f, s, p.strike());
		}
		// Vol matching option value using Newton-Raphson.
		template<class K>
		X implied(F f, K vp, payoff::put<K> p, S s0 = 0, size_t n = 0, S eps = 0) const
		{
			static S epsilon = std::numeric_limits<S>::epsilon();

			if (s0 == 0) {
				s0 = S(0.1);
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

			S s_ = s0 + 2 * eps; // loop at least once
			while (fabs(s_ - s0) > eps) {
				s_ = s0 - (value(f, s0, p) - vp) / vega(f, s0, p);
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
			K k = c.strike();

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
			K k = c.strike();

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
			return gamma(f, s, c.strike());
		}
		template<class K>
		X vega(F f, S s, const payoff::call<K>& c) const
		{
			return vega(f, s, c.strike());
		}
		template<class K>
		X implied(F f, K vc, payoff::call<K> c, S s0 = 0, size_t n = 0, S eps = 0) const
		{
			K k = c.strike();

			// c - p = f - k so p = c - f + k
			return implied(f, vc - f + k, payoff::put<K>(k), s0, n, eps);
		}

		//!!! digital_call, digital_put
	};

	/// <summary>
	/// Option value and greeks
	/// </summary>
	/// <typeparam name="M">Model of underlying</typeparam>
	/// <typeparam name="F">Forward</typeparam>
	/// <typeparam name="S">Vol</typeparam>
	/// <typeparam name="K">Strike</typeparam>
	template<class M,
		class F = typename M::type, class S = typename M::ctype, class K = typename M::type,
		class X = std::common_type_t<F, S, K>>
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

		// F = f exp(s X - kappa(s)) so
		// E[F] = f, and Var[log F] = s^2.
		// F <= k iff X <= (log(k/f) + cumulant(s))/s
		X moneyness(F f, S s, K k) const
		{
			ensure (f > 0);
			ensure (s > 0);
			ensure (k > 0);

			return (::log(k / f) + m.cumulant(s)) / s;
		}

		// p = k P(F <= k) - f P^s(F <= k)
		// call (k > 0) or put (k < 0) value
		// c = p + f - k
		X value(F f, S s, K k) const
		{
			X f_k;

			if (k < 0) { // put
				k = -k; 
				f_k = 0;
			}
			else { // call
				f_k = f - k;
			}

			if (f == 0) {
				return X(0);
			}
			if (s == 0) { // intrinsic
				return (std::max)(f_k == 0 ? k - f : f - k, X(0));
			}
			if (k == 0) {
				return f_k == 0 ? 0 : f;
			}

			X x = moneyness(f, s, k);

			return k * m.cdf(x) - f * m.cdf(x, s) + f_k;
		}
		X put_value(F f, S s, K k) const
		{
			return value(f, s, -k);
		}
		X call_value(F f, S s, K k) const
		{
			return value(f, s, k);
		}

		// dp/df = -P_s(F <= k)
		// call (k > 0) or put (k < 0) value
		// dc = dp + 1
		X delta(F f, S s, K k) const
		{
			X one = 1;

			if (k < 0) { // put
				k = -k;
				one = 0;
			}

			if (f == 0) {
				return X(0);
			}
			if (s == 0) {
				return one == 0 ? X(f <= k) : X(f > k);
			}
			if (k == 0) {
				return one;
			}

			X x = moneyness(f, s, k);

			return -m.cdf(x, s) + one;
		}
		X put_delta(F f, S s, K k) const
		{
			return delta(f, s, -k);
		}
		X call_delta(F f, S s, K k)
		{
			return delta(f, s, k);
		}

		// c = p + f - k so d^2c/df^2 = d^2p/df^2
		X gamma(F f, S s, K k) const
		{
			if (k < 0) {
				k = -k;
			}

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
		X vega(F f, S s, K k) const
		{
			if (k < 0) {
				k = -k;
			}

			auto x = moneyness(f, s, k);

			return f * f * m.cdf(x, s, 1)/k; // /sigma //!!! only for normal model
		}

		// Vol matching option price using Newton-Raphson.
		X implied(F f, S v, K k, S s0 = 0, size_t n = 0, S eps = 0) const
		{
			static S epsilon = std::numeric_limits<S>::epsilon();

			if (k < 0) { // put
				k = -k;
				v = v - f + k;
			}
			if (s0 == 0) {
				s0 = S(0.1);
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

			S s_ = s0 + 2*eps; // loop at least once
			while (fabs(s_ - s0) > eps) {
				s_ = s0 - (value(f, s0, k) - v) / vega(f, s0, k);
				std::swap(s_, s0);
				if (--n == 0) {
					break;
				}
			}
			ensure(n != 0);

			return s_;
		}
		X put_implied(F f, S p, K k, S s0 = 0, size_t n = 0, S eps = 0) const
		{
			return implied(f, p, -k, s0, n, eps);
		}
		X call_implied(F f, S c, K k, S s0 = 0, size_t n = 0, S eps = 0) const
		{
			return implied(f, c, k, s0, n, eps);
		}

	};

}
