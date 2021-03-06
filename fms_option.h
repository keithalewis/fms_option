﻿// fms_option.h - option valuation and greeks
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
		class F = typename M::xtype, class S = typename M::stype,
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

#pragma region value

		template<class K>
		X value(F f, S s, const payoff::call<K>& c) const
		{
			K k = ::abs(c.strike);

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
		X value(F f, S s, const payoff::put<K>& p) const
		{
			auto k = abs(p.strike);

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
		// use negative strike for put
		template<class K>
		X value(F f, S s, K k) const
		{
			return k > 0 ? value(f, s, payoff::call(k)) : value(f, s, payoff::put(-k));
		}

		template<class K>
		X value(F f, S s, const payoff::digital_call<K>& c) const
		{
			K k = c.strike;

			if (f == 0) {
				return X(0);
			}
			if (s == 0) { // intrinsic
				return X(1)*(f > k);
			}
			if (k == 0) {
				return 1;
			}

			X x = moneyness(f, s, k);

			return 1 - m.cdf(x);
		}

		template<class K>
		X value(F f, S s, const payoff::digital_put<K>& p) const
		{
			auto k = p.strike;

			if (f == 0) {
				return X(0);
			}
			if (s == 0) { // intrinsic
				return X(1)*(f <= k);
			}
			if (k == 0) {
				return X(0);
			}

			X x = moneyness(f, s, k);

			return m.cdf(x);
		}

#pragma endregion // value

#pragma region delta

		template<class K>
		X delta(F f, S s, const payoff::call<K>& c) const
		{
			K k = ::abs(c.strike);

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
		X delta(F f, S s, const payoff::put<K>& p) const
		{
			K k = ::abs(p.strike);

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
		// negative strike indicates put
		template<class K>
		X delta(F f, S s, K k) const
		{
			return k > 0 ? delta(f, s, payoff::call(k)) : delta(f, s, payoff::put(-k));
		}

		template<class K>
		X delta(F f, S s, const payoff::digital_call<K>& c) const
		{
			return -delta(f, s, payoff::digital_put(c.strike));
		}
		template<class K>
		X delta(F f, S s, const payoff::digital_put<K>& p) const
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

			return -m.cdf(x, 0, 1)/(f*s);
		}

#pragma endregion // delta

#pragma region gamma

		// c = p + f - k so d^2c/df^2 = d^2p/df^2
		template<class K>
		X gamma(F f, S s, K k) const
		{
			k = ::abs(k);

			if (f == 0 or k == 0) {
				return X(0);
			}
			if (s == 0) {
				return f == k ? std::numeric_limits<X>::infinity() : X(0);
			}

			X x = moneyness(f, s, k);

			return m.cdf(x, s, 1) / (f * s);
		}
		template<class K>
		X gamma(F f, S s, const payoff::call<K>& c) const
		{
			return gamma(f, s, c.strike);
		}
		template<class K>
		X gamma(F f, S s, const payoff::put<K>& p) const
		{
			return gamma(f, s, p.strike);
		}

		template<class K>
		X gamma(F f, S s, const payoff::digital_call<K>& c) const
		{
			return -gamma(f, s, payoff::digital_put(c.strike));
		}
		template<class K>
		X gamma(F f, S s, const payoff::digital_put<K>& p) const
		{
			K k = p.strike;

			if (f == 0) {
				return X(0);
			}
			if (s == 0) {
				return f == k ? std::numeric_limits<X>::infinity() : 0;
			}
			if (k == 0) {
				return X(0);
			}

			X x = moneyness(f, s, k);

			return (m.cdf(x, 0, 1) * s - m.cdf(x, 0, 2)) / (f * f * s * s);
		}

#pragma endregion // gamma

#pragma region vega

		// f - k = c - p so dc/ds = dp/ds
		template<class K>
		X vega(F f, S s, K k) const
		{
			k = ::abs(k);

			auto x = moneyness(f, s, k);

			return -f * m.edf(s, x); 
		}
		template<class K>
		X vega(F f, S s, const payoff::call<K>& c) const
		{
			return vega(f, s, c.strike);
		}
		template<class K>
		X vega(F f, S s, const payoff::put<K>& p) const
		{
			return vega(f, s, p.strike);
		}

		template<class K>
		X vega(F f, S s, const payoff::digital_call<K>& c) const
		{
			return -vega(f, s, payoff::digital_put(c.strike));
		}
		template<class K>
		X vega(F f, S s, const payoff::digital_put<K>& p) const
		{
			K k = p.strike;

			if (f == 0) {
				return X(0);
			}
			if (k == 0) {
				return X(0);
			}

			X x = moneyness(f, s, k);

			return -m.cdf(x, 0, 1)*(x - s)/s;
		}


#pragma endregion // vega

		/*
		// If we know the implied vol is s then if v > v0 where v0 is
		// the at-the-money value it must be a call if f > k and a put
		// if f < k. If v < v0 then it must be a put if f > k and a call
		// if f < k. The exclusive or of v > v0 and f > k is false for
		// a call and true for a put.
		// We don't know s so we use the current best guess.
		*/

		template<class K>
		inline S improve(S s, F f, S v, K k) const
		{
			auto dvs = vega(f, s, k);
			S vc = value(f, s, payoff::call(k));
			return s - (vc - v) / dvs;
			//S sp = s - (vc - f + k - v) / dvs;

			//return abs(sc - s) < abs(sp - s) ? sc : sp;
		}

		// Vol matching put or call option value using Newton-Raphson.
		// There is no need to specify if the value is for a put or a call.
		template<class K>
		K implied(F f, S v, K k, S s = 0, size_t n = 0, S eps = 0) const
		{
			ensure(v > 0);
			static constexpr S epsilon = std::numeric_limits<S>::epsilon();

			if (s == 0) {
				F c = value(f, s, f);
				s = 2.5 * c / f;
			}
			if (k < 0) {
				k = -k;
				v = f - k + v;
			}
			if (n == 0) {
				n = 100;
			}
			if (eps == 0) {
				eps = sqrt(epsilon);
			}
			else if (eps <= epsilon) {
				eps = 10 * epsilon;
			}

			S s_ = s + 2 * eps;
			while (fabs(s_ - s) > eps) {
				s_ = improve(s, f, v, k);
				std::swap(s_, s);
				if (--n == 0) {
					break;
				}
			}
			ensure(n != 0);

			return s_;
		}
	};

}
