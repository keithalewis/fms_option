// fms_option.h - option pricing and greeks
// The underlying at expiration is F = f exp(s X - kappa(s)),
// where kappa(s) = log E[s X] is the cumulant of X.
// Note E[F] = f and Var(log F) = s^2 if E[X] = 0, E[X^2] = 1.
// The Black model is X standard normal and s = sigma sqrt(t).
#pragma once
#include <concepts>
#include <functional>
#include <limits>
#include <cassert>

#define ensure assert

namespace fms {

	// F <= k iff X <= (log(k/f) + kappa(s))/s
	template<class D, class F, class S, class K>
	inline auto moneyness(F f, S s, K k)
	{
		ensure (f > 0);
		ensure (s > 0);
		ensure (k > 0);

		return (log(k / f) + D::kappa(s)) / s;
	}

	template<class D, class F, class S, class K>
	inline auto put_value(F f, S s, K k)
	{
		auto z = moneyness<D>(f, s, k);

		return k * D::cdf(z) - f * D::cdf(z, s);
	}
	// Put-call parity: c - p = f - k
	template<class D, class F, class S, class K>
	inline auto call_value(F f, S s, K k)
	{
		return put_value<D>(f, s, k) + f - k;
	}

	// (d/df) E[max{k - F, 0}] = E[-exp(s X - kappa(s)) 1(F <= k)] = -P_s(F <= k)
	template<class D, class F, class S, class K>
	inline auto put_delta(F f, S s, K k)
	{
		auto z = moneyness<D>(f, s, k);

		return -D::cdf(z, s);
	}
	// c = p + f - k so dc/df = dp/df + 1
	template<class D, class F, class S, class K>
	inline auto call_delta(F f, S s, K k)
	{
		return put_delta<D>(f, s, k) + 1;
	}

	// c = p + f - k so d^2c/df^2 = d^2p/df^2
	template<class D, class F, class S, class K>
	inline auto option_gamma(F f, S s, K k)
	{
		auto z = moneyness<D>(f, s, k);

		return D::pdf(z, s)/(f*s);
	}

	// f - k = c - p so dc/ds = dp/ds
	template<class D, class F, class S, class K>
	inline auto option_vega(F f, S s, K k)
	{
		auto z = moneyness<D>(f, s, k);

		return D::pdf(z, s) / (f * s);
	}

	// Newton-Raphson
	// eps = sqrt(machine epsilon), better: use ULP
	// Volatility matching put price.
	template<class D, class F, class S, class K>
	inline auto option_implied(F f, S p, K k, S s0 = S(0.1), S eps = S(1e-8), size_t n = 100)
	{
		S s_ = s0 - (p - put_value(f, s0, k)) / option_vega(f, s0, k);

		while (n-- && fabs(s_ - s0) > eps) {
			s0 = s_;
			s_ = s0 - (p - put_value(f, s0, k)) / option_vega(f, s0, k);
		}
		//!!! warn if n == 0

		return s_;
	}
}
