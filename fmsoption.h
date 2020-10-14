// fmsoption.h - option pricing
#pragma once
#include <concepts>
#include <functional>
#include <limits>
#include "fmserror.h"

namespace fms {

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
}
