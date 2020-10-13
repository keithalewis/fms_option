// fmsoption.h - option pricing
#pragma once
#include <concepts>
#include <functional>
#include <limits>
#include "fmserror.h"

namespace fms {

	// Assume E[X] = 0, E[X^2] = 1.
	// The cumulant is kappa(s) = log E[exp(s X)] = s^2/2 + ...
	// Define F = f exp(s X - kappa(s)).
	// Note E[F] = f, Var(log F) = s^2, and F <= k iff X <= (log(k/f) + kappa(s))/s.
	// E[max{k - F, 0}] = E[(k - F) 1(F <= k)] = k P(F <= k) - f P_s(F <= k),
	// where dP_s/dP = exp(s X - kappa(s)).

	using X = double;

	// Dependency injection
	template<class D, class T = X>
	concept Distribution = requires (T x, T s) {
		D::pdf(x, s) -> T;  // Esscher transform of density
		D::cdf(x, s) -> T;  // Esscher transform of distribution
		D::kappa(x) -> T; // cumulant
	};

	//template<Distribution D>
	template<class D>
	inline X moneyness(X f, X s, X k)
	{
		FLOAT_ENSURE(f > 0);

		return (log(k / f) + D::kappa(s)) / s;
	}

	//template<Distribution D>
	template<class D>
	inline X put_value(X f, X s, X k)
	{
		X z = moneyness<D>(f, s, k);

		return k * normal::cdf(z, 0) - f * normal::cdf(z, s);
	}
	// Put-call parity: c - p = f - k
	template<Distribution D>
	inline X call_value(X f, X s, X k)
	{
		return put_value<D>(f, s, k) + f - k;
	}

	// (d/df) E[max{k - F, 0}] = E[-exp(s X - kappa(s)) 1(F <= k)] = -P_s(F <= k)
	template<Distribution D>
	inline X put_delta(X f, X s, X k)
	{
		X z = moneyness<D>(f, s, k);

		return -D::cdf(z, s);
	}
	// c = p + f - k so dc/df = dp/df + 1
	template<Distribution D>
	inline X call_delta(X f, X s, X k)
	{
		return put_delta<D>(f, s, k) + X(1);
	}

	// (d^2/df^2) E[max{k - F, 0}] = (d/df) -P_s(F <= k) = -P_s'(k)
	template<Distribution D>
	inline X put_gamma(X f, X s, X k)
	{
		X z = moneyness<D>(f, s, k);

		return -D::pdf(z, s);
	}
	// c = p + f - k so dc/df = dp/df + 1
	template<Distribution D>
	inline X call_gamma(X f, X s, X k)
	{
		return put_gamma<D>(f, s, k);
	}
}
