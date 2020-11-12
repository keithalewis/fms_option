// fms_option.h - option valuation and greeks
/// # FMS Option
///
/// European _option valuation_ involves calculating the expected value of
/// the _option payoff_ at _expiration_. Greeks are derivatives of the value.
/// Payoff is a function of the _underlying_ at expiration.
/// The underlying at expiration is $F = f \exp(s X - κ(s))$,
/// where $κ(s) = \log E[s X]$ is the cumulant of $X$.
/// Note $E[F] = f$ and $\Var(\log F) = s^2$ if $E[X] = 0$ and $E[X^2] = 1$.
/// For example, the Black model is $X$ standard normal and $s = σ \sqrt(t)$
/// where $σ$ is the volatilty and $t$ is time in years to expiration.
/// 
/// The (forward) value of an option paying $π(F)$ at expiration is $v = E[π(F)]$ so 
/// _delta_ is $dv/df = E[π'(F) dF/df] = E[π'(F)\exp(s X - κ(s)) ]$, 
/// _gamma_ is $d^2v/df^2 = E[π''(F)\exp(s X - κ(s))^2 ]$, and 
/// _vega_ is $dv/ds = E[π'(F) dF/ds] = E[π'(F)F(X - κ'(s)]$.
/// 
/// A _put option_ pays $\max\{k - F,0\}$ and a _call option_ pays $\max\{F - k, 0\}$ at expiration.
/// Note $max\{F - k, 0\} - \max\{k - F,0\} = F - k$ is a _forward_ with _strike_ $k$.
/// Define _moneyness_ $x$ by $F = k$ iff $X = x = (\log(k/F) + κ(s))/s$.
/// The value of a put is 
/// $p = E[\max\{k - F, 0\}] = E[(k - F)1(F\le k)] = k P(X \le x) - f P^s(X \le x)$
/// where $dP^s/dP = \exp(s X - κ(s))$.
/// Put delta is $dp/df = E[1(F \le k)exp(s X - κ(s))] = P^s(X \le d)$ and
/// put gamma is $d^2p/df^2 = E^s[δ_k(F)] = ...$
#pragma once
#include "fms_ensure.h"
#include <concepts>
#include <functional>
#include <limits>

namespace fms {

	// F <= k iff X <= (log(k/f) + kappa(s))/s
	template<class D, class F = double, class S = double, class K = double>
	inline auto moneyness(F f, S s, K k)
	{
		ensure (f > 0);
		ensure (s > 0);
		ensure (k > 0);

		return (log(k / f) + D::kappa(s)) / s;
	}

	template<class D, class F = double, class S = double, class K = double>
	inline auto put_value(F f, S s, K k)
	{
		using X = std::common_type_t<F, S, K>;

		if (f == 0 or k == 0) {
			return X(0);
		}
		if (s == 0) {
			return X(k - f) * X(f <= k);
		}

		auto x = moneyness<D>(f, s, k);

		return k * D::cdf(x) - f * D::cdf(x, s);
	}
	// Put-call parity: c - p = f - k
	template<class D, class F = double, class S = double, class K = double>
	inline auto call_value(F f, S s, K k)
	{
		return put_value<D>(f, s, k) + f - k;
	}

	// (d/df) E[max{k - F, 0}] = E[-exp(s X - kappa(s)) 1(F <= k)] = -P_s(F <= k)
	template<class D, class F = double, class S = double, class K = double>
	inline auto put_delta(F f, S s, K k)
	{
		using X = std::common_type_t<F, S, K>;

		if (f == 0 or k == 0) {
			return X(0);
		}
		if (s == 0) {
			return X(f <= k);
		}

		X x = moneyness<D>(f, s, k);

		return -D::cdf(x, s);
	}
	// c = p + f - k so dc/df = dp/df + 1
	template<class D, class F = double, class S = double, class K = double>
	inline auto call_delta(F f, S s, K k)
	{
		return put_delta<D>(f, s, k) + 1;
	}

	// c = p + f - k so d^2c/df^2 = d^2p/df^2
	template<class D, class F = double, class S = double, class K = double>
	inline auto option_gamma(F f, S s, K k)
	{
		using X = std::common_type_t<F, S, K>;

		if (f == 0 or k == 0) {
			return X(0);
		}
		if (s == 0) {
			return f == k ? std::numeric_limits<X>::infinity() : X(0);
		}

		X x = moneyness<D>(f, s, k);

		return D::pdf(x, s) / (f * s); //!!! only for normal model
	}

	// f - k = c - p so dc/ds = dp/ds
	template<class D, class F = double, class S = double, class K = double>
	inline auto option_vega(F f, S s, K k)
	{
		auto x = moneyness<D>(f, s, k);

		return D::pdf(x, s) * f * s; // /sigma //!!! only for normal model
	}

	// Newton-Raphson
	// eps = sqrt(machine epsilon), better: use ULP
	// Volatility matching put price.
	template<class D, class F = double, class S = double, class K = double>
	inline auto put_implied(F f, S p, K k, S s0 = S(0.1), S eps = S(1e-8), size_t n = 100)
	{
		S s_ = s0 - (p - put_value(f, s0, k)) / option_vega(f, s0, k);

		while (n-- and fabs(s_ - s0) > eps) {
			s0 = s_;
			s_ = s0 - (p - put_value(f, s0, k)) / option_vega(f, s0, k);
		}
		//!!! warn if n == 0

		return s_;
	}
	template<class D, class F = double, class S = double, class K = double>
	inline auto call_implied(F f, S c, K k, S s0 = S(0.1), S eps = S(1e-8), size_t n = 100)
	{
		return put_implied(f, c - f + k, k, s0, eps, n);
	}

	// Parameterize by model for add-ins.
	template<class D, class F = double, class S = double, class K = double>
	struct option {
		static auto moneyness(F f, S s, K k)
		{
			return fms::moneyness<D>(f, s, k);
		}
		static auto put_value(F f, S s, K k)
		{
			return fms::put_value<D>(f, s, k);
		}
		static auto call_value(F f, S s, K k)
		{
			return fms::call_value<D>(f, s, k);
		}
		static auto value(F f, S s, K k)
		{
			return k > 0 ? call_value(f, s, k) : put_value(f, s, -k);
		}
		static auto put_delta(F f, S s, K k)
		{
			return fms::put_delta<D>(f, s, k);
		}
		static auto call_delta(F f, S s, K k)
		{
			return fms::call_delta<D>(f, s, k);
		}
		static auto delta(F f, S s, K k)
		{
			return k > 0 ? call_delta(f, s, k) : put_delta(f, s, -k);
		}
		static auto gamma(F f, S s, K k)
		{
			return fms::option_gamma<D>(f, s, k);
		}
		static auto vega(F f, S s, K k)
		{
			return fms::option_vega<D>(f, s, k);
		}
		static auto put_implied(F f, S p, K k)
		{
			return fms::put_implied<D>(f, p, k);
		}
		static auto call_implied(F f, S c, K k)
		{
			return fms::call_implied<D>(f, c, k);
		}
		static auto implied(F f, S s, K k)
		{
			return k > 0 ? call_implied(f, s, k) : put_implied(f, s, -k);
		}
		// useful for valuing portfolios
		template<class X = double>
		static void greeks(X f, X s, X k, X* v, X* d = 0, X* g = 0, X* e = 0)
		{
			if (v) {
				*v += value(f, s, k);
			}
			if (d) {
				*d += delta(f, s, k);
			}
			if (g) {
				*g += gamma(f, s, k);
			}
			if (e) {
				*e += vega(f, s, k);
			}
		}
	};
}
