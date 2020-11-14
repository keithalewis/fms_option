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
/// _gamma_ is $d^2v/df^2 = E[π''(F)\exp(s X - κ(s))^2]$, and 
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
#include <algorithm>
#include <cmath>
#include <concepts>
#include <functional>
#include <limits>
#include "fms_ensure.h"

namespace fms {

	template<class D,
		class F = D::type, class S = D::type, class K = D::type,
		class X = std::common_type_t<F, S, K>>
	struct option {

		// F <= k iff X <= (log(k/f) + kappa(s))/s
		static X moneyness(F f, S s, K k)
		{
			ensure (f > 0);
			ensure (s > 0);
			ensure (k > 0);

			return (::log(k / f) + D::kappa(s)) / s;
		}

		// p = k P(F <= k) - f P^s(F <= k)
		// call (k > 0) or put (k < 0) value
		// c = p + f - k
		static X value(F f, S s, K k)
		{
			X f_k = f - k;

			if (k < 0) { // put
				k = -k; 
				f_k = 0;
			}

			if (f == 0) {
				return X(0);
			}
			if (s == 0) { // intrinsic
				return (std::max)(f_k == 0 ? -f_k : f_k, X(0));
			}
			if (k == 0) {
				return f_k == 0 ? 0 : f;
			}

			X x = moneyness(f, s, k);

			return k * D::cdf(x) - f * D::cdf(x, s) + f_k;
		}
		static X put_value(F f, S s, K k)
		{
			return value(f, s, -k);
		}
		static X call_value(F f, S s, K k)
		{
			return value<D>(f, s, k);
		}

		// (dp/df) = -P_s(F <= k)
		// call (k > 0) or put (k < 0) value
		// dc = dp + 1
		static X delta(F f, S s, K k)
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

			return -D::cdf(x, s) + one;
		}
		static X put_delta(F f, S s, K k)
		{
			return delta(f, s, -k);
		}
		static X call_delta(F f, S s, K k)
		{
			return delta(f, s, k);
		}

		// c = p + f - k so d^2c/df^2 = d^2p/df^2
		static X gamma(F f, S s, K k)
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

			return D::pdf(x, s) / (f * s); //!!! only for normal model
		}

		// f - k = c - p so dc/ds = dp/ds
		static X vega(F f, S s, K k)
		{
			if (k < 0) {
				k = -k;
			}

			auto x = moneyness(f, s, k);

			return D::pdf(x, s) * f * s; // /sigma //!!! only for normal model
		}

		// Volatility matching option price using Newton-Raphson.
		static X implied(F f, S v, K k, S s0 = S(0.1), size_t n = 10, S eps = 0)
		{
			static S epsilon = std::numeric_limits<S>::epsilon();

			if (k < 0) { // put
				v = v - f + k;
				k = -k;
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
				s_ = s0 - (v - value(f, s0, k)) / vega(f, s0, k);
				s0 = s_;
				if (--n == 0) {
					break;
				}
			}
			ensure(n != 0);

			return s_;
		}
		static X put_implied(F f, S p, K k, S s0 = S(0.1), size_t n = 100, S eps = 0)
		{
			return implied(f, p, -k, s0, n, eps);
		}
		static X call_implied(F f, S c, K k, S s0 = S(0.1), size_t n = 100, S eps = 0)
		{
			return implied(f, c, k, s0, n, eps);
		}

	};

	// Parameterize by model for add-ins.
	template<class F = double, class S = double, class K = double>
	struct option_nvi {
		using X = std::common_type_t<F, S, K>;

		virtual ~option_nvi() = default;

		X moneyness(F f, S s, K k)
		{
			return moneyness_(f, s, k);
		}
		X value(F f, S s, K k)
		{
			return value_(f, s, k);
		}
		X put_value(F f, S s, K k)
		{
			return value_(f, s, -k);
		}
		X call_value(F f, S s, K k)
		{
			return value_(f, s, k);
		}
		X delta(F f, S s, K k)
		{
			return delta_(f, s, k);
		}
		X put_delta(F f, S s, K k)
		{
			return delta_(f, s, -k);
		}
		X call_delta(F f, S s, K k)
		{
			return delta_(f, s, k);
		}
		X gamma(F f, S s, K k)
		{
			return gamma_(f, s, k);
		}
		X vega(F f, S s, K k)
		{
			return vega_(f, s, k);
		}
		X implied(F f, S s, K k)
		{
			return implied_(f, s, k);
		}
		X put_implied(F f, S p, K k)
		{
			return implied_(f, p, -k);
		}
		X call_implied(F f, S c, K k)
		{
			return implied_(f, c, k);
		}

	private:
		
		virtual X moneyness_(F f, S s, K k) = 0;
		virtual X value_(F f, S s, K k) = 0;
		virtual X delta_(F f, S s, K k) = 0;
		virtual X gamma_(F f, S s, K k) = 0;
		virtual X vega_(F f, S s, K k) = 0;
		virtual X implied_(F f, S v, K k) = 0;
	};

	template<class D, class F = D::type, class S = D::type, class K = D::type>
	struct option_model : public option_nvi<F,S,K> {
		using X = std::common_type_t<F, S, K>;

		~option_model() = default;

		X moneyness_(F f, S s, K k) //override
		{
			return option<D, F, S, K>::moneyness(f, s, k);
		}
		X value_(F f, S s, K k) //override
		{
			return option<D, F, S, K>::value(f, s, k);
		}
		X delta_(F f, S s, K k) //override
		{
			return option<D, F, S, K>::delta(f, s, k);
		}
		X gamma_(F f, S s, K k) //override
		{
			return option<D, F, S, K>::gamma(f, s, k);
		}
		X vega_(F f, S s, K k) //override
		{
			return option<D, F, S, K>::vega(f, s, k);
		}
		X implied_(F f, S v, K k) //override
		{
			return option<D, F, S, K>::implied(f, v, k);
		}
	};
}
