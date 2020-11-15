// fms_option.h - option valuation and greeks
// https://keithalewis.github.io/math/op.html
#pragma once
#include <algorithm>
#include <cmath>
#include <concepts>
#include <functional>
#include <limits>
#include "fms_ensure.h"

namespace fms {

	template<class M,
		class F = typename M::type, class S = typename M::ctype, class K = typename M::type,
		class X = std::common_type_t<F, S, K>>
	class option {
		const M& m;
	public:
		option(const M&m)
			: m(m)
		{ }
		~option()
		{ }

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
			X f_k = f - k;

			if (k < 0) { // put
				k = -k; 
				f_k = 0;
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
			return value<M>(f, s, k);
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

			return m.cdf(x, s, 1) / (f * s); //!!! only for normal model
		}

		// f - k = c - p so dc/ds = dp/ds
		X vega(F f, S s, K k) const
		{
			if (k < 0) {
				k = -k;
			}

			auto x = moneyness(f, s, k);

			return m.cdf(x, s, 1) * f * s; // /sigma //!!! only for normal model
		}

		// Volatility matching option price using Newton-Raphson.
		X implied(F f, S v, K k, S s0 = 0, size_t n = 0, S eps = 0) const
		{
			static S epsilon = std::numeric_limits<S>::epsilon();

			if (k < 0) { // put
				v = v - f + k;
				k = -k;
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
				s_ = s0 - (v - value(f, s0, k)) / vega(f, s0, k);
				s0 = s_;
				if (--n == 0) {
					break;
				}
			}
			ensure(n != 0);

			return s_;
		}
		X put_implied(F f, S p, K k, S s0 = S(0.1), size_t n = 100, S eps = 0) const
		{
			return implied(f, p, -k, s0, n, eps);
		}
		X call_implied(F f, S c, K k, S s0 = S(0.1), size_t n = 100, S eps = 0) const
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

	template<class M, class F = M::type, class S = M::type, class K = M::type>
	struct option_model : public option_nvi<F,S,K> {
		using X = std::common_type_t<F, S, K>;

		~option_model() = default;

		X moneyness_(F f, S s, K k) //override
		{
			return option<M, F, S, K>::moneyness(f, s, k);
		}
		X value_(F f, S s, K k) //override
		{
			return option<M, F, S, K>::value(f, s, k);
		}
		X delta_(F f, S s, K k) //override
		{
			return option<M, F, S, K>::delta(f, s, k);
		}
		X gamma_(F f, S s, K k) //override
		{
			return option<M, F, S, K>::gamma(f, s, k);
		}
		X vega_(F f, S s, K k) //override
		{
			return option<M, F, S, K>::vega(f, s, k);
		}
		X implied_(F f, S v, K k) //override
		{
			return option<M, F, S, K>::implied(f, v, k);
		}
	};
}
