// fms_variate_discrete.h - discrete random variate
#pragma once
#include <cmath>
#include <algorithm>
#include <compare>
#include <numeric>
#include <valarray>
#include "fms_ensure.h"

namespace fms::variate {

	template<class X = double, class S = X>
	class discrete {
		std::valarray<X> x;
		std::valarray<X> p;
	public:
		typedef X type;
		typedef S ctype;

		// zero
		discrete()
			: x({ 0 }), p({1})
		{ }
		discrete(size_t n, const X* _x, const X* _p)
			: x(_x, n), p(_p, n)
		{
			if (n == 1) {
				p[0] = 1;
			}
			ensure(0 <= p.min());
			ensure(fabs(p.sum() - X(1)) <= std::numeric_limits<X>::epsilon());
		}
		discrete(const std::initializer_list<X>& x, const std::initializer_list<X>& p)
			: discrete(x.size(), x.begin(), p.begin())
		{
			ensure(x.size() == p.size());
		}
		discrete(const discrete&) = default;
		discrete& operator=(const discrete&) = default;
		~discrete()
		{ }

		//auto operator<=>(const discrete&) const = default;

		X cdf(X x_, S s = 0, size_t n = 0) const noexcept
		{
			if (n == 0) {
				X P = 0;
				S ks = cumulant(s);

				// return sum(exp(s*x[x <= x_] - cumulant(s))*p[x <= x_]);
				for (size_t i = 0; i < x.size(); ++i) {
					P += (x[i] <= x_) * ::exp(s*x[i] - ks)*p[i];
				}

				return P;
			}

			// return infinity at point masses
			return std::end(x) == std::find(std::begin(x), std::end(x), x_) ? X(0) : std::numeric_limits<X>::infinity();
		}
		S cumulant(S s, size_t n = 0) const noexcept
		{
			S e0 = e(s, 0);
			if (n == 0) {
				return ::log(e0);
			}

			S e1 = e(s, 1);
			if (n == 1) {
				return e1 / e0;
			}

			S e2 = e(s, 2);
			if (n == 2) {
				return (e0 * e2 - e1 * e1) / (e0 * e0);
			}

			return std::numeric_limits<S>::quiet_NaN();
		}
	private:
		// (d/ds)^n sum_i exp(s x_i) p_i = sum_i exp(s x_i) x_i^n p_i
		S e(S s, size_t n) const
		{
			S E = 0;

			for (size_t i = 0; i < x.size(); ++i) {
				E += ::exp(s * S(x[i])) * ::pow(S(x[i]), S(n)) * S(p[i]);
			}

			return E;
		}
	};

}

