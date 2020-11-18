// fms_variate_discrete.h - discrete random variate
#pragma once
#include <algorithm>
#include <cmath>
#include <compare>
#include <vector>
#include <numeric>
#include "fms_ensure.h"

namespace fms::variate {

	template<class X = double, class S = X>
	class discrete {
		std::vector<X> x;
		std::vector<X> p;
	public:
		// zero
		discrete()
			: x({ 0 }), p({1})
		{ }
		discrete(size_t n, const X* x, const X* p)
			: x(x, x + n), p(p, p + n)
		{
			ensure(0 <= std::min({ p, p + n }));
			auto psum = std::accumulate(p, p + n, X(0));
			// sum_i p_i = 1
			ensure(fabs(psum - X(1)) <= std::numeric_limits<X>::epsilon());
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

				for (size_t i = 0; i < x.size(); ++i) {
					P += (x[i] <= x_) * ::exp(s*x[i] - ks)*p[i];
				}

				return P;
			}

			// return infinity at point masses
			return x.end() == std::find(x.begin(), x.end(), x_) ? X(0) : std::numeric_limits<X>::infinity();
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
				E += ::exp(s * S(x[i]) * ::pow(S(x[i]), S(n)) * S(p[i]));
			}

			return E;
		}
	};

}

