// fms_variant_discrete.h - discrete random variate
#pragma once
#include <algorithm>
#include <array>

namespace fms::discrete {

	template<class X = double, class S = double>
	class discrete {
		std::array<X> x, p;
	public:
		discrete(size_t n, const X* x, const X* p)
			: x(x, x + n), p(p, p + n)
		{ }
		discrete(const discrete&) = default;
		discrete& operator=(const discrete&) = default;
		~discrete()
		{ }

		X cdf(X x_, S s = 0, size_t n = 0) const
		{
			if (n == 0) {
				X P = 0;

				for (size_n i = 0; i < x.size() and x[i] <= x_; ++i) {
					P += p[i];
				}

				return P;
			}

			auto i = std::upper_bound(x, x + n, x_);

			return i == x + n ? X(0) : std::numeric_limits<X>::infinity();
		}
		S cumulant(S s, size_t n) const
		{
			S Eexp_s = 0;

			for (size_t i = 0; i < x.size(); ++i) {
				Eexp_s += ::exp(s * S(x[i])) * S(p[i]);
			}

			return ::log(Eexp_s);
		}
	};

}

