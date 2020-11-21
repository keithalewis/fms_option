// fms_variate_normal.h - normal distribution

#pragma once
#include <cmath>

namespace fms::variate {

	template<class X = double, class S = double>
	class normal {
		static constexpr X M_SQRT2 = X(1.41421356237309504880);
		static constexpr X M_SQRT2PI = X(2.50662827463100050240);
		X mu, sigma;
	public:
		typedef X type;
		typedef S ctype;

		normal(X mu = 0, X sigma = 1)
			: mu(mu), sigma(sigma == 0 ? 1 : sigma)
		{ }
		normal(const normal&) = default;
		normal& operator=(const normal&) = default;
		~normal()
		{ }

		// Normal mean 0 variance 1
		X cdf0(X x, size_t n = 0) const noexcept
		{
			if (n == 0) {
				return (1 + ::erf(x / X(M_SQRT2))) / 2;
			}

			X phi = ::exp(-x * x / X(2)) / X(M_SQRT2PI);

			return phi * H(n - 1, x) * (n % 2 == 0 ? -1 : 1) ;
		}

		X cdf(X x, S s = 0, size_t n = 0) const noexcept
		{
			return cdf0(((x - mu) / sigma) - s, n)/::pow(sigma,n);
		}

		// cumulant kappa(s) = mu s + sigma^2 s^2/2
		S cumulant(S s, size_t n = 0) const noexcept
		{
			if (n == 0) {
				return S(mu) * s + S(sigma) * S(sigma) * s * s/ 2;
			}
			if (n == 1) {
				return S(mu) + S(sigma) * s;
			}
			if (n == 2) {
				return S(sigma);
			}

			return S(0);
		}
	private:
		// Hermite polynomials H_0(x) = 1, H_1(x) = x, H_{n+1}(x) = x H_n(x) - n H_{n-1}(x)
		static constexpr X H(size_t n, X x) noexcept
		{
			if (n == 0) {
				return X(1);
			}
			if (n == 1) {
				return x;
			}

			return x * H(n - 1, x) - X(n - 1) * H(n - 2, x);
		}
	};


}