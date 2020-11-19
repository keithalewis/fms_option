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

		X cdf(X x, S s = 0, size_t n = 0) const noexcept
		{
			X z = (x - mu) / sigma;
			X z_ = z - X(s);

			if (n == 0) {
				return (1 + ::erf(z_ / X(M_SQRT2))) / 2;
			}

			X phi = ::exp(-z_ * z_ / X(2)) / (sigma*X(M_SQRT2PI));
			
			return phi * H(n - 1, z_) / ::pow(-sigma, X(n - 1));
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