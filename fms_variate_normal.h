// fms_variate_normal.h - normal distribution
#pragma once
#include <cmath>

namespace fms::variate {

	template<class X = double, class S = X>
	class normal_impl
	{
		static constexpr X M_SQRT2 = X(1.41421356237309504880);
		static constexpr X M_SQRT2PI = X(2.50662827463100050240);
		X mu, sigma;
	public:
		typedef X xtype;
		typedef S stype;

		normal_impl(X mu = 0, X sigma = 1)
			: mu(mu), sigma(sigma == 0 ? 1 : sigma)
		{ }
		normal_impl(const normal_impl&) = default;
		normal_impl& operator=(const normal_impl&) = default;
		normal_impl(normal_impl&&) = default;
		normal_impl& operator=(normal_impl&&) = default;
		~normal_impl()
		{ }

		// Normal mean 0 variance 1
		static X cdf01(X x, size_t n = 0) noexcept
		{
			if (n == 0) {
				return (1 + ::erf(x / X(M_SQRT2))) / 2;
			}

			X phi = ::exp(-x * x / X(2)) / X(M_SQRT2PI);

			return phi * H(n - 1, x) * (n % 2 == 0 ? -1 : 1);
		}

		X cdf(X x, S s = 0, size_t n = 0) const noexcept
		{
			return cdf01(((x - mu) / sigma) - s, n)/::pow(sigma,X(n));
		}

		static S cumulant01(S s, size_t n = 0)
		{
			if (n == 0) {
				return s * s / 2;
			}
			if (n == 1) {
				return s;
			}
			if (n == 2) {
				return 1;
			}

			return S(0);
		}

		// cumulant kappa(s) = mu s + sigma^2 s^2/2
		S cumulant(S s, size_t n = 0) const noexcept
		{
			S m_ = S(mu);
			S s_ = S(sigma);

			if (n == 0) {
				return m_ * s + cumulant01(s_ * s, 0);
			}
			if (n == 1) {
				return m_ + cumulant01(s_ * s, 1) * s_;
			}
			if (n == 2) {
				return s_ * s_;
			}

			return 0;
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

	template<class X = double, class S = X>
	using normal = variate_model<normal_impl<X,S>>;
}