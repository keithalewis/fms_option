// fms_normal.h - normal distribution
/// # Normal
///
/// The density function of a standard normal variate $X$ is $\phi(x) = \exp(-x^2/2)/\sqrt{2π}$.
/// and its cumulant is $κ(s) = s^2/2$.
/// For any normal variate $N$, $E[\exp(N)] = \exp(E[N] + \Var(N)/2)$
/// and $E[\exp(N) g(N)] = \exp(E[N]) E[g(N + \Var(N))]$.
/// If $N$ and $M_j$ are jointly normal $[\exp(N) g(M_j)] = \exp(E[N]) E[g(M_j + \Cov(N, M_j))]$.
/// If $\Phi$ is the cumulative distribution of $X$ then 
/// $\Phi^s(x) = P^s(X \le x) = E[\exp(s X - s^2/2) 1(X \le x)] = P(X + s \le x) = \Phi(x - s)$.
#pragma once
#define _USE_MATH_DEFINES 
#include <cmath>
#include "fms_variate.h"

namespace fms::variate {

	template<class X = double>
	class normal : public variate_base<X> {
		static constexpr X M_SQRT2 = X(1.41421356237309504880);
		static constexpr X M_SQRT2PI = X(6.28318530717958647688);
		X mu, sigma;
	public:
		typedef X type;

		normal(X mu = 0, X sigma = 1)
			: mu(mu), sigma(sigma == 0 ? 1 : sigma)
		{ }
		normal(const normal&) = default;
		normal& operator=(const normal&) = default;
		~normal()
		{ }

		X cdf(X x, X s = 0, size_t n = 0) const
		{
			X x_ = x - s;
			X z = (x_ - mu) / sigma;

			if (n == 0) {
				return (1 + ::erf(z / X(M_SQRT2))) / 2;
			}

			X phi = ::exp(-z * z / X(2)) / (sigma*X(M_SQRT2PI));
			
			return phi; // * hermite(z,n) ...
		}
		// cumulant
		X cumulant(X s, size_t n = 0) const
		{
			if (n == 0) {
				return mu * s + sigma * sigma * s * s/ 2;
			}
			if (n == 1) {
				return mu + sigma * s;
			}

			return 0;
		}
	};
}