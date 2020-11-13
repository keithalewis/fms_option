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

namespace fms::variate {

	template<class X = double>
	struct normal {
		typedef X type;

		static constexpr X M_SQRT2 = X(1.41421356237309504880);
		static constexpr X M_SQRT2PI = X(6.28318530717958647688);

		static X pdf(X x, X s = 0)
		{
			X x_ = x - s;

			return ::exp(-x_ * x_ / X(2)) / X(M_SQRT2PI);
		}
		static X cdf(X x, X s = 0)
		{
			X x_ = x - s;

			return (1 + ::erf(x_ / X(M_SQRT2))) / 2;
		}
		// cumulant
		static X kappa(X s)
		{
			return s * s / 2;
		}
	};
}