// fms_normal.h - normal distribution
/// # Normal
///
/// The density function of a standard normal variate $X$ is $\phi(x) = \exp(-x^2/2)/\sqrt{2π}$.
/// and its cumulant is $κ(s) = s^2/2$.
/// For any normal variate $N$, $E[\exp(N)] = \exp(E[N] + \Var(N)/2)$
/// and $E[\exp(N) g(N)] = \exp(E[N]) E[g(N + \Var(N))]$.
/// If $N$ and $N_j$ are jointly normal $[\exp(N) g(N_j)] = \exp(E[N]) E[g(N_j + \Cov(N, N_j))]$.
/// If $\Phi$ is the cumulative distribution of $X$ then 
/// $\Phi^s(x) = P^s(X \le x) = E[\exp(s X - s^2/2) 1(X \le x)] = P(X + s \le x) = \Phi(x - s)$ so
/// $P^s(N \le x) = P^s(X \le (x - μ)/σ) = \Phi(z - s)$ where $z = (x - μ)/σ$.
/// The $n$-th derivative is $(d/dx)^n P^s(N \le x) = \Phi^{(n)}(z - s)/\σ^n$.
/// 
/// Let $\phi(x) = \Phi'(x) = \exp(-x^2)/sqrt{2 π}$.
/// We have $\phi^(n)(x) = (-1)^n \phi(x) H_n(x)$ where $H_n(x)$ is the Hermite polynomial of degree $n$.
/// They satisfy $H_0(x) = 1$, $H_1(x) = x$, and $H_{n+1}(x) = x H_n(x) - n H_{n-1}(x)$, $n > 1$.

#pragma once
#include <cmath>

namespace fms::variate {

	template<class X = double, class S = double>
	class normal {
		static constexpr X M_SQRT2 = X(1.41421356237309504880);
		static constexpr X M_SQRT2PI = X(6.28318530717958647688);
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

		X cdf(X x, S s = 0, size_t n = 0) const
		{
			X z = (x - mu) / sigma;
			X z_ = z - X(s);

			if (n == 0) {
				return (1 + ::erf(z_ / X(M_SQRT2))) / 2;
			}

			X phi = ::exp(-z_ * z_ / X(2)) / (sigma*X(M_SQRT2PI));
			
			return phi * H(n - 1, z_) / ::pow(-sigma, X(n - 1));
		}
		// cumulant
		S cumulant(S s, size_t n = 0) const
		{
			if (n == 0) {
				return S(mu) * s + S(sigma) * S(sigma) * s * s/ 2;
			}
			if (n == 1) {
				return S(mu) + S(sigma) * s;
			}

			return S(0);
		}
	private:
		// Hermite polynomials
		static constexpr X H(size_t n, X x)
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