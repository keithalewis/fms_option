// fms_logistic.h - Logistic distribution
// F(x;a) = 1/(1 + exp(-ax))
// kappa(s;a) = log E[exp(s X)] = \int_R exp(sx) dF(x)
// Let u = F(x) so exp(x) = u^a (1 - u)^{-a}
// \int_R exp(sx) dF(x) = \int_0^1 u^{sa} (1 - u)^{-sa} du = Beta(1 + sa, 1 - sa).
// Using Beta(alpha, beta) = Gamma(alpha) Gamma (beta)/Gamma(alpha + beta) we have
// kappa(s;a) = log Gamma(1 + sa) + log Gamma(1 - sa) since Gamma(2) = 1.
// kappa'(s;a) = a Gamma'(1 + sa)/Gamma(1 + sa) - a Gamma'(1 - sa)/Gamma(1 - sa)
// kappa''(s;a) = a {a G(1 + sa)G''(1 + sa) - a G'(1 + sa)^2}/G(1 + sa)^2
//               -a {-a G(1 - sa)G''(1 - sa) + a G'(1 - sa)^2}/G(1 - sa)^2.
// log Gamma(1 + z) = - gamma z + sum{k >= 2} zeta(k)/k (-z)^k.
#pragma once
#define _USE_MATH_DEFINES 
#include <cmath>

namespace fms {

	template<class X = double>
	struct logistic {
		// scale parameter for variance 1
		static constexpr X a = sqrt(X(3)) / X(M_PI);
		static X pdf(X x, X s = 0)
		{
			X e = exp(-x / a);

			return e / (a * (1 + e) * (1 + e));
		}
		// Use incomplete beta function.
		static X cdf(X x, X s = 0)
		{
			return 1 / (1 + exp(-x / a));
		}
		// cumulant
		static X kappa(X s)
		{
			return s * s / 2;
		}

	};

}
