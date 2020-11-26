// fms_logistic.h - Logistic distribution
// F(x;a) = 1/(1 + exp(-ax))
// kappa(s;a) = log E[exp(s X)] = \int_R exp(sx) dF(x)
// Let u = F(x) so exp(x) = u^a (1 - u)^{-a}
// \int_R exp(sx) dF(x) = \int_0^1 u^{sa} (1 - u)^{-sa} du = Beta(1 + sa, 1 - sa).
// Using Beta(alpha, beta) = Gamma(alpha) Gamma (beta)/Gamma(alpha + beta) we have
// kappa(s;a) = log Gamma(1 + sa) + log Gamma(1 - sa) since Gamma(2) = 1.
// log Gamma(1 + z) = - gamma z + sum{k >= 2} zeta(k)/k (-z)^k.
#pragma once
#define _USE_MATH_DEFINES 
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include "fms_ensure.h"

namespace fms::variate {

	template<class X = double, class S = X>
	struct logistic {
		// scale parameter for variance 1
		static constexpr X a = X(M_SQRT3) / X(M_PI);

		typedef X xtype;
		typedef S stype;

		// Use incomplete beta function.
		static X cdf(X x, S s = 0, size_t n = 0)
		{
			ensure(-1 < s and s < 1);

			if (n == 0) {
				X u = 1 / (1 + ::exp(-x / a));
				if (s == 0) {
					return u;
				}
				else {
					return gsl_sf_beta_inc(1 + s, 1 - s, u);
				}
			}

			X e = ::exp(-x / a);
			X du = (e / a) / ((1 + e) * (1 + e));
			if (n == 1) {
				if (s == 0) {
					return du;
				}
				else {
					return::exp(s * x - cumulant(s)) * du;
				}
			}
			
			return std::numeric_limits<X>::quiet_NaN();
	
		}
		// cumulant
		static S cumulant(S s, size_t n = 0)
		{
			ensure(-1 < s and s < 1);

			if (n == 0) {
				return gsl_sf_lngamma(1 + s) + gsl_sf_lngamma(1 - s);
			}

			int n_ = static_cast<int>(n - 1);

			return gsl_sf_psi_n(n_, 1 + s) - gsl_sf_psi_n(n_, 1 - s);
		}

	};

}
