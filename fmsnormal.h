// fmsnormal.h - normal distribution
#pragma once
#define _USE_MATH_DEFINES 
#include <cmath>
#include "fmserror.h"

// sqrt(2 pi)
#define M_SQRT2PI (2 * M_SQRT2 / M_2_SQRTPI)

namespace fms {

	template<class X = double>
	struct normal {
		static X pdf(X x, X s = 0)
		{
			X x_ = x - s;

			return exp(-x_ * x_ / 2) / M_SQRT2PI;
		}
		static X cdf(X x, X s = 0)
		{
			X x_ = x - s;

			return (X(1) + ::erf(x_ / M_SQRT2)) / X(2);
		}
		// cumulant
		static X kappa(X s)
		{
			return s * s / 2;
		}
	};
}