// fms_variate_discrete.t.cpp - test discrete variate
#include <cassert>
#include "fms_variate_discrete.h"

using namespace fms;

template<class X = double>
int test_variate_discrete()
{
	{
		variate::discrete<X, X> x;
		variate::discrete<X, X> x2(x);
		/*
		assert(x == x2);
		assert(!(x != x2));
		assert(!(x < x2));
		assert(x >= x2);
		*/
		x = x2;

		assert(x.cdf(-1) == 0);
		assert(x.cdf(0) == 1);
		assert(x.cdf(1) == 1);

		assert(x.cumulant(0) == 0);
	}
	{
		variate::discrete<X, X> x({ -1,1 }, { 0.5, 0.5 });
		variate::discrete<X, X> x2(x);
		/*
		assert(x == x2);
		assert(!(x != x2));
		assert(!(x < x2));
		assert(x >= x2);
		*/
		x = x2;

		assert(x.cdf(-2) == 0);
		assert(x.cdf(-1) == 0.5);
		assert(x.cdf(0) == 0.5);
		assert(x.cdf(1) == 1);
		assert(x.cdf(2) == 1);

		assert(x.cumulant(0, 0) == 0);
		assert(x.cumulant(0, 1) == 0);
		assert(x.cumulant(0, 2) == 1);

		for (X s : {X(-1), X(0), X(0.1), X(1)}) {
			X err = x.cumulant(s);
			err -= ::log(::cosh(s));
			assert(fabs(err) <= std::numeric_limits<X>::epsilon());
		}

	}


	return 0;
}
int test_variate_discrete_d = test_variate_discrete<double>();
int test_variate_discrete_f = test_variate_discrete<float>();
