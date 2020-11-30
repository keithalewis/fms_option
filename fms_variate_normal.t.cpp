// fms_option.t.cpp - test fms option
#include <cassert>
#include <functional>
#include <iostream>
#include <utility>
#include "fms_test.h"
#include "fms_variate_normal.h"
#include "fms_variate_handle.h"

using namespace fms;
using namespace fms::variate;

template<class X>
int test_variate_normal()
{
	X eps = std::numeric_limits<X>::epsilon();
	X dx = X(0.001);

	{
		variate::normal<X, X> n;

		assert(n.cumulant(0) == 0); // true for all cumulants
		assert(n.cumulant(0, 1) == 0); // mean
		assert(n.cumulant(0, 2) == 1); // variance
		assert(n.cumulant(0, 3) == 0);

		test_variate(n, dx);
	}
	{
		X mu = 2;
		X sigma = 3;
		variate::normal<X, X> n(2,3);

		assert(n.cumulant(0) == 0); // true for all cumulants
		assert(n.cumulant(0, 1) == mu); // mean
		assert(n.cumulant(0, 2) == sigma*sigma); // variance
		assert(n.cumulant(0, 3) == 0);

		test_variate(n, dx);
	}
	{
		X mu = 2;
		X sigma = 3;
		variate::normal<X, X> nms(mu, sigma);
		variate_standard n(nms);

		assert(n.cdf(0) == X(0.5));

		assert(n.cumulant(0) == 0); // true for all cumulants
		assert(n.cumulant(0, 1) == 0); // mean
		assert(fabs(n.cumulant(0, 2) - 1) <= eps); // variance
		assert(n.cumulant(0, 3) == 0);
	}
	{
		X mu = 2;
		X sigma = 3;
		variate::normal<X> nms(mu, sigma);
		variate_handle n(nms);

		assert(n.cdf(mu) == X(0.5));

		assert(n.cumulant(0) == 0); // true for all cumulants
		assert(n.cumulant(0, 1) == mu); // mean
		assert(fabs(n.cumulant(0, 2) - sigma*sigma) <= eps); // variance
		assert(n.cumulant(0, 3) == 0);
	}
	
	return 0;
}
int test_variate_normal_f = test_variate_normal<float>();
int test_variate_normal_d = test_variate_normal<double>();
