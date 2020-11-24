// fms_variate_logistic.t.cpp - test fms option
#include <cassert>
#include <functional>
#include <iostream>
#include <utility>
#include "fms_test.h"
#include "fms_variate_logistic.h"
#include "fms_variate.h"

using namespace fms;
using namespace fms::variate;

template<class X>
int test_variate_logistic()
{
	X dx = X(0.001);

	{
		variate::logistic<X> n;

		assert(n.cumulant(0) == 0); // true for all cumulants
		assert(n.cumulant(0, 1) == 0); // mean
		assert(n.cumulant(0, 2) == 1); // variance
		assert(n.cumulant(0, 3) == 0);

		test_variate(n, dx);
	}

	return 0;
}
//int test_variate_logistic_f = test_variate_logistic<float>();
int test_variate_logistic_d = test_variate_logistic<double>();
