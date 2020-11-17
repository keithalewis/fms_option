// fms_option.t.cpp - test fms option
#include <cassert>
#include <functional>
#include <iostream>
#include <utility>
#include "fms_test.h"
#include "fms_option.h"
#include "fms_variate_normal.h"
#include "fms_variate.h"

using namespace fms;
using namespace fms::variate;

template<class X>
int test_option_normal_value()
{
	X eps = std::numeric_limits<X>::epsilon();
	X f = X(100);
	X s = X(0.1); // 3-month 20% volatility
	X k = X(100);
	X p = X(3.9877611676744920);

	X x;
	{
		variate::normal<X,X> n;
		option m(n);
		x = m.moneyness(f, s, k);
		x -= X(0.05);
		ensure (fabs(x) < eps);
		x = m.put_value(f, s, k);
		x -= p;
		ensure(fabs(x) <= 10 * eps);
		X cp = m.call_value(f, s, k) - m.put_value(f, s, k);
		ensure(cp == f - k);
	}
	
	{
		variate_model n(variate::normal<X, X>{});
		option m(n);
		x = m.moneyness(f, s, k);
		x -= X(0.05);
		ensure(fabs(x) < eps);
		x = m.put_value(f, s, k);
		x -= p;
		ensure(fabs(x) <= 10 * eps);
	}
	

	x = x;

	return 0;
}
int test_option_normal_value_f = test_option_normal_value<float>();
int test_option_normal_value_d = test_option_normal_value<double>();

int main()
{
	return 0;
}