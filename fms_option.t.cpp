// fms_option.t.cpp - test fms option
#include <cassert>
#include <iostream>
#include "fms_normal.h"
#include "fms_option.h"

using namespace fms::option;
using namespace fms::variate;

template<class X>
int test_option()
{
	X eps = std::numeric_limits<X>::epsilon();
	X f = X(100);
	X s = X(0.1); // 3-month 20% vol
	X k = X(100);

	X x;
	x = moneyness<normal<X>>(f, s, k);
	x -= X(0.05);
	ensure (fabs(x) < eps);
	x = put_value<normal<X>>(f, s, k);
	x -= X(3.9877611676744920);
	ensure(fabs(x) <= 10 * eps);

	fms::option_model<normal<X>,X,X,X> N;
	x = N.moneyness(f, s, k);
	x -= X(0.05);
	ensure(fabs(x) < eps);
	x = N.put_value(f, s, k);
	x -= X(3.9877611676744920);
	ensure(fabs(x) <= 10 * eps);

	x = x;

	return 0;
}

int main()
{
	try {
		test_option<float>();
		test_option<double>();
	}
	catch (const std::exception& ex) {
		std::cerr << ex.what() << std::endl;
	}

	return 0;
}