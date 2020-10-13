// fmsoption.t.cpp - test fms option
#include <cassert>
#include "fmsnormal.h"
#include "fmsoption.h"

using namespace fms;

int test_float_error()
{
	char msg[] = "message";
	double x = float_error_set(msg);
	assert(std::isnan(x));
	const char* s = float_error_get(x);
	assert(0 == strcmp(s, msg));
	assert(s == msg);

	return 0;
}

int test_option()
{
	double f = 100;
	double s = 0.2;
	double k = 100;

	double p;
	p = moneyness<normal>(f, s, k);
	p = put_value<normal>(f, s, k);
	f = -1;
	p = moneyness<normal>(f, s, k);
	if (std::isnan(p)) {
		const char* s = float_error_get(p);
		s = s;
	}

	p = p;

	return 0;
}

int main()
{
	test_float_error();
	test_option();

	return 0;
}