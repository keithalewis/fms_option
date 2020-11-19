// fms_sequence.t.cpp - test sequences
#include <cassert>
#include "fms_sequence.h"

using namespace fms::sequence;

int test_counted()
{
	int i[] = { 1,2,3 };
	{
		counted c(i, 3);
		counted c2{ c };
		c = c2;
		assert(c);
		assert(*c == i[0]);
		++c;
		assert(c);
		assert(*c == i[1]);
		++c;
		assert(c);
		assert(*c == i[2]);
		++c;
		assert(!c);
	}
	{
		counted c(i, 3);
		auto f = [](int i) -> double { return 0.1 * i; };
		auto ac = apply(f, c);
		assert(ac);
		assert(*ac == 0.1*i[0]);
		++ac;
		assert(ac);
		assert(*ac == 0.1*i[1]);
		++ac;
		assert(ac);
		assert(*ac == 0.1*i[2]);
		++ac;
		assert(!ac);
	}
	{
		counted c(i, 3);
		auto pc = mask(c >= 2, c);
		assert(pc);
		assert(*pc == 2);
		++pc;
		assert(pc);
		assert(*pc == 3);
		++pc;
		assert(!pc);
	}
	{
		counted c(i, 3);
		auto pc = c | c >= 2;
		assert(pc);
		assert(*pc == 2);
		++pc;
		assert(pc);
		assert(*pc == 3);
		++pc;
		assert(!pc);
	}

	{
		counted c(i, 3);
		auto fc = fold(std::plus<int>{}, c, 0);
		assert(fc);
		assert(*fc == 1);
		++fc;
		assert(fc);
		assert(*fc == 3);
		++fc;
		assert(*fc == 6);
		++fc;
		assert(!fc);
	}
	{
		counted c(i, 3);
		auto fc = sum(c);
		assert(fc);
		assert(*fc == 1);
		++fc;
		assert(fc);
		assert(*fc == 3);
		++fc;
		assert(*fc == 6);
		++fc;
		assert(!fc);
	}
	{
		counted c(i, 3);
		auto fc = sum(c);
		assert(back(fc) == 6);
	}
	{
		counted c(i, 3);
		auto fc = sum(c | c > 1);
		assert(back(fc) == 5);
	}
	{
		counted c(i, 3);
		auto fc = sum(c | c <= 2);
		assert(back(fc) == 3);
	}
	{
		counted c(i, 3);
		auto fc = sum(c | c == 2);
		assert(back(fc) == 2);
	}
	{
		array c(i);
		auto fc = sum(c | c != 2);
		assert(back(fc) == 4);
	}
	{
		array c(i);
		assert(c);
		assert(c.size() == 3);
		++c;
		assert(c);
		assert(c.size() == 2);
	}

	return 0;
}
int test_counted_ = test_counted();