// fms_variate.h - Variate concept
#pragma once
#include <concepts>

namespace fms {

	// base class for variate models
	template<class X = double>
	struct variate_base {
		virtual ~variate_base()
		{ }
	};

	/*
	template<typename M, typename X = M::type>
	concept variate_model = requires (M m, X x, X s, size_t n) {
		typename M::type;
		// n-th derivative of cumulative distribution function
		{ m.cdf(x, s, n) } -> std::convertible_to<X>;
		// n-th derivative of cumulant
		{ m.cumulant(s, n) } -> std::convertible_to<X>;
	};
	*/
}
