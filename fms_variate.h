// fms_variate.h - Concept for random variates.
// A random variable is determined by its cumulative distribution function
// F(x) = P(X <= x). Its cumulant is k(s) = log E[exp(s X)] and its
// Esscher transformed distribution is dF_s(x) = exp(s x - k(s)) dF(x) 
#pragma once
#include <concepts>

namespace fms {

	template<typename M, class X = typename M::xtype, class S = typename M::stype>
	concept variate_concept = requires (M m, X x, S s, size_t n) {
		std::semiregular<M>;
		{ m.cumulant(s, n) } -> std::convertible_to<S>;
		{ m.cdf(x, s, n) } -> std::convertible_to<X>;
	};

	template<variate_concept M>
	struct variate_model : public M {
		using M::M;
	};

	//template<class X = double, class S = X>
	//using XXX = variate_model<XXX_impl>;

}
