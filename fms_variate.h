// fms_variate.h - Variate concept
#pragma once
#include <concepts>

namespace fms {

	// NVI base class for variate models
	template<class X = double, class S = double>
	struct variate_base {
		virtual ~variate_base()
		{ }
		X pdf(X x, S s = 0) const
		{
			return cdf_(x, s, 1);
		}
		X cdf(X x, S s = 0, size_t n = 0) const
		{
			return cdf_(x, s, n);
		}
		S cumulant(S s, size_t n = 0) const
		{
			return cumulant_(s, n);
		}
	private:
		virtual X cdf_(X x, S s, size_t n) const = 0;
		virtual S cumulant_(S s, size_t n) const = 0;
	};
	
	template<class M, class X = M::type, class S = M::ctype>
	class variate_model : public variate_base<X, S> {
		M m;
	public:
		typedef X type;
		typedef S ctype;
		variate_model(const M& m)
			: m(m)
		{ }
		variate_model(const variate_model&) = default;
		variate_model& operator=(const variate_model&) = default;
		~variate_model()
		{ }
		X cdf_(X x, S s = 0, size_t n = 0) const override
		{
			return m.cdf(x, s, n);
		}
		S cumulant_(S s, size_t n = 0) const override
		{
			return m.cumulant(s, n);
		}
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
