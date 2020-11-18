// fms_variate.h - Interface class for random variates.
#pragma once
#include <concepts>

namespace fms {

	// NVI base class for variates
	template<class X = double, class S = double>
	struct variate_base {
		typedef typename X type;
		typedef typename S ctype;

		variate_base()
		{ }
		variate_base(const variate_base&) = delete;
		variate_base& operator=(const variate_base&) = delete;
		virtual ~variate_base()
		{ }

		X pdf(X x, S s = 0) const
		{
			return cdf_(x, s, 1);
		}
		// transformed cumulative distribution function and derivatives
		X cdf(X x, S s = 0, size_t n = 0) const
		{
			return cdf_(x, s, n);
		}
		// log E[exp(sX)]
		S cumulant(S s, size_t n = 0) const
		{
			return cumulant_(s, n);
		}
	private:
		virtual X cdf_(X x, S s, size_t n) const = 0;
		virtual S cumulant_(S s, size_t n) const = 0;
	};
	
	// implement for a specific model and make a copy
	template<class M, class X = M::type, class S = M::ctype>
	class variate_model : public variate_base<X, S> {
		M m;
	public:
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
