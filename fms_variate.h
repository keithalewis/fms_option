// fms_variate.h - Interface class for random variates.
// A random variable is determined by its cumulative distribution function
// F(x) = P(X <= x). Its cumulant is k(s) = log E[exp(s X)] and its
// Esscher transform is dP^s = exp(s X - k(s)) dP with
// E^s[g(X)] = E[g(X) exp(s X - k(s))]
#pragma once
#include <concepts>

namespace fms {

	// NVI base class for variates
	template<class X = double, class S = double>
	struct variate_base {
		typedef typename X xtype;
		typedef typename S stype;

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
	template<class M, class X = M::xtype, class S = M::stype>
		requires std::semiregular<M>
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

	// Model with mean 0 variance 1
	template<class M, class X = M::xtype, class S = M::stype>
		requires std::semiregular<M>
	class variate_standard : public variate_base<X, S> {
		M m;
		X mu, sigma;
	public:
		variate_standard(const M& m)
			: m(m), mu(m.cumulant(0,1)), sigma(::sqrt(m.cumulant(0,2)))
		{ }
		variate_standard(const variate_standard&) = default;
		variate_standard& operator=(const variate_standard&) = default;
		~variate_standard()
		{ }

		X cdf_(X x, S s = 0, size_t n = 0) const override
		{
			return m.cdf(mu + sigma * x, s, n) * ::pow(sigma, X(n));
		}
		S cumulant_(S s, size_t n = 0) const override
		{
			return m.cumulant(s / sigma, n) / ::pow(sigma, X(n)) - (n == 0 ? s * mu / sigma : n == 1 ? mu / sigma : X(0));
		}
	};

	/*
	template<typename M, typename X = M::xtype>
	concept variate_model = requires (M m, X x, X s, size_t n) {
		typename M::xtype;
		// n-th derivative of cumulative distribution function
		{ m.cdf(x, s, n) } -> std::convertible_to<X>;
		// n-th derivative of cumulant
		{ m.cumulant(s, n) } -> std::convertible_to<X>;
	};
	*/
}
