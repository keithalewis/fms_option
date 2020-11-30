// fms_variate_handle.h - Interface class for random variates.
#pragma once
#include "fms_variate.h"

namespace fms {

	// NVI base class for variates
	template<class X = double, class S = double>
	struct variate_base_impl {
		typedef typename X xtype;
		typedef typename S stype;

		variate_base_impl()
		{ }
		variate_base_impl(const variate_base_impl&) = delete;
		variate_base_impl& operator=(const variate_base_impl&) = delete;
		virtual ~variate_base_impl()
		{ }

		X pdf(X x, S s = 0, size_t n = 0) const
		{
			return cdf_(x, s, n - 1);
		}
		// transformed cumulative distribution function and derivatives
		X cdf(X x, S s = 0, size_t n = 0) const
		{
			return cdf_(x, s, n);
		}
		// (d/ds)^n log E[exp(sX)]
		S cumulant(S s, size_t n = 0) const
		{
			return cumulant_(s, n);
		}
	private:
		virtual X cdf_(X x, S s, size_t n) const = 0;
		virtual S cumulant_(S s, size_t n) const = 0;
	};

	template<class X = double, class S = X>
	using variate_base = variate_model<variate_base_impl<X,S>>;

	// implement for a specific model and make a copy
	template<class M, class X = M::xtype, class S = M::stype>
	class variate_handle : public variate_base<X, S>
	{
		M m;
	public:
		variate_handle(const M& m)
			: m(m)
		{ }
		variate_handle(const variate_handle&) = default;
		variate_handle& operator=(const variate_handle&) = default;
		~variate_handle()
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
	class variate_standard : public variate_base<X, S> 
	{
		M m;
		X mu, sigma;
	public:
		variate_standard(const M& m)
			: m(m), mu(m.cumulant(0, 1)), sigma(::sqrt(m.cumulant(0, 2)))
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

}
#pragma once
