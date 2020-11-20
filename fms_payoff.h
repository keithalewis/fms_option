// fms_option_payoff.h - standard option payoffs
#pragma once

namespace fms::payoff {

	// NVI base class for standard option payoffs
	template<class X = double>
	struct payoff_base {
		typedef typename X type;

		payoff_base()
		{ }
		payoff_base(const payoff_base&) = delete;
		payoff_base& operator=(const payoff_base&) = delete;
		virtual ~payoff_base()
		{ }
		X strike() const
		{
			return strike_();
		}
	private:
		virtual X strike_() const = 0;
	};

	template<class X = double>
	class call : public payoff_base<X> {
		X k;
	public:
		call(X k)
			: k(k)
		{ }
		X strike_() const override
		{
			return k;
		}
	};

	template<class X = double>
	class put : public payoff_base<X> {
		X k;
	public:
		put(X k)
			: k(k)
		{ }
		X strike_() const override
		{
			return k;
		}
	};

	template<class X = double>
	class digital_call : public payoff_base<X> {
		X k;
	public:
		digital_call(X k)
			: k(k)
		{ }
		X strike_() const override
		{
			return k;
		}
	};

	template<class X = double>
	class digital_put : public payoff_base<X> {
		X k;
	public:
		digital_put(X k)
			: k(k)
		{ }
		X strike_() const override
		{
			return k;
		}
	};
}
