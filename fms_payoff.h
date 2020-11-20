// fms_option_payoff.h - standard option payoffs
#pragma once

namespace fms::payoff {

	// NVI base class for standard option payoffs
	template<class K = double>
	class payoff_base {
		K k;
	public:
		typedef typename K type;

		payoff_base(K k = 0)
			: k(k)
		{ }
		payoff_base(const payoff_base&) = delete;
		payoff_base& operator=(const payoff_base&) = delete;
		virtual ~payoff_base()
		{ }
		K strike() const
		{
			return k;
		}
	};

	template<class K = double>
	class call : public payoff_base<K> { };

	template<class K = double>
	class put : public payoff_base<K> { };

	template<class K = double>
	class digital_call : public payoff_base<K> { };

	template<class K = double>
	class digital_put : public payoff_base<K> { };

}
