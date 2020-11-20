// fms_option_payoff.h - standard option payoffs
#pragma once
#include <concepts>

namespace fms::payoff {

	// base class for standard option payoffs
	template<class K = double>
	//requires = std::is_floating_point_v<K>
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
	struct call : public payoff_base<K> {
		using payoff_base<K>::payoff_base;
	};

	template<class K = double>
	struct put : public payoff_base<K> { 
		using payoff_base<K>::payoff_base;
	};

	template<class K = double>
	struct digital_call : public payoff_base<K> { 
		using payoff_base<K>::payoff_base;
	};

	template<class K = double>
	struct digital_put : public payoff_base<K> { 
		using payoff_base<K>::payoff_base;
	};

}
