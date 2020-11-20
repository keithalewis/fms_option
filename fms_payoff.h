// fms_option_payoff.h - standard option payoffs
#pragma once
#include <concepts>

namespace fms::payoff {

	// base class for standard option payoffs
	template<class K = double>
	//requires = std::is_floating_point_v<K>
	struct payoff_base {
		K strike;
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
