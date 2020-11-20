// fms_option_payoff.h - standard option payoffs
#pragma once
#include <concepts>

namespace fms::payoff {

	// base class for standard option payoffs
	template<class K = double>
		requires std::is_floating_point_v<K>
	struct option { K strike; };

	template<class K = double>
	struct call : public option<K> {
		call(K k) : option<K>{ k } { }
	};

	template<class K = double>
	struct put : public option<K> { 
		put(K k) : option<K>{ k } { }
	};

	template<class K = double>
	struct digital_call : public option<K> { 
		digital_call(K k) : option<K>{ k } { }
	};

	template<class K = double>
	struct digital_put : public option<K> { 
		digital_put(K k) : option<K>{ k } { }
	};

}
