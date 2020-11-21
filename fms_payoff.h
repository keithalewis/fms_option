// fms_option_payoff.h - standard option payoffs
#pragma once
#include <concepts>

namespace fms::payoff {

	// base class for standard option payoffs
	template<template<class> class O, class K = double>
		requires std::is_floating_point_v<K>
	struct option {
		typedef O<K> type;
		K strike; 
	};

	template<class K = double>
	struct put : public option<put, K> {
		typedef put<K> type;
		put(K k) : option<put, K>{ k } { }
	};

	template<class K = double>
	struct call : public option<call, K > {
		typedef call<K> type;
		call(K k) : option<call, K>{ k } { }
	};

	/*
	template<class K = double>
	struct digital_call : public option<K> { 
		digital_call(K k) : option<K>{ k } { }
	};

	template<class K = double>
	struct digital_put : public option<K> { 
		digital_put(K k) : option<K>{ k } { }
	};
	*/
}
