# General Option Pricing and Greeks

This library implements general European option pricing and greeks.
The underlying payoff _F_ is parameterized by _forward_ _f_ and _vol_ _s_
via _F_ = _f e<sup>sX - κ(s)</sup_ where _κ(s)_ = log _E_[exp(_s X_)]
 is the _cumulant of _X_.
