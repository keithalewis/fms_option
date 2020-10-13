# General Option Pricing and Greeks

The forward value of a _put option_ with strike _k_ is  _p_ = E[max{_k_ - _F_, 0}],
where _F_ is the price of the underlying at expiration. Every positive random variable can be
parameterized by _F_ = _f_ _e_<sup>_s X_ - _κ_(_s_)</sup> where _κ(s)_ = log E[_e_<sup>_s X_</sup>] is the cumulant of _X_.
We can assume _X_ has mean 0 and variance 1 by adjusting _f_ and _s_
appropriately.
Note _f_ is the expected value of _F_ and _s_<sup>2</sup> is the variance of  log _F_. 

We have _p_ = E[max{_k_ - _F_, 0}] = E[(_k_ - _F_)1(_F_ &le; _k_)] = _k_ _P_(_F_ &le; _k_) - _f_ _P<sub>s</sub>_(_F_ &le; _k_),
where _dP<sub>s</sub>_/_dP_ = _e_<sup>_s X_ - _κ_(_s)_</sup> is the Esscher transform of _P_.
The formula for the option value is

&emsp;_p_ = _k_ _P_(_X_ &le; _z_) - _f_ _P<sub>s</sub>_(_X_ &le; _z_),

where _z_ = (log(_k_/_f_) + _κ_(_s_))/_s_ is the _moneyness_.

The forward value of a _call option_ with strike _k_ is  _c_ = E[max{_F_ - _k_, 0}].
Since max{_F_ - _k_, 0} - max{_k_ - _F_, 0} = _F_ - _k_ we
have 

&emsp;_c_ = _p_ + _f_ - _k_ = _f_ _P<sub>s</sub>_(_X_ &gt; _z_) - _k_ _P_(_X_ &gt; _z_).

## Delta
The _delta_ of a put option is _dp_/_df_ = E[-_e_<sup>_s X_ - _κ_(_s_)</sup> 1(_F_ &le; _k_)] = -_P<sub>s</sub>_(_X_ &le; _z_). The delta of
a call option is _dc_/_df_ = _dp_/_df_ + 1 = _P<sub>s</sub>_(_X_ &gt; _z_).

## Gamma
The _gamma_ of a put option is _d_<sup>2</sup>p/_df_<sup>2</sup>.
