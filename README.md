# General Option Pricing and Greeks

This library implements general European option (forward) pricing and greeks.
The underlying payoff _F_ is parameterized by _forward_ _f_ and _vol_ _s_
via _F_ = _f e<sup>sX - κ(s)</sup_ where _κ(s)_ = log _E_[exp(_s X_)]
is the _cumulant_ of _X_. We can, and do, assume _X_ has mean 0 and variance 1
so _E_[_F_] = _f_ and Var(log(_F_)) = _s_<sup>2</sup>.
See [Option Pricing](https://keithalewis.github.io/math/op.html) for details.
 
To implement a model of the variate _X_
write a (value type) class with member functions `cumulant(S s, size_t n)` and
`cdf(X x, S s, size_t n)` 
that implement the derivatives of the cumulant of _X_ and the derivatives of the cumulative distribution
function of the _Esscher transform_ _X<sub>s</sub>_.

If _X_ is normal then _κ(s)_ = _s<sup>2</sup>_/2 and _X<sub>s</sub>_ = _X_ + _s_.
See [normal_variate.h](https://github.com/keithalewis/fmsoption/blob/master/fms_variate_normal.h)
for the implementation.

European value and greeks of puts and calls can be calculated using the `option` class.
```C++
normal<> N;
option o(N);
o.value(f, s, k); // value of call with forward s, vol s, and strike k
o.value(f, s, call(k)) // same
o.delta(f, s, put(k)); // delta of put with forward s, vol s, and strike k
```

Implied vol is calculated using `option::implied(f, s, k)`
