# General option pricing and greeks

The forward value of a put option with strike $k$ is $p = E[\max\{k - F, 0\}]$,
where $F$ is underlying at expiration. Every positive random variable can be
parameterized by $F = f\exp(s X - κ(s))$ where $X$ has mean 0 and variance 1, 
and $κ(s) = \log E[\exp(s X)] is the cumulant of $X$.
Note $f$ is the expected value of $F$ and $s^2$ is the variance of $\log F$. 

We have $p = E[max\{k - F, 0\}] = E[(k - F)1(F \le k)] = k P(F\le k) - f P_s(F \le k)$
where $dP_s/dP = \exp(s X - κ(s))$ is the Esscher transform.
The formula for the option value is
$$
	p = k P(X \le z) - f P_s(X \le z),
$$
where $z = (\log(k/f) + κ(s))/s$ is the _moneyness_.

## Delta
The _delta_ of the option is $dp/df = E[-exp(s X - κ(s)) 1(F\le k)] = -P_s(z)$.