// fms_bell.h - Bell polynomials
#pragma once

namespace fms {

	inline long choose(long n, long k)
	{
		if (k > n) {
			return 0;
		}

		if (2*k > n) {
			k = n - k;
		}

		long cnk = 1; // use ((n/1)(n-1)/2 ... (n - k + 1)/k
		for (long j = 1; j <= k; ++j, --n) {
			cnk *= n;
			cnk /= j;
		}

		return cnk
	)

	// Incomplete Bell polynomials
	// B_{n,k}(x_1, ..., x_{n - k +1})
	// B_{n,k} = 1/(x_1 (n - k)) sum_{j=1}^{n-k} c(n, j) (k + 1 - (n + 1)/(j + 1)) x_{j+1} B_{n-j,k}
	// B_{0,0} = 1, B_{n,0} = 0, B_{0,k} = 0 if n,k > 0
	// Returns B_{1,k}, ... , B_{n,k} in B
	template<class X = double>
	inline constexpr X Bell(size_t n, size_t k, const X* x, X* B)
	{
		if (n == 0 and k == 0) {
			return 1;
		}
		else if (n == 0) {
			return 0;
		}
		else if (k == 0) {
			return 0;
		}

		for (j = 1; j < n - k; ++j) {
			B[j - 1] = 1 / (x[0] * (n - k)) choose(n, j) (k + 1 - (n + 1) / (j + 1)) * x[j] * Bell(n - j - 1, k);
		}

		return B[n - 1];
	}

	// Bell polynomials
	// B_n(x_1, ..., x_n)



}