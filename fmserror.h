// fmserror.h - error handling using NaNs
// Use the macros FLOAT_ERROR_SET("message") to set NaN with pointer to "message"
// and FLOAT_ERROR_GET(x) to get message pointer from NaN
#pragma once
#include <limits>
#include <cstdint>

namespace fms {

	// return NaN with message pointer in significand
	inline double float_error_set(const char* s)
	{
		union {
			uint64_t u;
			double d;
		} ud = { .d = std::numeric_limits<double>::quiet_NaN() };
		static_assert(sizeof(ud) == sizeof(double));

		static uint64_t m = (UINT64_MAX >> 12); // significand mask
		ud.u = (ud.u & ~m) | (reinterpret_cast<uint64_t>(s) & m);

		return ud.d;
	}
	// extract message pointer from NaN
	inline const char* float_error_get(double x)
	{
		union {
			uint64_t u;
			double d;
		} ud = { .d = x };

		static uint64_t m = (UINT64_MAX >> 12); // significand mask
		ud.u &= m;

		return reinterpret_cast<const char*>(ud.u);
	}
}