// fmserror.h - error handling using NaNs
// Use the macro FLOAT_ENSURE to return a NaN with pointer to "message"
#pragma once
#include <limits>
#include <cstdint>

#ifndef FLOAT_ENSURE
#define ENSURE_HASH_(x) #x
#define ENSURE_STRZ_(x) ENSURE_HASH_(x)
#define ENSURE_FILE "file: " __FILE__
#ifdef __FUNCTION__
#define ENSURE_FUNC "\nfunction: " __FUNCTION__
#else
#define ENSURE_FUNC ""
#endif
#define ENSURE_LINE "\nline: " ENSURE_STRZ_(__LINE__)
#define ENSURE_SPOT ENSURE_FILE ENSURE_LINE ENSURE_FUNC
#define FLOAT_ENSURE(e) if (!(e)) { return float_error_set(ENSURE_SPOT "\nensure: \"" #e "\" failed"); }
#endif

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

		return reinterpret_cast<const char*>(ud.u & m);
	}
}