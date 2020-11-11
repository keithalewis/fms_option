// fms_error.h - error handling using NaNs
// Use the macro FLOAT_ENSURE to return a NaN with pointer to "message"
#pragma once
#include <limits>
#include <cstdint>

// git config --get remote.origin.url
// git rev-parse HEAD

#ifndef FLOAT_ENSURE
#define ENSURE_HASH_(x) #x
#define ENSURE_STRZ_(x) ENSURE_HASH_(x)
#define ENSURE_FILE "file: " __FILE__
#define ENSURE_LINE "\nline: " ENSURE_STRZ_(__LINE__)
#ifdef __FUNCTION__
#define ENSURE_FUNC "\nfunction: " __FUNCTION__
#else
#define ENSURE_FUNC ""
#endif
#define ENSURE_SPOT ENSURE_FILE ENSURE_LINE ENSURE_FUNC
#define ENSURE_MSG(e) ENSURE_SPOT "\nensure: \"" #e "\" failed"
#define FLOAT_ENSURE(e) if (!(e)) {                 \
	static const char* msg = ENSURE_MSG(e);         \
	return fms::float_error_set(msg); } else (void)0
#endif

namespace fms {

	// IEEE 64-bit significand mask
	inline const uint64_t FLOAT_SIG_MASK = (UINT64_MAX >> 12);
	union UD { uint64_t u; double d; };
	static_assert(sizeof(UD) == sizeof(double));


	// return NaN with message pointer in significand
	inline double float_error_set(const char* s)
	{
		UD ud = { .d = std::numeric_limits<double>::quiet_NaN() };
		ud.u = (ud.u & ~FLOAT_SIG_MASK) | (reinterpret_cast<uint64_t>(s) & FLOAT_SIG_MASK);

		return ud.d;
	}
	// extract message pointer from NaN
	inline const char* float_error_get(double x)
	{
		UD ud = { .d = x };

		return reinterpret_cast<const char*>(ud.u & FLOAT_SIG_MASK);
	}
}