// fms_sequence.h - iterators with operator bool() const
#pragma once
#include <functional>
#include <iterator>
#include <type_traits>

namespace fms::sequence {

	struct end {};

	template<class I, class T = std::iterator_traits<I>::value_type>
	class base {
		I i;
	public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = T;
		base(const I& i)
			: i(i)
		{ }
		base(const base&) = default;
		base& operator=(const base&) = default;
		virtual ~base()
		{ }
		operator==(const base&) const
		{
			return !operator bool();
		}
		operator bool() const
		{
			return true;
		}
		value_type operator*() const
		{
			return *i;
		}
		base& operator++()
		{
			++i;

			return *this;
		}
	};

	// right fold using binop
	template<class Binop, class S, class T = S::value_type>
	class fold : public base<S,T> {
		T t;
	public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = T;
		fold(const S& s, value_type t0)
			: base<S, T>(s), t(t0)
		{
			if (s) {
				t = Binop{}(t0, S::operator*());
			}
		}
		fold(const fold&) = default;
		fold& operator=(const fold&) = default;
		virtual ~fold()
		{ }
		using S::operator bool() const;
		T operator*() const override
		{
			return *t;
		}
		fold& operator++() override
		{
			S::operator bool() && t = Binop{}(t, S::operator*());
			S::operator++();

			return *this;
		}
	};
	template<class S>
	inline auto sum(S s)
	{
		return fold<std::plus<S::value_type>>(s, S::value_type(0));
	}
	template<class S>
	inline auto product(S s)
	{
		return fold<std::multiplies<S::value_type>>(s, S::value_type(1));
	}

	template<class S>
	inline S::value_type back(S s)
	{
		S::value_type t;

		while (s) {
			t = *s;
			++s;
		}

		return t;
	}
	template<class S>
	inline auto sum(S s)
	{
		return fold<std::plus<S::value_type>>(s, 0);
	}
}

/*
int s[] = {1,2,3}
auto f = fold<std::plus<int>, ...>(sequence(i), 0);
*f == 1;
++f;
*f == 3;
++f;
*f == 6;
++f;
!f;


*/