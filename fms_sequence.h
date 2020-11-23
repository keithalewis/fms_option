// fms_sequence.h - iterators with operator bool() const
#pragma once
#include <concepts>
#include <functional>
#include <iterator>
#include <type_traits>


template<class S>
concept sequence = requires (const S s) {
		{ s.operator bool() } -> std::same_as<bool>;
	} && std::forward_iterator<S>;


namespace fms::sequence {

	struct end {};

	// NVI interface
	template<class I, class T = std::iterator_traits<I>::value_type>
	struct base {
		typedef typename std::forward_iterator_tag iterator_category;
		typedef typename T value_type;
		typedef typename I type;

		base()
		{ }
		base(const base&) = default;
		base& operator=(const base&) = default;
		virtual ~base()
		{ }
		/*
		operator==(const base&) const
		{
			return !operator bool(); // ???
		}
		*/
		operator bool() const
		{
			return op_bool();
		}
		T operator*() const
		{
			return op_star();
		}
		base& operator++()
		{
			return op_incr();
		}
	private:
		virtual bool op_bool() const = 0;
		virtual value_type op_star() const = 0;
		virtual base& op_incr() = 0;
	};
	
	// apply function to sequence
	template<class F, class I, 
		class T = I::value_type,
		class U = std::invoke_result_t<F, T>>
	class apply : public base<I, U>
	{
		F f;
		I i;
	public:
		apply(const F& _f, const I& _i)
			: f(_f), i(_i)
		{ }
		apply(const apply&) = default;
		apply& operator=(const apply&) = default;
		~apply()
		{ }

		bool op_bool() const override
		{
			return i;
		}
		U op_star() const override
		{
			return f(*i);
		}
		apply& op_incr() override
		{
			++i;

			return *this;
		}
	};

	/*
	// filter sequence based on predicte
	template<class P, class I, class T = I::value_type>
	class filter : public base<I, T>
	{
		P p;
		I i;
	public:
		filter(const P& _p, const I& _i)
			: p(_p), i(_i)
		{
			while (i and !p(*i)) {
				++i;
			}
		}
		filter(const filter&) = default;
		filter& operator=(const filter&) = default;
		~filter()
		{ }

		bool op_bool() const override
		{
			return i;
		}
		T op_star() const override
		{
			return *i;
		}
		filter& op_incr() override
		{
			++i;
			while (i and !p(*i)) {
				++i;
			}

			return *this;
		}
	};
	*/

	// mask sequence based on predicte
	template<class M, class I, class T = I::value_type>
	class mask : public base<I, T>
	{
		M m;
		I i;
		void next()
		{
			while (i and !*m) {
				++m;
				++i;
			}
		}
	public:
		mask(const M& _m, const I& _i)
			: m(_m), i(_i)
		{
			next();
		}
		mask(const mask&) = default;
		mask& operator=(const mask&) = default;
		~mask()
		{ }

		bool op_bool() const override
		{
			return i;
		}
		T op_star() const override
		{
			return *i;
		}
		mask& op_incr() override
		{
			++m;
			++i;
			next();

			return *this;
		}
	};

	// right fold using binop
	template<class Op, class I, class T = I::value_type>
	class fold : public base<I,T> {
		Op op;
		I i;
		T t;
	public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = T;

		fold()
		{ }
		fold(const Op& op, const I& _i, value_type t0)
			: op(op), i(_i), t(t0)
		{
			if (i) {
				t = op(t0, *i);
			}
		}
		fold(const fold&) = default;
		fold& operator=(const fold&) = default;
		~fold()
		{ }

		bool op_bool() const override
		{
			return i;
		}
		T op_star() const override
		{
			return t;
		}
		fold& op_incr() override
		{
			if (++i) {
				t = op(t, *i);
			}

			return *this;
		}
	};
	template<class S>
	inline auto sum(S s)
	{
		return fold(std::plus<S::value_type>{}, s, S::value_type(0));
	}

	template<class S>
	inline auto product(S s)
	{
		return fold(std::multiplies<S::value_type>, s, S::value_type(1));
	}

	template<class S>
	inline auto back(S s)
	{
		typename S::value_type t(0);

		while (s) {
			t = *s;
			++s;
		}

		return t;
	}

	template<class I, class T = typename std::iterator_traits<I>::value_type>
	class counted : public base<I, T> {
		I i;
		size_t n;
	public:
		counted()
			: n(0)
		{ }
		counted(I i, size_t n)
			: i(i), n(n)
		{ }
		counted(const counted&) = default;
		counted& operator=(const counted&) = default;
		~counted()
		{ }

		// remaining size
		size_t size() const
		{
			return n;
		}

		bool op_bool() const override
		{
			return n != 0;
		}
		T op_star() const override
		{
			return *i;
		}
		counted& op_incr() override
		{
			if (n) {
				--n;
				++i;
			}

			return *this;
		}
	};
	template<class T>
	struct array : public counted<const T*, T> {
		template<size_t N>
		array(const T(&i)[N])
			: counted<const T*,T>(i, N)
		{ }
	};
}

template<class S, class T = S::value_type>
	requires std::is_base_of_v<fms::sequence::base<typename S::type, T>, S>
inline auto operator==(S s, T t)
{
	return fms::sequence::apply([s, t](T si) { return si == t; }, s);
}
template<class S, class T = S::value_type>
requires std::is_base_of_v<fms::sequence::base<typename S::type, T>, S>
inline auto operator!=(S s, T t)
{
	return fms::sequence::apply([s, t](T si) { return si != t; }, s);
}
template<class S, class T = S::value_type>
requires std::is_base_of_v<fms::sequence::base<typename S::type, T>, S>
inline auto operator<(S s, T t)
{
	return fms::sequence::apply([s, t](T si) { return si < t; }, s);
}
template<class S, class T = S::value_type>
requires std::is_base_of_v<fms::sequence::base<typename S::type, T>, S>
inline auto operator<=(S s, T t)
{
	return fms::sequence::apply([s, t](T si) { return si <= t; }, s);
}
template<class S, class T = S::value_type>
requires std::is_base_of_v<fms::sequence::base<typename S::type, T>, S>
inline auto operator>(S s, T t)
{
	return fms::sequence::apply([s, t](T si) { return si > t; }, s);
}
template<class S, class T = S::value_type>
requires std::is_base_of_v<fms::sequence::base<typename S::type, T>, S>
inline auto operator>=(S s, T t)
{
	return fms::sequence::apply([s,t](T si) { return si >= t; }, s);
}

template<class M, class S>
//requires std::is_base_of_v<fms::sequence::base<typename M::type, T>, S>
inline auto operator|(S s, M m)
{
	return fms::sequence::mask(m, s);
}
