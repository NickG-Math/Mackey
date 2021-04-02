#pragma once
#include "../Z_n.hpp"

///@file
///@brief Contains the class for Z/N coefficients

namespace mackey
{

	namespace implementation_details
	{

		template <int64_t N, typename T>
		std::array<std::array<Z_mod<N, T>, N>, N> get_add_table()
		{
			std::array<std::array<T, N>, N> adder;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					adder[i][j].x = (i + j) % N;
			return adder;
		}

		template <int64_t N, typename T>
		std::array<std::array<Z_mod<N, T>, N>, N> get_sub_table()
		{
			std::array<std::array<T, N>, N> suber;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					suber[i][j].x = (i + N - j) % N;
			return suber;
		}

		template <int64_t N, typename T>
		std::array<std::array<Z_mod<N, T>, N>, N> get_mul_table()
		{
			std::array<std::array<Z_mod<N, T>, N>, N> muler;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					muler[i][j].x = (i * j) % N;
			return muler;
		}

		template <int64_t N, typename T>
		std::array<std::array<Z_mod<N, T>, N>, N> get_div_table()
		{
			std::array<std::array<Z_mod<N, T>, N>, N> diver;
			for (int j = 1; j < N; j++)
			{
				int k = 1;
				for (; (j * k) % N != 1; k++)
				{
				}
				for (int i = 0; i < N; i++)
				{
					diver[i][j].x = (i * k) % N;
				}
			}
			return diver;
		}
		template <int64_t N, typename T>
		const static std::array<std::array<Z_mod<N, T>, N>, N> add_table = get_add_table<N, T>();
		template <int64_t N, typename T>
		const static std::array<std::array<Z_mod<N, T>, N>, N> sub_table = get_sub_table<N, T>();
		template <int64_t N, typename T>
		const static std::array<std::array<Z_mod<N, T>, N>, N> mul_table = get_mul_table<N, T>();
		template <int64_t N, typename T>
		const static std::array<std::array<Z_mod<N, T>, N>, N> div_table = get_div_table<N, T>();

	}

	template <int64_t N, typename T>
	Z_mod<N, T>::Z_mod(int x) : x(abs(x) % N) {}

	template <int64_t N, typename T>
	Z_mod<N, T>::Z_mod(int64_t x) : x(abs(x) % N) {}

	template <int64_t N, typename T>
	Z_mod<N, T>::operator char() const
	{ //used by basis element/ kernel mod image
		return static_cast<char>(x);
	}

	template <int64_t N, typename T>
	Z_mod<N, T>::operator short() const
	{ //used by basis element/ kernel mod image
		return static_cast<short>(x);
	}

	template <int64_t N, typename T>
	Z_mod<N, T>::operator int() const
	{ //used by basis element/ kernel mod image
		return static_cast<int>(x);
	}

	template <int64_t N, typename T>
	Z_mod<N, T>::operator int64_t() const
	{ //used by basis element/ kernel mod image
		return static_cast<int64_t>(x);
	}

	template <int64_t N, typename T>
	Z_mod<N, T>::operator unsigned char() const
	{ //used by basis element/ kernel mod image
		return static_cast<unsigned char>(x);
	}

	template <int64_t N, typename T>
	Z_mod<N, T>::operator unsigned short() const
	{ //used by basis element/ kernel mod image
		return static_cast<unsigned short>(x);
	}

	template <int64_t N, typename T>
	Z_mod<N, T>::operator unsigned int() const
	{ //used by basis element/ kernel mod image
		return static_cast<unsigned int>(x);
	}

	template <int64_t N, typename T>
	Z_mod<N, T>::operator uint64_t() const
	{ //used by basis element/ kernel mod image
		return static_cast<uint64_t>(x);
	}

	template <int64_t N, typename T>
	Z_mod<N, T> Z_mod<N, T>::operator+(Z_mod<N, T> b) const
	{
		if constexpr (std::is_same_v<T, bool>)
			return Z_mod<N, T>(x != b.x);
		else
			return implementation_details::add_table<N, T>[x][b.x];
	}

	template <int64_t N, typename T>
	Z_mod<N, T> Z_mod<N, T>::operator-(Z_mod<N, T> b) const
	{
		if constexpr (std::is_same_v<T, bool>)
			return Z_mod<N, T>(x != b.x);
		else
			return implementation_details::sub_table<N, T>[x][b.x];
	}

	template <int64_t N, typename T>
	Z_mod<N, T> &Z_mod<N, T>::operator+=(Z_mod<N, T> b)
	{
		*this = *this + b;
		return *this;
	}

	template <int64_t N, typename T>
	Z_mod<N, T> &Z_mod<N, T>::operator-=(Z_mod<N, T> b)
	{
		*this = *this - b;
		return *this;
	}

	template <int64_t N, typename T>
	Z_mod<N, T> &Z_mod<N, T>::operator*=(Z_mod<N, T> b)
	{
		*this = *this * b;
		return *this;
	}

	template <int64_t N, typename T>
	Z_mod<N, T> &Z_mod<N, T>::operator/=(Z_mod<N, T> b)
	{
		*this = *this / b;
		return *this;
	}

	template <int64_t N, typename T>
	bool Z_mod<N, T>::operator==(Z_mod<N, T> a) const
	{
		return x == a.x;
	}

	template <int64_t N, typename T>
	bool Z_mod<N, T>::operator!=(Z_mod<N, T> a) const
	{
		return x != a.x;
	}

	template <int64_t N, typename T>
	bool Z_mod<N, T>::operator<=(Z_mod<N, T> a) const
	{
		return x <= a.x;
	}

	template <int64_t N, typename T>
	Z_mod<N, T> operator-(Z_mod<N, T> a)
	{
		if constexpr (std::is_same_v<T, bool>)
			return a;
		else
			return implementation_details::sub_table<N, T>[0][a.x];
	}

	template <int64_t N, typename T>
	Z_mod<N, T> operator*(Z_mod<N, T> a, Z_mod<N, T> b)
	{ //Eigen needs this to be non member
		if constexpr (std::is_same_v<T, bool>)
			return Z_mod<N, T>(a.x && b.x);
		else
			return implementation_details::mul_table<N, T>[a.x][b.x];
	}

	//N better be prime for this
	template <int64_t N, typename T>
	Z_mod<N, T> operator/(Z_mod<N, T> a, Z_mod<N, T> b)
	{
		if constexpr (std::is_same_v<T, bool>)
		{
			if (b.x)
				return a;
			else
			{
				std::cerr << "Division by 0:\n"
						  << a << " / " << b << "\n\n";
				abort();
			}
		}
		else
			return implementation_details::div_table<N, T>[a.x][b.x];
	}

	template <typename T>
	T abs(T a)
	{
		static_assert(std::is_integral_v<T> || SFINAE::is_finite_cyclic<T>::value, "Only implemented for integer and Z/N coefficients");
		if constexpr (std::is_integral_v<T>)
			return std::abs(a);
		else
			return a;
	}

	template <int64_t N, typename T>
	std::ostream &operator<<(std::ostream &out, const Z_mod<N, T> a)
	{ //Eigen needs this to be non member
		return out << a.x;
	}

}