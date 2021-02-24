#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <array>


///@file
///@brief Contains the class for Z/N coefficients

///////////////////////////////////
///The class of Z/N coefficients

///The operators are self explanatory. N should be prime for division to work.
/////////////////////////////////
template<char N>
class Z {
public:
	static const char order = N; ///<The N of Z/N
	char x; ///< A modulo N number
	Z<N>() : x(0) {} ///<Default value 0
	Z<N>(int x) : x(abs(x)% N) {} ///<Take residue mod N
	Z<N>(float x) : Z<N>(static_cast<int>(x)) {}  ///<Take residue mod N
	explicit operator char() const { 
		return static_cast<char>(x);
	}
	explicit operator short() const {
		return static_cast<short>(x);
	}
	explicit operator int() const {
		return static_cast<int>(x);
	}
	explicit operator long() const {
		return static_cast<long>(x);
	}
	explicit operator float() const {
		return static_cast<float>(x);
	}
	explicit operator double() const {
		return static_cast<double>(x);
	}
	Z<N> operator +(Z<N> b) const {
		return add_table[x][b.x];
	}
	Z<N>& operator ++() {
		if (x == N - 1)
			x = 0;
		else
			x++;
		return *this;
	}
	Z<N> operator ++(int) {
		if (x == N - 1)
			x = 0;
		else
			x++;
		return *this;
	}
	Z<N> operator -(Z<N> b) const {
		return sub_table[x][b.x];
	}
	Z<N> operator +=(Z<N> b) {
		x = add_table[x][b.x].x;
		return *this;
	}
	Z<N> operator -=(Z<N> b) {
		x = sub_table[x][b.x].x;
		return *this;
	}
	Z<N> operator %(Z<N> b) const {
		return Z<N>(x % b.x);
	}
	bool operator ==(Z<N> b) const {
		return (x == b.x);
	}
	bool operator !=(Z<N> b) const {
		return (x != b.x);
	}
	bool operator <(Z<N>& b) const {
		return (x < b.x);
	}
	bool operator <=(Z<N>& b) const {
		return (x <= b.x);
	}

	const static std::array<std::array<Z<N>, N>, N> add_table;
	const static std::array<std::array<Z<N>, N>, N> sub_table;
	const static std::array<std::array<Z<N>, N>, N> mul_table;
	const static std::array<std::array<Z<N>, N>, N> div_table;

	static std::array<std::array<Z<N>, N>, N> make_addition_table() {
		std::array<std::array<Z<N>, N>, N> adder;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				adder[i][j].x = (i + j) % N;
		return adder;
	}

	static std::array<std::array<Z<N>, N>, N> make_subtraction_table() {
		std::array<std::array<Z<N>, N>, N> subber;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				subber[i][j].x = (i+ N- j) % N;
		return subber;
	}

	static std::array<std::array<Z<N>, N>, N> make_mult_table() {
		std::array<std::array<Z<N>, N>, N> multer;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				multer[i][j].x = (i * j) % N;
		return multer;
	}

	static std::array<std::array<Z<N>, N>, N> make_div_table() {
		std::array<std::array<Z<N>, N>, N> diver;
		for (int j = 1; j < N; j++) {
			int k = 1;
			for (; (j * k) % N != 1;k++) {}
			for (int i = 0; i < N; i++) {
				diver[i][j].x = (i * k) % N;
			}
		}
		return diver;
	}

};



template<char N>
const std::array<std::array<Z<N>, N>, N> Z<N>::add_table = Z<N>::make_addition_table();
template<char N>
const std::array<std::array<Z<N>, N>, N> Z<N>::sub_table = Z<N>::make_subtraction_table();
template<char N>
const std::array<std::array<Z<N>, N>, N> Z<N>::mul_table = Z<N>::make_mult_table();
template<char N>
const std::array<std::array<Z<N>, N>, N> Z<N>::div_table = Z<N>::make_div_table();

template<char N>
Z<N> operator -(Z<N> a) {
	return Z<N>::sub_table[0][a.x];
}

template<char N>
Z<N> operator *(Z<N> a, Z<N> b) { //Eigen needs this to be non member
	return Z<N>::mul_table[a.x][b.x];
}

//N better be prime for this
template<char N>
Z<N> operator /(Z<N> a, Z<N> b) {
	return Z<N>::div_table[a.x][b.x];
}
template<char N>
Z<N> abs(Z<N> a) {
	return a;
}
template<char N>
std::ostream& operator<<(std::ostream& out, const Z<N> a) { //Eigen needs this to be non member
	return out << a.x;
}


///////////////////////////////////
///The class of Z/2 coefficients

///Specialization of the Z/N template class for better performance. Operators are self explanatory
/////////////////////////////////
template<>
class Z<2> {
public:
	static const int order = 2; ///<The 2 of Z/2
	bool x; ///<A modulo 2 value is Boolean
	Z() :x(0) {} ///<Initialize by 0
	Z(bool x) :x(x) {}
	Z(int x) :x(x % 2) {} ///<Take residue mod 2
	Z(long x) :x(x % 2) {}
	Z(float x) : Z<2>(static_cast<int>(x)) {} ///<Take residue mod 2
	explicit operator char() const {
		return static_cast<char>(x);
	}
	explicit operator short() const {
		return static_cast<short>(x);
	}
	explicit operator int() const {
		return static_cast<int>(x);
	}
	explicit operator long() const {
		return static_cast<long>(x);
	}
	explicit operator float() const {
		return static_cast<float>(x);
	}
	explicit operator double() const {
		return static_cast<double>(x);
	}
	Z<2>& operator ++() {
		x = !x;
		return *this;
	}
	Z<2> operator ++(int) {
		x = !x;
		return *this;
	}
	Z<2> operator +(Z<2> b) const {
		return Z<2>(x != b.x);
	}
	Z<2> operator -(Z<2> b) const {
		return *this+b;
	}
	Z<2> operator +=(Z<2> b) {
		x = (x != b.x);
		return *this;
	}
	Z<2> operator -=(Z<2> b) {
		x = (x != b.x);
		return *this;
	}
	Z<2> operator %(Z<2> b) const { //the only way for a%b=1 is b=0 and a=1
		return Z<2>(x && !b.x);
	}
	bool operator ==(Z<2> b) const {
		return (x == b.x);
	}
	bool operator !=(Z<2> b) const {
		return (x != b.x);
	}
	bool operator <(Z<2> b) const {
		return (x < b.x);
	}
	bool operator >(Z<2> b) const {
		return (x > b.x);
	}
	bool operator <=(Z<2> b) const {
		return (x <= b.x);
	}
	bool operator >=(Z<2> b) const {
		return (x >= b.x);
	}
	Z<2> operator /(Z<2> b) const {
		if (b.x)
			return *this;
		else
			throw("Division by 0 in Z/2 coefficients");
	}
};

Z<2> operator -(Z<2> a) { //-v for vector v needs this to be non member
	return a;
}

Z<2> operator *(Z<2> a, Z<2> b) { //Eigen needs this to be non member
	return Z<2>(a.x && b.x);
}
 Z<2> abs(Z<2> a) {
	return a;
}
 std::ostream& operator<<(std::ostream& out, Z<2> a) {
	return out << a.x;
}



///See the Eigen documentation for this. 
namespace Eigen {
	///Specializing NumTraits to Z/nZ coefficients
	template<char N>
	struct NumTraits<Z<N>>
	{
		typedef Z<N> Real;
		typedef Z<N> Nested;
		typedef int Literal;
		enum {
			IsComplex = 0,
			IsInteger = 0,
			IsSigned = 0,
			RequireInitialization = 0,
			ReadCost = 1,
			AddCost = 1,
			MulCost = 1
		};
		static inline int dummy_precision() { return 0; }
		static inline int digits10() { return 0; }
	};

	///Specializing NumTraits to Z/2 coefficients
	template<>
	struct NumTraits<Z<2>>
	{
		typedef Z<2> Real;
		typedef Z<2> Nested;
		typedef int Literal;
		enum {
			IsComplex = 0,
			IsInteger = 0,
			IsSigned = 0,
			RequireInitialization = 0,
			ReadCost = 1,
			AddCost = 1,
			MulCost = 1
		};
		static inline int dummy_precision() { return 0; }
		static inline int digits10() { return 0; }
	};

}