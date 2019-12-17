#pragma once
#include <iostream>
#include <Eigen/Dense>


///@file
///@brief Contains the class for Z/N coefficients

///////////////////////////////////
///The class of Z/N coefficients

///The operators are self explanatory
/////////////////////////////////
template<int N>
class Z {
public:
	static const int order = N; ///<The N of Z/N
	int x; ///< A modulo N number
	Z<N>() : x(0) {} ///<Default value 0
	Z<N>(int x) : x(x% N) {} ///<Take residue mod N
	Z<N>(long x) : x(x% N) {} ///<Take residue mod N
	Z<N>(float x) : Z<N>(static_cast<int>(x)) {}  ///<Take residue mod N
	inline explicit operator char() const { 
		return static_cast<char>(x);
	}
	inline explicit operator short() const {
		return static_cast<short>(x);
	}
	inline explicit operator int() const {
		return static_cast<int>(x);
	}
	inline explicit operator long() const {
		return static_cast<long>(x);
	}
	inline explicit operator float() const {
		return static_cast<float>(x);
	}
	inline explicit operator double() const {
		return static_cast<double>(x);
	}
	inline Z<N> operator +(const Z<N>& b) const {
		return Z<N>(x + b.x);
	}
	inline Z<N> operator -(const Z<N>& b) const {
		return Z<N>(x - b.x);
	}
	inline Z<N> operator +=(const Z<N>& b) {
		x = (x + b.x) % N;
		return *this;
	}
	inline Z<N> operator -=(const Z<N>& b) {
		x = (x - b.x) % N;
		return *this;
	}
	inline Z<N> operator -() const {
		return Z<N>(-x);
	}
	inline Z<N> operator %(const Z<N>& b) const {
		return Z<N>(x % b.x);
	}
	inline bool operator ==(const Z<N>& b) const {
		return (x == b.x);
	}
	inline bool operator !=(const Z<N>& b) const {
		return (x != b.x);
	}
	inline bool operator <(const Z<N>& b) const {
		return (x < b.x);
	}
	inline bool operator <=(const Z<N>& b) const {
		return (x <= b.x);
	}
};

template<int N>
inline Z<N> operator -(const Z<N>& a) {
	return Z<N>(-a.x);
}
template<int N>
inline Z<N> operator *(const Z<N>& a, const Z<N>& b) { //Eigen needs this to be non member
	return Z<N>(a.x * b.x);
}
template<int N>
inline Z<N> operator /(const Z<N>& a, const Z<N>& b) {
	return Z<N>(a.x / b.x);
}
template<int N>
inline Z<N> abs(const Z<N>& a) {
	return Z<N>(abs(a.x));
}
template<int N>
inline std::ostream& operator<<(std::ostream& out, const Z<N>& a) { //Eigen needs this to be non member
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
	Z(int x) :x(x % 2) {} ///<Take residue mod N
	Z(long x) :x(x%2){}
	Z(float x) : Z<2>(static_cast<int>(x)) {} ///<Take residue mod N
	inline explicit operator char() const {
		return static_cast<char>(x); 
	}
	inline explicit operator short() const {
		return static_cast<short>(x);
	}
	inline explicit operator int() const {
		return static_cast<int>(x);
	}
	inline explicit operator long() const {
		return static_cast<long>(x);
	}
	inline explicit operator float() const {
		return static_cast<float>(x);
	}
	inline explicit operator double() const {
		return static_cast<double>(x);
	}
	inline Z<2> operator +(const Z<2>& b) const {
		return Z<2>(x ^ b.x);
	}
	inline Z<2> operator -(const Z<2>& b) const {
		return Z<2>(x ^ b.x);
	}
	inline Z<2> operator +=(const Z<2>& b) {
		x = x ^ b.x;
		return *this;
	}
	inline Z<2> operator -=(const Z<2>& b) {
		x = x ^ b.x;
		return *this;
	}
	inline Z<2> operator -() {
		return *this;
	}
	inline Z<2> operator %(const Z<2>& b) const {
		return Z<2>(x % b.x);
	}
	inline bool operator ==(const Z<2>& b) const {
		return (x == b.x);
	}
	inline bool operator !=(const Z<2>& b) const {
		return (x != b.x);
	}
	inline bool operator <(const Z<2>& b) const {
		return (x < b.x);
	}
	inline bool operator <=(const Z<2>& b) const {
		return (x <= b.x);
	}
};

inline Z<2> operator -(const Z<2>& a) {
	return a;
}

inline Z<2> operator *(const Z<2>& a, const Z<2>& b) { //Eigen needs this to be non member
	return Z<2>(a.x * b.x);
}
inline Z<2> operator /(const Z<2>& a, const Z<2>& b) {
	return Z<2>(a.x / b.x);
}
inline Z<2> abs(const Z<2>& a) {
	return a;
}
inline std::ostream& operator<<(std::ostream& out, const Z<2>& a) {
	return out << a.x;
}

///See the Eigen documentation for this. 
namespace Eigen {
	///Specializing NumTraits to Z/nZ coefficients
	template<int N>
	struct NumTraits<Z<N>>
	{
		typedef Z<N> Real;
		typedef Z<N> NonInteger;
		typedef Z<N> Nested;
		typedef int Literal;
		enum {
			IsComplex = 0,
			IsInteger = 0,
			IsSigned = 0,
			RequireInitialization = 0,
			ReadCost = 1,
			AddCost = 3,
			MulCost = 3
		};
		static inline Z<N> digits10() { return 0; }
		static inline Z<N> dummy_precision() { return 0; }
	};
}


//We use SFINAE to detect if a given type is Z<N>
namespace {
	template <typename T, typename = void>
	struct SFINAE : std::false_type {	};

	template <typename T>
	struct SFINAE<T, decltype(T::order, void())> : std::true_type {	};

	template <class T>
	constexpr bool detector(std::false_type) {
		return 0;
	}

	template <class T>
	constexpr bool detector(std::true_type) {
		return 1;
	}
}

namespace Mackey {
	///Detects if T=Z<N>
	template <class T>
	constexpr bool is_finite_cyclic() {
		return detector<T>(SFINAE<T>{});
	}
}
